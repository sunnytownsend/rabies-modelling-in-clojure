clear all
close all
clc
tic


%% load density grid, generate landgrid and grid for model cases
%% and set where the location of incursion will be
load BaliDensityGrid.mat %contains density grid
%densitygrid = round(unifrnd(0,10,nrows,ncols));
landgrid = (densitygrid>-1);
nrows = size(densitygrid,1);
ncols = size(densitygrid,2);
incursion_row = 89;
incursion_col = 93;
VCgrid = zeros(nrows,ncols);%

%% probability observe any given rabid dog
p_obs = 0.2; %WHO-20% for confirmed cases ??

%% parameterisation of serial interval gamma distribution
% m = 25; v = 50;A = m^2/v; B = v/m; %plot(gampdf(0:100,A,B))
shapeA = 1.463;
scaleB = 1/0.06207;% plot(gampdf(0:100,shapeA,scaleB))

%% Parameterisation of negative binomial dist for number of bites per dog
%% and set relationship between R0 and dog density
r = 1.2;%mean for distribution is based on R0
R0string = '1pt2';
k = 1.33;
p=k/(r+k);
R0density = 0;%0=no relationship

%% parameterisation of local movement spatial kernal (negative binomial dist)
shape = 0.215;%from Katie's gamma dist
rate = 0.245;
mdist=shape/rate;%convert into parameters for negative binomial distribution
vdist=shape/(rate^2);
kdist=mdist^2/(vdist-mdist);%0.28
pdist=kdist/(mdist+kdist);%0.62
%                 y=nbinrnd(kdist,pdist,1000,1);
%                 hist(y,max(y))

%% Parameterisation of rate of loss of vaccination coverage
%% N(t+1) = lambda*N(t)
%% lambda=0.9965 causes VC to drop 50% in a year
%% if vaccine_immunity_duration = 2 years, lambda = 0.9967
vaccine_immunity_duration = 2; %vd='1y';%years.
daily_waning_rate = 1/(vaccine_immunity_duration*360);%rabisin lasts 3 years
annual_dog_turnover = 0.5; %turnoverstr='pt5';%can't be 1
annual_whatsleft=1-annual_dog_turnover;
daily_whatsleft = annual_whatsleft^(1/360);
lambda = daily_whatsleft - daily_waning_rate;
%I had just got whatsleft and turnover the wrong way round

%% Set when vaccination starts, the coverage,and whether repeated in following years
vacc_cov_achieved = 0.7;%vaccination coverage acheived at time of vaccinationvacc.
chosen_scenarios = 3; %see scenarios.m for options
id = 'Test';
no_rand_sqs_unvacc = 24*10;%must be multiple of 24 because 24 blocks e.g.3 blocks worth of random unvacc squares would be 35*24 (=~same area)
num_missing_blocks=2; %ONLY WORKS IF chosen_scenarios = 12; 
months_to_erad=[];
cases_vacc_start = 2000; startstring = num2str(cases_vacc_start);
cases_stop_epidemic = cases_vacc_start*10;
num_campaigns = 3;%0=campaign happens once(not repeated), 1=campaign happens twice
vacc_duration = 6;%duration of campaign (months), must divide into 12 i.e. 1,2,3,4,6,8,12 or 24
gap_btw_campaigns = 6;%months
number_scenarios = max(chosen_scenarios);
reactive_period = 30;%no.days use cases from to determine where to vaccinate - for reactive scenario (10)
delay_reaction = 0;%no.days delay in cases use to determine where to vacc -for reactive scenario (10)
db = 2;%doggy bus probability%

%COMPLIANCE - IMPROVING PERFORMANCE OF LOW COVERAGE CAMPAIGNS
%every time revisit, get more unvaccinated dogs
compliance=0;%1=high compliance of dogs/owners - repeat campaigns in same year reach more dogs

%% Set saving and plotting outputs
saveoutputs = 0;
save_cases=0;
runname = id;
% strcat('VC',num2str(vacc_cov_achieved*100),...
%     '_gap',num2str(gap_btw_campaigns),...
%     '_camps',num2str(num_campaigns));
plot_timeseries = 1; %plot timeseries at end of run (not map)
plotBalimap_prevacc = 0;
plotBalimap = 1;
plotBalimap_start = 0;
plotBalimap_end = 0;

%% initialisation for main code
epirun = 0;
seed = 2;


%% For multiple runs
while epirun < 1,
    %     seed = seed+1;
    
    % need to reset for each scenario
    rand('state',seed)
    randn('state',seed)
    
    %% initialising for main code
    % initialise with single infected
    clear parents caselist totalcases
    parents = index_case(incursion_row,incursion_col);
    caselist = parents;
    
    % initialise time-related info
    epidemic_launch = 0;
    epidemic_failure = 0;
    today = 0;
    VCgrid = zeros(nrows,ncols);% initialise VC grid
    
    
    %% intial plot of rabies-free
    if plotBalimap_start == 1,
        plotbali_talk(landgrid,incursion_col,incursion_row,today,caselist,VCgrid)
        title(gca,strcat('seed=',num2str(seed),', msv=',num2str(today/30)))
        if saveoutputs==1, print('-djpeg','-r300',strcat(runname,'_',num2str(ceil(today/30)),'.jpeg')); end
    end
    
    %%%%%%% PHASE 1 - EPIDEMIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% main code for initial epidemic (prior to control)
    while epidemic_launch == 0 && epidemic_failure == 0,
        today = today+1;
        
        %% classify parents as rabid today or incubating
        inftodayrows = find(parents(:,6)==today);
        parentstoday = parents(inftodayrows,:);
        notinftodayrows = find(parents(:,6)~=today);
        incubating = parents(notinftodayrows,:);
        nparentstoday = size(parentstoday,1);
        alldaughterstoday = [];
        
        %% for dogs that are 'live' today - calculate infectious offspring
        %% number, movement and incubation period
        if nparentstoday > 0,
            for i = 1:nparentstoday,
                
                %% get VC
                row = parentstoday(i,2);
                col = parentstoday(i,3);
                VC = VCgrid(row,col);%VC for grid born into
                
                % get R0 if there is a R0-density relationship
                if R0density == 1,
                    dogdensity = densitygrid(row,col);
                    r = (0.4/31)*dogdensity + 1;
                    p=k/(r+k);
                end
                
                %% number of infecteds born to a rabid dog (parent today)
                parentID=parentstoday(i,1);
                npotentialdaughters = nbinrnd(k,p,1);
                ndaughters = sum(binornd(1,(1-VC),npotentialdaughters,1));
                row=find(parents(:,1)==parentID);
                parents(row,7) = ndaughters; %number of infectious daughters
                row2=find(caselist(:,1)==parentID);
                caselist(row2,7) = ndaughters; %number of infectious daughters
                
                % rabid dog case observed?
                obs = binornd(1,p_obs);
                parents(row,8) = obs;
                caselist(row2,8) = obs;
                
                %% incubation period and movement for each daughter
                if ndaughters > 0,
                    daughters = daughters_incubation_movement(ndaughters,caselist,parentID,parentstoday,i,...
                        shapeA,scaleB,densitygrid,kdist,pdist,nrows,ncols,db);
                    
                    caselist = vertcat(caselist, daughters);
                    if i == 1,
                        alldaughterstoday = daughters;
                    else alldaughterstoday = vertcat(alldaughterstoday, daughters);
                    end
                    clear daughters
                    
                end % end to 'if ndaughters>0'
            end
            
            %% construct new parents table based on incubating and new daughter
            %% infections from today
            %% determine whether eradicated rabies via vaccination,
            %% fade out as epidemic never took off
            if size(incubating,1)>0, %still got incubating dogs
                parents = vertcat(incubating,alldaughterstoday);
            elseif size(alldaughterstoday,1)>0,
                parents = alldaughterstoday;
            else epidemic_failure = 1;%initial epidemic did not take off
            end
            
        end %end of 'any parents today?'
        
        
        %% count up total and daily new cases
        nr = find(caselist(:,6)==today);%new cases today
        or = find(caselist(:,6)<today);
        dailynewcases(today,1) = size(nr,1);
        totalcases(today,1) = size(or,1) + size(nr,1);
        if totalcases(today,1) >= cases_vacc_start,
            epidemic_launch = 1;%epidemic has taken off and reached number cases for starting control
            epirun = epirun+1;
        end
        
        if plotBalimap_prevacc == 1 && rem(today,30) == 0,
            plotbali_talk(landgrid,incursion_col,incursion_row,today,caselist,VCgrid)
            title(gca,strcat('seed=',num2str(seed),', msv=',num2str(today/30)))
            if saveoutputs==1, print('-djpeg','-r300',strcat(runname,'_',num2str(ceil(today/30)),'.jpeg')); end
        end
        
        
    end % end of while loop - know whether initial_epidemic took off or did not
    
    
    
    
    %%%%%%% PHASE 2 - DISEASE CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if epidemic_launch == 1,
        %if plot_timeseries == 1, figure('windowstyle','docked'); plot((1:today)/30,totalcases); end
        
        % store
        vacc_start_day = today+1;
        caselist_start_vacc = caselist;
        totalcases_start_vacc = totalcases;
        dailynewcases_start_vacc = dailynewcases;
        parents_start_vacc = parents;
        
        
        for scenario = chosen_scenarios,
            
            % need to reset for each scenario
            rand('state',seed)
            randn('state',seed)
            
            [scenariostring] = scenarios(scenario);
            
            % reset scenario specific variables
            eradication = 0;
            epidemic_not_controlled = 0;
            vacc_months_since_start = 0;
            VCgrid = zeros(nrows,ncols);% initialise VC grid
            clear VCisland
            today = vacc_start_day-1;
            caselist = caselist_start_vacc;
            totalcases = totalcases_start_vacc;
            dailynewcases = dailynewcases_start_vacc;
            parents = parents_start_vacc;
            
            
            while eradication == 0 && epidemic_not_controlled == 0,
                today = today+1;
                
                %% classify parents as rabid today or incubating
                inftodayrows = find(parents(:,6)==today);
                parentstoday = parents(inftodayrows,:);
                notinftodayrows = find(parents(:,6)~=today);
                incubating = parents(notinftodayrows,:);
                nparentstoday = size(parentstoday,1);
                alldaughterstoday = [];
                
                %% for dogs that are 'live' today - calculate infectious offspring
                %% number, movement and incubation period
                if nparentstoday > 0,
                    for i = 1:nparentstoday,
                        
                        %% get VC (and R0 if there is a R0-density relationship)
                        if today > 1, %use initalise code above for today=1
                            row = parentstoday(i,2);
                            col = parentstoday(i,3);
                            VC = VCgrid(row,col);%VC for grid born into
                            if R0density == 1,
                                dogdensity = densitygrid(row,col);
                                r = (0.4/31)*dogdensity + 1;
                                p=k/(r+k);
                            end
                        end
                        
                        %% number of infecteds born to a rabid dog (parent today)
                        parentID=parentstoday(i,1);
                        npotentialdaughters = nbinrnd(k,p,1);
                        ndaughters = sum(binornd(1,(1-VC),npotentialdaughters,1));
                        row=find(parents(:,1)==parentID);
                        parents(row,7) = ndaughters; %number of infectious daughters
                        row2=find(caselist(:,1)==parentID);
                        caselist(row2,7) = ndaughters; %number of infectious daughters
                        
                        % rabid dog case observed?
                        obs = binornd(1,p_obs);
                        parents(row,8) = obs;
                        caselist(row2,8) = obs;
                        
                        %% incubation period and movement for each daughter
                        if ndaughters > 0,
                            
                            daughters = daughters_incubation_movement(ndaughters,caselist,parentID,parentstoday,i,...
                                shapeA,scaleB,densitygrid,kdist,pdist,nrows,ncols,db);
                            
                            caselist = vertcat(caselist, daughters);
                            if i == 1,
                                alldaughterstoday = daughters;
                            else alldaughterstoday = vertcat(alldaughterstoday, daughters);
                            end
                            clear daughters
                            
                        end % end to 'if ndaughters>0'
                    end
                    
                    %% construct new parents table based on incubating and new daughter
                    %% infections from today
                    %% determine whether eradicated rabies via vaccination,
                    %% fade out as epidemic never took off
                    if size(incubating,1)>0, %still got incubating dogs
                        parents = vertcat(incubating,alldaughterstoday);
                    elseif size(alldaughterstoday,1)>0,
                        parents = alldaughterstoday;
                    else eradication = 1;%success!
                        months_to_erad(epirun,scenario) = ceil((today-vacc_start_day)/30);
                    end
                    
                end %end of 'any parents today?'
                
                
                % keep track of number of cases
                nr = find(caselist(:,6)==today);%new cases today
                or = find(caselist(:,6)<today);
                dailynewcases(today,1) = size(nr,1);
                totalcases(today,1) = size(or,1) + size(nr,1);
                
                %% stop if epidemic out of control
                if totalcases(today,1) >= cases_stop_epidemic, ...
                        epidemic_not_controlled = 1;
                    months_to_erad(epirun,scenario) = -1;
                end
                
                %%%% VACCINATE (at monthly increments)
                if today==vacc_start_day || rem(today-vacc_start_day,30) == 0,
                    
                    % keep track of months from start of vacc and which
                    % campaign
                    vacc_months_since_start = vacc_months_since_start+1;
                    vacc_campaign = floor((vacc_months_since_start-1)/...
                        (vacc_duration+gap_btw_campaigns) +1);
                    vacc_camp_month = vacc_months_since_start -...
                        (vacc_campaign-1)*(vacc_duration+gap_btw_campaigns);
                    seed_scenario_vaccy_vaccm_camp_campm = [seed scenario...
                        floor(vacc_months_since_start/12) vacc_months_since_start vacc_campaign vacc_camp_month];
                    
                    if scenario >= 2, %scenarios based on blocks
                        if vacc_campaign <= num_campaigns ...% clause for repeat vaccination and
                                && vacc_camp_month <= vacc_duration, % clause for month in which vacc occurs (not every month there is vaccination)
                            VCgrid = vaccinatingblocks5(vacc_camp_month,vacc_campaign,scenario,VCgrid,...
                                vacc_cov_achieved,vacc_duration,caselist,today,reactive_period,...
                                delay_reaction,no_rand_sqs_unvacc,landgrid,num_missing_blocks,compliance);
                        end
                    end
                    
                    %% plot Bali map with new and old cases, and vaccination coverage
                    %% as a white transparant layer that fades as declines with time
                    if plotBalimap == 1 & epirun == 1 %
                        %                         plotbali(densitygrid,incursion_col,incursion_row,today,caselist,...
                        %                             VCgrid)
                        plotbali_talk(landgrid,incursion_col,incursion_row,today,caselist,VCgrid)
                        title(gca,strcat('seed=',num2str(seed),', msv=',num2str(ceil(today/30))))
                        if saveoutputs==1, print('-djpeg','-r300',strcat(runname,'_',num2str(ceil(today/30)),'.jpeg')); end
                        
                    end
                    
                    %% else any other day we update VCgrid with waning coverage
                else VCgrid = VCgrid * lambda; % daily waning based on yearly waning of 50%
                end
                
                %% record mean vaccination coverage for whole island
                VCislandgrid = VCgrid.*landgrid;
                VCisland(today,1) = sum(sum(VCislandgrid))/6648;%sum of VC for all squares/no.squares in grid
                
                
                
            end
            
            %% VISUALISE timeseries - for given scenario and seed
            %% calculate time series based on caselist alone
            %% plot cases times series
            %                 plot_timeseries_cases(totalcases, dailynewcases, VCisland,...
            %                     today, caselist)
            [monthlynewcases,monthlynewcases_obs] = monthlytimeseries(...
                caselist, today, plot_timeseries, VCisland);
            
            %% collect data on cases too for generating 95% CI of trajectories
            if save_cases==1,
                MNC{epirun,1}=monthlynewcases;
                MNCO{epirun,1}=monthlynewcases_obs;
            end
            
            
        end % scenario increment
        
        if saveoutputs == 1,
            save(runname,'runname','months_to_erad','chosen_scenarios','seed')
        end
    end % end of check that using launched epidemics only
    
end %seed increment

%%
%subsitute 0s for NaNs for scenarios not included so no confusion from when
%eradication not successful(-1)
%one line for each epirun, one col for each scenario
if isempty(months_to_erad)==0; months_to_erad(find(months_to_erad==0))=NaN;end

%% final plot of rabies-free
if plotBalimap_end == 1,
    plotbali_talk(landgrid,incursion_col,incursion_row,today,caselist,VCgrid)
    title(gca,strcat('seed=',num2str(seed),', msv=',num2str(today/30)))
    if saveoutputs==1, print('-djpeg','-r300',strcat(runname,'_',num2str(ceil(today/30)),'.jpeg')); end
end

if saveoutputs == 1,
    save(strcat(runname,'_all'))
end
toc/60/60
