%% Initialization parameters
% This section needs to be run regardless of the figure that you want to
% generate. The Other sections can be run in any order.
dosave = true; % If you want to save figures as eps files. Set to false otherwise
basic = ProfherSetParameters(); % Parameters for the ProFHER application

%% Setup cost function (Figure 3)
bestfit_ab = ProfherSetup_cost(dosave);

%% Base parameters maximization (Figure 4)
ProfherBaseCase(dosave)

%% Marginal value vs. marginal cost (Figure 5)
% make sure the output of ProfherSetup_cost is correctly used inside the
% following script
ProfherWTPrecrate_power(dosave)

%% Benefits of optimizing recruitment rate with fixed r (Figure 6 left panel)
ProfherRecruitmentRate(dosave)

%% Benefits of optimizing recruitment rate assuming no time-dependence (Figure 6 right panel)
ProfherLossIgnoreTimeDependence(dosave)

%% Base case power and CPCS plots (Figure 7)
ProfherPower_CPCS_Tr(dosave)

%% Sequential delay and recruitment rate (Figure 8)
% warning: this takes in the order of ~20 minutes in a modern computer
Profher_sequential_delay_recruitmentrate(dosave)

%% Numerical comparative statics (Appendix D)
% warning: the following code takes many hours to run because of the
% explore_parameter_space_*.m scripts. The output of those scripts saves
% the explore_*.mat files, which are necessary before you can run the
% corresponding SimCompStats_*.m scripts. The SimCompStats_*.m save the
% output to comp_stats_*.txt files.
    %% Case I (Table EC.4)
        %% Optimal trials
        explore_parameter_space_caseI
        %% Derivative computation
        SimCompStats_caseI
    %% Case II (Table EC.5)
        %% Optimal trials
        explore_parameter_space_caseII
        %% Derivative computation
        SimCompStats_caseII
    %% Case III (Table EC.6)
        %% Optimal trials
        explore_parameter_space_caseIII
        %% Derivative computation
        SimCompStats_caseIII
    %% Case IV (Table EC.3)
        %% Optimal trials
        explore_parameter_space_caseIV
        %% Derivative computation
        SimCompStats_caseIV

%% Profher application supplement (Appendix G)
%% Post-trial population plots (Figure EC.1)
ProfherTimeHorizon_Tr_disc(dosave)

%% Discount rate (Figure EC.2)
ProfherDiscountRate(dosave)

%% Fraction in new treatment (Figure EC.3)
Profherpn_Tr(dosave)

%% Online Offline (Figure EC.4)
ProfherOnlineOffline(dosave)

%% Asymptotes as of post-adoption population cases I-III (Figures EC.5 and EC.6)
ProfherTimeHorizon(dosave)

%% Asymptotes as of post-adoption population case IV (Figure EC.7)
ProfherTimeHorizon_Tr_caseIV(dosave)
