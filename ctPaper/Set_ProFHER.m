function[ basic, advanced] = Set_ProFHER(TMax,tryMC)
% SetProFHER: the input parameters for the ProFHER Application paper with
% adjusted parameters to account for mix of treatments at the start of the
% trial. Adapted from the JRSSB paper of Chick, Forster, Pertile
%(alpha order) and extended for the ProFHER application paper (available as
% a discussion paper here:
% http://ideas.repec.org/p/yor/yorken/19-01.html )
%
% (c) 2014, 2018 S Chick, P Pertile, M Forster, A Alban
% Created: 17 April 2014
% Last touched: 11 April 2018
% Updated for ProFHER parameter values and to account for a 
% mix of treatments: 10 April 2018 
%
% For sources of parameter values consult Table 1 and Appendix A 
% of Application paper. 
% 
%%%%%% Case study parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[basic, advanced] = DelayInputConstructor() ;
basic.lambda = 20000 ; % WTP for one unit of effectiveness 
basic.pN = 0.39 ; % proportion using New at the start of the study
%basic.pN = 0.0000 ; % proportion using New at the start of the study
basic.P = 42000 ; % number of patients to benefit from the adoption decision
basic.sigma = 4400 ; % sampling sd for INMB (in £) BASE CASE
basic.t0 = 2 ; % effective sample size of the prior distribution ( must be greater than 1 )
delayinyears = 12 /12 ; % delay in observing INMB in the study ( in years )   
basic.numpairs = 125 ; % number of pairwise allocations in the study ( NB this is NOT TMax )
accrualperiod =  32 / 12 ; % actual duration of the recruitment period of trial ( in years )             
basic.ICostN = 0 ; % cost of switching all patients to New ( in £ )
basic.ICostS = 0 ; % cost of switching all patients to Standard ( in £ )
basic.cResearch = 4080 ; % research cost per pairwise allocation ( in £ )
annualdiscountrate = 0 ; % set to zero if no discounting required
patientsperyear = basic.numpairs / accrualperiod  ;  % maximum number of patients put into trial per year
basic.tau = delayinyears * patientsperyear ; % delay as measured by the number of pairwise allocations
%basic.numpairs = 100 ; % R&R check for clinical trials (blank if you want to run the proper analysis for 125 and 250)
basic.TMax = TMax ; % assumption - the maximum number of pairwise allocations that can be made in the trial
basic.tryMC = tryMC ; % set to 1 if you would like to run the Monte Carlo simulation instead of the bootstrap
basic.Compare125250 = 1 ; % set to 1 if you would like to compare the TMax = 250 and TMax = 125 operating characteristics
basic.linewid = 1 ; % set linewidth for plots from bootstrap analysis
basic.fixedCost = 0 ; % estimated fixed costs of trial
%basic.mumax = 10 * 50 * 250000 / (basic.t0 + basic.TMax) ; % 500000 was the max s_Y in Martin-Paolo code for the same application
basic.mumax = 20000 ; 
basic.mumin = -1 * basic.mumax ; 
basic.mu0 = 0; % prior mean 
basic.z  = 0 ; % critical value for classical statistical inference test ( e.g. 1.96 )
basic.theta = exp( - log(1+annualdiscountrate) / patientsperyear ); %  per patient discount rate
%%%%%%%%%
basic.online = false;  % set to true if per patient reward is gained ( online learning )
% MIX OF TREATMENT ADJUSTMENTS. Adjust terms in objective function to fit
% special cases from Alban, Chick, Forster (2018), "Extending a Bayesian
% decision-theoretic approach to value-based sequential clinical trial
% design" (Proceedings of the 2018 Winter Simulation Conference). Hereafter
% the "WSC paper". Note that when there is no online learning, it is not
% necessary to divide through by (1-2pN) and this is accounted for below 
% for the case of 0 < basic.pN < 1/2. 
if basic.pN == 0                % sec 3.2.1 of WSC paper: current practice is standard treatment
   basic.factor = 1.0 ;        % adjustment factor in denominator
   basic.c = basic.cResearch ; % adjusted research cost per pairwise allocation
   basic.PPatients = basic.P ; % adjusted patients to benefit
   basic.ICost = basic.ICostN ;% adjusted switching cost
   basic.addon = 0 ;           % adjusted add-on at end of trial
elseif basic.pN < 1/2                       % sec 3.2.2 of WSC paper: mix of treatments is not an option
   basic.factor = basic.online * ( 1 - 2 * basic.pN ) + ( 1 - basic.online ) ;   % adjustment factor in denominator ( refer to WSC paper ) 
   basic.c = basic.cResearch / basic.factor ;    % adjusted research cost per pairwise allocation
   basic.PPatients = basic.P / basic.factor ;          % adjusted patients to benefit
   basic.ICost = ( basic.ICostN - basic.ICostS ) / basic.factor  ; % adjusted switching cost
elseif basic.pN == 1/2
   basic.factor = 1.0 ; % adjustment factor in denominator ( refer to WSC paper ) 
   basic.c = basic.cResearch ; % adjusted research cost per pairwise allocation
   basic.PPatients = basic.P ; % adjusted patients to benefit
   basic.ICost = basic.ICostN - basic.ICostS  ; % adjusted switching cost
else
   warning('WARNING: pN should be between 0 and 1/2, or both I_N and I_S should be 0');
   basic.factor = 1.0 ; % adjustment factor in denominator ( refer to Theory paper ) 
   basic.c = basic.cResearch ; % adjusted research cost per pairwise allocation
   basic.PPatients = basic.P ; % adjusted patients to benefit
   basic.ICost = 0.0 ;  % adjusted switching cost
   basic.addon = 0.0 ;  % adjusted add-on at end of trial
end
advanced.fixedP = true; % set to true if there are exactly PPatients to follow once a treatment or alternative is selected
% set to false if the patients which are not included in trial due
% to early stopping are also included in the treatment arm, that
% is, if fixedP is false, then there are P + (TMAX - T) patients
% which are included, where T is the stopping time (number of
% patients) included in the trial (and T \in 0, 1, ..., TMAX).
advanced.nochangeallowed = false;    % set to false (the better default) meaning that it is ok
% to change what you think is best after the added tau data points
% come in. set this to true if you want to ignore the added tau
% pending data points to make a decision.

%%%%%% Numerical computation parameters %%%%%%

%advanced.MinGridPerStdev = (basic.TMax-basic.t0)/s_step_t;    % if this is negative, code tries to use its negative as the delta-t
%advanced.MinGridPerStdev = 200; 
advanced.MinGridPerStdev = 20 ;    % minimum number of grid points per stdev of output, or negative to use dt=1
% For plotting
advanced.smoothed = true ;   % true or false if curve smoothing is desired.
advanced.NumGridsForContours = 100; % number of time values for generating contour plots (also num plots in animation to use during recursion over t)
advanced.bigfontsize=16 ;
advanced.smallfontsize=14 ;
advanced.fontname='Helvetica';
advanced.titlestring='Surgery versus sling';
advanced.filestring='ProFHER';
advanced.dirstring='ProFHER';
advanced.PLOTSIMS = false ;
advanced.saveplot = true;           % set to true to save plots in files, false to not save files automatically
advanced.DOPLOT= false ;
advanced.keepAllOutput = true ;
advanced.MinGridPerStdev = 30; %200;
advanced.numinsimFreqDeltaVec=200;  % 200 probably enough, or 350 is more than enough
advanced.simNumReps = 10;  % set to specific number of replications required for this specific experiment. 
end

