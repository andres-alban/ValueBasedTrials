% Routine currently in use to produce plots of the paper: Forster, Brealey,
% Chick et al. (2021) Cost-effective clinical trial design: Application of 
% a Bayesian sequential model to the ProFHER pragmatic trial
% https://journals.sagepub.com/doi/10.1177/17407745211032909

close all hidden;
clear;

% Set the parameter values 
[basic125, advanced125] = Set_ProFHER(125,0);
[basic125, advanced125, rval, msgs] = DelayInputValidator( basic125, advanced125 );
if ~rval, disp(msgs), end
[basic250, advanced250] = Set_ProFHER(250,0);
[basic250, advanced250, rval, msgs] = DelayInputValidator( basic250, advanced250 );
if ~rval, disp(msgs), end
ProductionReps = 1500;
fignum = 1;

% Compute boundaries for TMax=125
[~, mat125] = DelayCurvesRecur(basic125, advanced125);
[mat125] = DelayStageOneProFHER(basic125, advanced125, mat125 ); 
toneshot = (0:basic125.TMax);
[mat125]  = DelayOptimalBayesOneStage( basic125, advanced125, toneshot, mat125 ); 

% Compute boundaries for TMax=250
[~, mat250] = DelayCurvesRecur(basic250, advanced250);
[mat250] = DelayStageOneProFHER(basic250, advanced250, mat250 ); 
toneshot = (0:basic250.TMax);
[mat250]  = DelayOptimalBayesOneStage( basic250, advanced250, toneshot, mat250 );

% We need the following code to set the seed to the state it was used for
% the results in the paper
advanced125.CRN = true;
advanced125.CRNAcrossExperiment = 37;
DelaySimOverview( fignum, basic125, advanced125, mat125 ) ;

% Bootstrap profher data and compute statistics
[ samplepaths, tpaths, profPath, tProf, basic125 ] = bootstrapProFHER(basic125); 
[ tboot, samplepathsboot, out125, outQmax125, trialSampleSize125, inmbStop125 ] = ComputeStatistics( basic125, mat125, samplepaths, tpaths ) ; 
DelaySimOverview( fignum, basic250, advanced250, mat250 ) ;
[ samplepaths, tpaths, ~, ~, basic250 ] = bootstrapProFHER(basic250); 
[ ~, ~, out250, outQmax250, trialSampleSize250, inmbStop250 ] = ComputeStatistics( basic250, mat250, samplepaths, tpaths ) ; 

% Figure 3
DoFigure3(basic125,advanced125,mat125,mat250,tProf,profPath,tboot,samplepathsboot,1);

% Figure 4a and 4b
DoFigure4ab(basic250,advanced250,trialSampleSize250,2)
DoFigure4ab(basic125,advanced125,trialSampleSize125,3)

% Figure 4c
DoFigure4c(basic125,advanced125,trialSampleSize125,trialSampleSize250,4)

% Figure 4d
DoFigure4d(basic125,advanced125,inmbStop125, inmbStop250, 5)

% Table 1 (bootstrap)
DoTable1(basic125,out125,outQmax125);
DoTable1(basic250,out250,outQmax250);

% Compute Monte Carlo statistics
basic125.tryMC = 1;
DelaySimOverview( fignum, basic125, advanced125, mat125 ) ;
[ samplepaths, tpaths, profPath, tProf, basic125 ] = bootstrapProFHER(basic125); 
[ tboot, samplepathsboot, out125, outQmax125, trialSampleSize125, inmbStop125 ] = ComputeStatistics( basic125, mat125, samplepaths, tpaths ) ; 
basic250.tryMC = 1;
DelaySimOverview( fignum, basic250, advanced250, mat250 ) ;
[ samplepaths, tpaths, ~, ~, basic250 ] = bootstrapProFHER(basic250); 
[ ~, ~, out250, outQmax250, trialSampleSize250, inmbStop250 ] = ComputeStatistics( basic250, mat250, samplepaths, tpaths ) ; 

% Table 1 (Monte Carlo)
DoTable1(basic125,out125,outQmax125);
DoTable1(basic250,out250,outQmax250);










