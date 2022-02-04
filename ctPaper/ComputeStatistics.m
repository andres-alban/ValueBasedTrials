function [ tfig3, samplepathsfig3, out, outQmax, trialSampleSize, inmbStop ] = ComputeStatistics( basic, mat, samplePathProf, tProfIncrease  ) 
% Calculate operating characteristics from the bootstrap

% Establish the length of the t vector (containing the discrete possible
% stopping times, measured in pairwise allocatons) and the length of the
% output matrix (one row for each run of the bootstrap)
lent = length( tProfIncrease )  ; 
lens = size( samplePathProf, 1 ) ; 
% Approximate the upper and lower stopping boundaries with a vector the same 
% dimension as the sample paths for the bootstrap
boundaryProfUpper = zeros( 1,  lent ) ; 
boundaryProfLower = zeros( 1,  lent ) ; 
% Select the discrete points on the boundary which will be compared with 
% the sample paths from the bootstrap ;
for i = 1:lent 
    t = tProfIncrease( i ) ;
    diff = abs( mat.tvec - t ) ;  
    minDiff = min( diff ) ; 
    indminDiff = minDiff == diff ;
    boundaryProfUpper( i ) =  max( indminDiff .* mat.bndupper )  ;
    boundaryProfLower( i ) =  min( indminDiff .* mat.bndlower ) ;
end
% Obtain the stopping times ;
% Vectors to hold the index for the time that the lower/upper boundaries are
% hit
crossLowerIndex = zeros( lens, 1 ) ; 
crossUpperIndex = zeros( lens, 1 ) ; 
% Vectors to hold the time of hitting the lower and upper boundary
crossLowerTime = zeros( lens, 1 ) ; 
crossUpperTime = zeros( lens, 1 ) ;
% = 1 if the lower/upper boundary is hit
crossLower = samplePathProf < boundaryProfLower ;
crossUpper = samplePathProf > boundaryProfUpper ;
% accumulate the number of hits of lower and upper boundary
cumCrossLower = cumsum( crossLower, 2 ) ; 
cumCrossUpper = cumsum( crossUpper, 2 ) ; 
% for each run of bootstrap, if never crosses lower/upper boundary,
% allocate final element of crossLower/crossUpper = 1 (otherwise next step
% does not work)
for j = 1:1:lens
    if ( cumCrossLower( j, lent ) == 0 ) 
        crossLower( j, lent ) = 1 ; 
    end
    if ( cumCrossUpper( j, lent ) == 0 ) 
        crossUpper( j, lent ) = 1 ;     
    end
end
% Establish index and crossing times of lower/upper boundaries
for j = 1:1:lens
    crossLowerIndex( j )  = find( crossLower( j, : ) > 0, 1 ) ; 
    crossLowerTime( j ) = tProfIncrease( crossLowerIndex( j ) ) ; 
    crossUpperIndex( j )  = find( crossUpper( j, : ) > 0, 1 ) ;
    crossUpperTime( j ) = tProfIncrease( crossUpperIndex( j ) ) ; 
end
% Establish whether crosses the lower or upper boundary first 
crossLowerFirst = crossLowerIndex <= crossUpperIndex ; 
crossUpperFirst = 1 - crossLowerFirst ; 
% Establish the stopping time in terms of the number of pairwise
% allocations
stoppingTime = ( crossLowerFirst .* crossLowerTime ) + ( crossUpperFirst .* crossUpperTime ) ;
% Replace paths that run to above TMax with the stopping time TMax
for j = 1:1:lens
    if( stoppingTime( j ) > basic.TMax + basic.t0 )
        stoppingTime( j ) = basic.TMax + basic.t0 ;
    end
end
% Finish time of trial is equal to stopping time plus number of pairwise
% allocations in pipeline
finishTime = stoppingTime + ceil( basic.tau ) ; 
% This is the actual number of pairwise allocations made in the sequential
% trail
trialSampleSize = finishTime - ceil( basic.tau ) - ceil( basic.t0 ) ; 
samplePathProf1 = samplePathProf .* ( tProfIncrease <= finishTime ) ;
% Saving 
saving = ( basic.numpairs - trialSampleSize ) .* basic.cResearch  ; 
% Plot some sample paths
tProf1 = tProfIncrease .* ones( lens, lent ) .* ( tProfIncrease <= finishTime ) ;  
% If stop first along lower boundary, calculate the stopping time for the
% block
crossLowerTimeFinal = crossLowerTime .* crossLowerIndex ;
crossLowerTimeFinal = crossLowerTimeFinal > 0 ;
% If stop first along lower boundary, calculate the stopping time for the
% block
crossUpperTimeFinal = crossUpperTime .* crossUpperIndex ;
crossUpperTimeFinal = crossUpperTimeFinal > 0 ;
% samplePathProf identifies the posterior mean by blocks of 10 patient
% pairs, but the final point in the sample path must be interpolated so as
% to take account of the true delay

inmbStop = zeros( lens, 1 ) ;
for j = 1:1:lens
    indicator = 0 ;
    for i = 2:1:lent
        if( samplePathProf1( j, i ) == 0 )
            cut = ( mod( basic.tau, basic.incr ) ) / basic.incr ;
            samplePathProf1( j, i ) = samplePathProf( j, i - 1 ) + ( samplePathProf( j, i ) - samplePathProf( j, i - 1 ) )  * cut ;
            inmbStop( j ) = samplePathProf1( j, i ) ;
            tProf1( j, i ) = tProfIncrease( i - 1 ) + ( tProfIncrease( i ) - tProfIncrease( i - 1 ) ) * ( mod( basic.tau, basic.incr ) ) / basic.incr ;
            indicator = 1 ; 
        else
            inmbStop( j ) = samplePathProf1( j, i ) ;
        end
        if( indicator == 1 )
            break
        end 
    end
end
indi = tProf1 ~= 0 ;
finalDecisionLower = inmbStop < 0 ; 
finalDecisionUpper = 1 - finalDecisionLower ; 
changeMindUpper2Lower = crossUpperFirst .* finalDecisionLower ; 
changeMindLower2Upper = crossLowerFirst .* ( 1 - finalDecisionLower ) ; 
%%%%%%%%%%%%%%%%%%%%%%OUTPUTS FOR BOOTSTRAPPED PATHS%%%%%%%%%%%%%%%%%%%%%%%
% Out vector: see below for definitions
out = horzcat( trialSampleSize, saving, inmbStop, finalDecisionLower, finalDecisionUpper, crossLowerFirst, crossUpperFirst, changeMindUpper2Lower, changeMindLower2Upper,  crossLowerTime, crossUpperTime ) ; 
% Output some sample paths for figure 3
samplePathProf2a = samplePathProf1( 1, : ) ; 
samplePathProf2b = samplePathProf1( 9, : ) ; 
samplePathProf2c = samplePathProf1( 45, : ) ;
tProf2 = tProf1 .* ( tProf1 <= finishTime ) ; 
tProf2a = tProf2( 1, : ) ; 
tProf2b = tProf2( 9, : ) ; 
tProf2c = tProf2( 45, : ) ; 
tfig3 = [tProf2a; tProf2b; tProf2c];
samplepathsfig3 = [samplePathProf2a;samplePathProf2b;samplePathProf2c];


%%%%%%%%%%%%%%%%%OUTPUTS FOR PATHS WHICH RUN TO QMAX%%%%%%%%%%%%%%%%%%%%%%%
inmbQMax = samplePathProf( :, lent ) ;
QMaxDecisionUpper = inmbQMax > 0 ;
QMaxDecisionLower = inmbQMax < 0 ;
outQmax = horzcat( inmbQMax, QMaxDecisionLower, QMaxDecisionUpper ) ; 
