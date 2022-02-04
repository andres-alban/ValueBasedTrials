function [ samplePaths, xIncrease, Profpath, x, basic ] = bootstrapProFHER( basic )
% Generates the path of the posterior mean from the ProFHER trial,
% with the number of blocks of patient pairs determined by the data that is
% made available from the ProFHER data set. Generates the resampled paths
% for the bootstrap

% number of replications required for bootstrap
n = 5000 ;
% Read the table of point estimates of QALYs and treatment cost data from
% the ProFHER trial reported in the supplemental material of Forster et al.
% (2021)
M = readmatrix('data_blocks.csv');
% dimension of M
m = size(M,1);
% select matrix of n recording the number of pairwise allocations for each
% outcome and cost for each block.
Mn = M(:,[1,6:9]);
% select matrix of outcomes from which resampling will be taken for the
% mean outcome and cost for each block. NB dimension of Mo is one less than
% that of Mn because resampling uses the actual observations and not the
% prior mean.
Mo = M(:,1:5);
% Increment of patient pairs (patient pairs per block)
basic.incr = 10 ;
% vector recording number of patient pairs for plots
x = ceil( basic.t0  + basic.tau ) : basic.incr : ceil( basic.t0 + basic.tau ) + ( ( m - 1 ) * basic.incr ) ;
% vector recording the path of the posterior mean
Profpath = zeros( 1, m ) ;

% create the true sample path for the cumulative mean INMB. Each average
%(for effectiveness new, effectiveness standard, cost new, cost standard)
% is a weighted average which is then used to calculate the estimate of
% INMB for each block
for j = 1:1:m
    if j == 1
        cumEN = Mo( j, 2 ) ;
        cumES = Mo( j, 3 ) ;
        cumCN = Mo( j, 4 ) ;
        cumCS = Mo( j, 5 ) ;
        Profpath( j ) = basic.lambda * ( cumEN - cumES ) - ( cumCN - cumCS ) ;
        cumN = Mn( j, : ) ;
    else
        cumNPrev = cumN ;
        cumN = cumNPrev + Mn( j, : ) ;
        obsEN = Mo( j, 2 ) ;
        obsES = Mo( j, 3 ) ;
        obsCN = Mo( j, 4 ) ;
        obsCS = Mo( j, 5 ) ;
        nEN = Mn( j, 2 ) ;
        nES = Mn( j, 3 ) ;
        nCN = Mn( j, 4 ) ;
        nCS = Mn( j, 5 ) ;
        cumEN = ( obsEN * nEN + cumEN * cumNPrev( 1, 2 ) ) / ( nEN + cumNPrev( 1, 2 ) ) ;
        cumES = ( obsES * nES + cumES * cumNPrev( 1, 3 ) ) / ( nES + cumNPrev( 1, 3 ) ) ;
        cumCN = ( obsCN * nCN + cumCN * cumNPrev( 1, 4 ) ) / ( nCN + cumNPrev( 1, 4 ) ) ;
        cumCS = ( obsCS * nCS + cumCS * cumNPrev( 1, 5 ) ) / ( nCS + cumNPrev( 1, 5 ) ) ;
        Profpath( j ) = basic.lambda * ( cumEN - cumES ) - ( cumCN - cumCS ) ;
    end
end
% define matrix to hold resampled paths of cumulative mean for the
% bootstrap. The dimension needs to be larger than the number of points
% used for the actual path of the posterior mean, to accommodate larger TMax
mIncrease = m + round( ( basic.TMax - m .* basic.incr ) / basic.incr ) * (round( ( basic.TMax - m .* basic.incr ) / basic.incr ) > 0) + 1 ;
% matrices to hold the differences in E and the differences in C
cumDiffE = zeros( n, mIncrease ) ;
cumDiffC = zeros( n, mIncrease ) ;
% matrix to hold the sample path for the replications
samplePaths = zeros( n, mIncrease ) ;
% t vector needs increasing too
% xIncrease = ceil( basic.t0  + basic.tau ) : basic.incr : ceil( basic.t0 + basic.tau ) + ( ( mIncrease - 1 ) * basic.incr ) ;
xIncrease = ceil( basic.t0  + basic.tau ) : basic.incr : ceil( basic.t0 + basic.tau ) + basic.TMax ;
xIncrease = horzcat( xIncrease, ceil( basic.t0 + basic.tau ) + basic.TMax ) ;
% drop the last element in xIncrease if it happens to be equal to the
% previous
if( xIncrease( 1, end ) == xIncrease( 1, end - 1 ) )
    xIncrease = xIncrease( 1, 1:end - 1 ) ;
end


% define matrix to hold cumulative sample size (dimension 1 x 5 to hold the
% index plus the four sample sizes for each of EN, ES, CN and CS)
cumN = zeros( 1, 5 ) ;
% create the resampled paths for the bootstrap
for i = 1:1:n
    % vector of random draws to select the rows of Mn and Mo for each
    % resampled path
    Y = horzcat( round( ( ( ( m - 1 ) + 0.5 ) - 0.5 ) .* rand( 1, mIncrease ) + 0.5 ) ) ;
    % create the resampled paths
    for j = 1:1:mIncrease
        % the first draw ( j = 1 ) should always be the prior mean and
        % effective sample size of the prior
        if j == 1
            cumEN = Mo( j, 2 ) ;
            cumES = Mo( j, 3 ) ;
            cumCN = Mo( j, 4 ) ;
            cumCS = Mo( j, 5 ) ;
            cumDiffC( i, j ) = cumCN - cumCS ;
            cumDiffE( i, j ) = cumEN - cumES ;
            samplePaths( i, j ) = basic.lambda * ( cumEN - cumES ) - ( cumCN - cumCS ) ;
            cumN = Mn( j, : ) ;
        else
            % for subsequent draws: first, exclude the prior mean and sample size
            % information from the matrices
            MnEx = Mn( 2 : size( Mn, 1), : ) ;
            MoEx = Mo( 2 : size( Mo, 1), : ) ;
            % to calculate the weighted means, we need the cumulative sample size
            % at j - 1
            cumNPrev = cumN ;
            cumN = cumNPrev + MnEx( Y( j ), : ) ;
            % Select the E and C for N and S at draw j > 1
            drawEN = MoEx( Y( j ), 2 ) ;
            drawES = MoEx( Y( j ), 3 ) ;
            drawCN = MoEx( Y( j ), 4 ) ;
            drawCS = MoEx( Y( j ), 5 ) ;
            % Select the sample size for each EN, ES, CN and CS
            nEN = MnEx( Y( j ), 2 ) ;
            nES = MnEx( Y( j ), 3 ) ;
            nCN = MnEx( Y( j ), 4 ) ;
            nCS = MnEx( Y( j ), 5 ) ;
            cumEN = ( drawEN * nEN + cumEN * cumNPrev( 1, 2 ) ) / ( nEN + cumNPrev( 1, 2 ) ) ;
            cumES = ( drawES * nES + cumES * cumNPrev( 1, 3 ) ) / ( nES + cumNPrev( 1, 3 ) ) ;
            cumCN = ( drawCN * nCN + cumCN * cumNPrev( 1, 4 ) ) / ( nCN + cumNPrev( 1, 4 ) ) ;
            cumCS = ( drawCS * nCS + cumCS * cumNPrev( 1, 5 ) ) / ( nCS + cumNPrev( 1, 5 ) ) ;
            % Enter the cumulative INMB into the cell for the sample paths
            cumDiffE( i, j ) = ( cumEN - cumES ) ;
            cumDiffC( i, j ) = ( cumCN - cumCS ) ;
            samplePaths( i, j ) = basic.lambda * ( cumEN - cumES ) - ( cumCN - cumCS ) ;
        end
    end
end
% for Clinical Trials R&R: try a different definition of sample path with
% no missing data and observations pairwise-by-pairwise. Set tryAlternative
% to 1 if you wish to do this
if basic.tryMC % basic.tryMC == 1
    mu = Profpath( end ) ;
    samplePathPointRR = mu + basic.sigma .* randn( n, round( basic.TMax, -1) ) ;
    sampleSize = 1:1:round( basic.TMax, -1) ;
    samplePathCumSum = cumsum( samplePathPointRR' ) ;
    samplePathCumSum = samplePathCumSum' ;
    samplePaths =  samplePathCumSum ./ ( basic.t0 + sampleSize ) ;
    % extract the cumulative estimate of the mean every ten observations up to
    % TMax and a bit
    samplePathBlocks = zeros( n, mIncrease-1 ) ;
    for i = 1:1:n
        for j = 1:1:mIncrease-1
            samplePathBlocks( i, j ) = samplePaths( i, j*10 ) ;
        end
    end
    % add the prior mean of zero
    samplePathBlocks = horzcat( zeros( n, 1 ),  samplePathBlocks ) ;
    % report the average of the cumulative sums at the end of the trial ;
    samplePaths = samplePathBlocks ;
end
end
