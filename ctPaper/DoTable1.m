function DoTable1(basic,out,outQmax)

% disp( 'Parameter values' ) ; 
% parVars = horzcat( basic.lambda, basic.pN, basic.P, basic.PPatients, basic.sigma, basic.t0, basic.ICostN, basic.ICostS, basic.cResearch,...
% basic.c, basic.tau, basic.TMax ) ; 
% U = table( parVars', 'VariableNames', { 'value' }, 'RowNames',...
%     { 'Max WTP for one unit of effectiveness (lambda)',  'Proportion on New (pN)',...
%     'Size of population to benefit (P)', 'Adjusted size of population to benefit', 'Sampling standard deviation (sigma)',...
%     'Effective sample size of prior (t0)', 'Cost of switching to New (I_N)', 'Cost of switching to Standard (I_S)',...
%     'Actual cost per pairwise allocation (c)', 'Adjusted cost per pairwise allocation', 'Delay in pairwise allocations (tau)',...
%     'Maximum number of pairwise allocations (TMax)' } ) ;
% disp( U ) ;

if (basic.tryMC)
    method = 'Monte Carlo';
else
    method = 'bootstrap';
end

mx = max( out ) ; 
mn = min( out ) ; 
av = mean( out ) ; 
s = std( out ) ; 
disp( ['Operating characteristics if trial uses boundaries with Q_{max}=',num2str(basic.TMax),' (',method,')'] ) ;
U = table( av', s', mn', mx', 'VariableNames', { 'mean', 'sd', 'min', 'max' }, 'RowNames',...
    { 'Sample size', 'Savings in budget', 'Posterior mean for cost-effectiveness', 'Final decision sling',...
    'Final decision surgery', 'Cross lower boundary', 'Cross upper boundary', 'Cross upper boundary but final decision sling',...
    'Cross lower boundary but final decision surgery',  'Time at crossing lower boundary', 'Time at crossing upper boundary' } ) ;
disp( U ) ;

mx = max( outQmax ) ; 
mn = min( outQmax ) ; 
av = mean( outQmax ) ; 
s = std( outQmax ) ; 
disp( ['Operating characteristics if trial runs to Q_{max}=', num2str(basic.TMax),' (',method,')']) ;
U = table( av', s', mn', mx', 'VariableNames', { 'mean', 'sd', 'min', 'max' }, 'RowNames',...
    { 'Posterior mean of cost-effectiveness', 'Final decision sling', 'Final decision surgery' } ) ;
disp( U ) ;