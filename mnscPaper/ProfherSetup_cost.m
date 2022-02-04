% Shape of setup cost assuming linear cost in sites
function bestfit_ab = ProfherSetup_cost(dosave)
close all
load 'Recruitment_by_site.mat' Recruitmentbysite
pat_per_site = Recruitmentbysite.Randomised;
pat_per_site = sort(pat_per_site,'descend');
% pat_per_site = pat_per_site(pat_per_site > 0);
% pat_per_site = pat_per_site(1:(end-5));
cum_pat = cumsum(pat_per_site);
cum_recrate = cum_pat/32; % the trial was 32 months long
s = (1:length(cum_recrate))';


bestfit2 = polyfit(cum_recrate,s,2);
% plot(cum_recrate,s,cum_recrate,polyval(bestfit2,cum_recrate))

bestfit_sq = sum(s.*cum_recrate.^2)/sum(cum_recrate.^4);
% plot(cum_recrate,s,cum_recrate,bestfit_sq*cum_recrate.^2)

bestfit_cube = sum(s.*cum_recrate.^3)/sum(cum_recrate.^6);
% plot(cum_recrate,s,cum_recrate,bestfit_cube*cum_recrate.^3)

bestfit3 = polyfit(cum_recrate,s,3);
% plot(cum_recrate,s,cum_recrate,polyval(bestfit3,cum_recrate))

bestfit_exp = sum(s.*exp(cum_recrate))/sum(exp(cum_recrate).^2);
% plot(cum_recrate,s,cum_recrate,bestfit_exp*exp(cum_recrate))

opt = optimoptions('fminunc','MaxIterations',1000,'MaxFunctionEvaluations',1000,...
    'StepTolerance',1e-12);
bestfit_ab = fminunc(@(x) sum((s - x(1)*cum_recrate.^x(2)).^2),[0.0032;2],opt);
% plot(cum_recrate,s,cum_recrate,bestfit_ab(1)*cum_recrate.^bestfit_ab(2))

cost_per_site = 480000/35; % There were 35 sites involved in ProFHER
r = 0:0.1:10;
ccap = cost_per_site*bestfit_ab(1)*r.^bestfit_ab(2);
figure(1)
plot(cum_recrate,s*cost_per_site,'o-r',r,ccap,'b','LineWidth',2,'MarkerSize',6,'MarkerFaceColor','r')
xlabel('Recruitment rate r (mo^{-1})')
ylabel('Variable cost of recruitment rate (£)')
legend('ProFHER data','Best fit of data','Location','NorthWest')
UtilStdizeFigureAbs(1);

if nargin < 1
    dosave = false;
end
if dosave
UtilSaveFigFile(1, 'Profher','ProfherCostRecRate', [],'eps')
end
