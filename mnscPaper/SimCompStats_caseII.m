% This script takes ~1' to process ~11,000 instances
load('explore_caseII.mat')
% Make sure to update basic definition from
% explore_parameter_space_caseII.m
basic = MixInputConstructor({'c',1000,'pN',0,'mu0',0,'IS',0,'IN',0,...
    'zeta',100,'H',200,'Tmax',200,'rmax',200,'Delta',0.1*200,'Population',10000,'n0',1,'sigmaX',1000,...
    'cfix',0,'cr',1,'cT',0,'r',1,'horizon', 'MktExcl','CTdiscrate',log(1+0.01)/100,...
    'online',0,'relevantdifference',1,'alpha',0.05,'power',0.8,...
    'UnkVariance',false,'UnkVarianceShape',3});
H = 200;
c = 1000;
sigmaX = 1000;
results_pos = results(results(:,15)~=0 & results(:,13)<H-results(:,10)-0.01,:);
combs = size(results_pos,1);
% initialize arrays to save derivative data
comp_stats_c = zeros(combs,5);
comp_stats_PH = zeros(combs,5);
comp_stats_CTdiscrate = zeros(combs,5);
comp_stats_n0 = zeros(combs,5);
comp_stats_sigmaX = zeros(combs,5);
comp_stats_Delta = zeros(combs,5);

tic()
for i = 1:size(results_pos,1)
    zeta = results_pos(i,11);
    basic.zeta = zeta;
    Delta = results_pos(i,10);
    basic.Delta = Delta;
    basic.Tmax = basic.H-Delta;
    if results_pos(i,1)
        basic = MixInputModifier(basic,{'horizon','MktExcl','Population',zeta*H});
    else
        basic = MixInputModifier(basic,{'horizon','Patent'});
    end
    cr = results_pos(i,2);
    theta = results_pos(i,3);
    basic.ccap = @(r) cr*r.^(theta);
    CTdiscrate = results_pos(i,4);
    basic.CTdiscrate = CTdiscrate;
    n0 = results_pos(i,5);
    basic.n0 = n0;
    mu0 = results_pos(i,6);
    basic.mu0 = mu0;
    IN = results_pos(i,7);
    basic.IN = IN*zeta*H*sigmaX/sqrt(n0);
    IS = results_pos(i,8);
    basic.IS = IS*zeta*H*sigmaX/sqrt(n0);
    pN = results_pos(i,9);
    basic.pN = pN;
    
    
    T = results_pos(i,13);
    r = results_pos(i,14);
    % compute derivatives and save to arrays
    % PH
    if results_pos(i,1)
        [dENG,dT,dr,dQ,VTT] = OneShotCompStats_caseIandII(basic,T,r,'Population',1,0.01);
    else
        [dENG,dT,dr,dQ,VTT] = OneShotCompStats_caseIandII(basic,T,r,'H',0.1,0.01);
    end
    comp_stats_PH(i,1) = dENG; comp_stats_PH(i,2) = dT; comp_stats_PH(i,3) = dr; comp_stats_PH(i,4) = dQ; comp_stats_PH(i,5) = VTT;
    % CTdiscrate
    [dENG,dT,dr,dQ,VTT] = OneShotCompStats_caseIandII(basic,T,r,'CTdiscrate',1e-6,0.01);
    comp_stats_CTdiscrate(i,1) = dENG; comp_stats_CTdiscrate(i,2) = dT; comp_stats_CTdiscrate(i,3) = dr; comp_stats_CTdiscrate(i,4) = dQ; comp_stats_CTdiscrate(i,5) = VTT;
    % n0
    [dENG,dT,dr,dQ,VTT] = OneShotCompStats_caseIandII(basic,T,r,'n0',1e-3,0.01);
    comp_stats_n0(i,1) = dENG; comp_stats_n0(i,2) = dT; comp_stats_n0(i,3) = dr; comp_stats_n0(i,4) = dQ; comp_stats_n0(i,5) = VTT;
    % sigmaX
    [dENG,dT,dr,dQ,VTT] = OneShotCompStats_caseIandII(basic,T,r,'sigmaX',1,0.01);
    comp_stats_sigmaX(i,1) = dENG; comp_stats_sigmaX(i,2) = dT; comp_stats_sigmaX(i,3) = dr; comp_stats_sigmaX(i,4) = dQ; comp_stats_sigmaX(i,5) = VTT;
    % c
    [dENG,dT,dr,dQ,VTT] = OneShotCompStats_caseIandII(basic,T,r,'c',1,0.01);
    comp_stats_c(i,1) = dENG; comp_stats_c(i,2) = dT; comp_stats_c(i,3) = dr; comp_stats_c(i,4) = dQ; comp_stats_c(i,5) = VTT;
    % Delta
    [dENG,dT,dr,dQ,VTT] = OneShotCompStats_caseIandII(basic,T,r,'Delta',0.1,0.01);
    comp_stats_Delta(i,1) = dENG; comp_stats_Delta(i,2) = dT; comp_stats_Delta(i,3) = dr; comp_stats_Delta(i,4) = dQ; comp_stats_Delta(i,5) = VTT;



    if mod(i,1e3) == 0
        fprintf('%d/%d\n',i,combs);
    end
end
toc()

% Save data
save('SimCompStats_caseII.mat','results_pos','comp_stats_PH','comp_stats_CTdiscrate','comp_stats_n0','comp_stats_sigmaX','comp_stats_c','comp_stats_Delta')

% print results
fileID = fopen('comp_stats_caseII.txt','w');
fprintf(fileID,' Parameter     dENG    dT    dr    dQ\n');
% P
x = comp_stats_PH(results_pos(:,1)==1,:);
x = x(~isnan(x(:,1)),:);
fprintf(fileID,'%10s %8.3f %5.3f %5.3f %5.3f\n','P', sum(x(:,1)>0)/size(x,1), sum(x(:,2)>0)/size(x,1), sum(x(:,3)>0)/size(x,1), sum(x(:,4)>0)/size(x,1));
% H
x = comp_stats_PH(results_pos(:,1)~=1,:);
x = x(~isnan(x(:,1)),:);
fprintf(fileID,'%10s %8.3f %5.3f %5.3f %5.3f\n','H', sum(x(:,1)>0)/size(x,1), sum(x(:,2)>0)/size(x,1), sum(x(:,3)>0)/size(x,1), sum(x(:,4)>0)/size(x,1));
% CTdiscrate
x = comp_stats_CTdiscrate;
x = x(~isnan(x(:,1)),:);
fprintf(fileID,'%10s %8.3f %5.3f %5.3f %5.3f\n','CTdiscrate', sum(x(:,1)>0)/size(x,1), sum(x(:,2)>0)/size(x,1), sum(x(:,3)>0)/size(x,1), sum(x(:,4)>0)/size(x,1));
% n0
x = comp_stats_n0;
x = x(~isnan(x(:,1)),:);
fprintf(fileID,'%10s %8.3f %5.3f %5.3f %5.3f\n','n0', sum(x(:,1)>0)/size(x,1), sum(x(:,2)>0)/size(x,1), sum(x(:,3)>0)/size(x,1), sum(x(:,4)>0)/size(x,1));
% sigmaX
x = comp_stats_sigmaX;
x = x(~isnan(x(:,1)),:);
fprintf(fileID,'%10s %8.3f %5.3f %5.3f %5.3f\n','sigmaX', sum(x(:,1)>0)/size(x,1), sum(x(:,2)>0)/size(x,1), sum(x(:,3)>0)/size(x,1), sum(x(:,4)>0)/size(x,1));
% c
x = comp_stats_c;
x = x(~isnan(x(:,1)),:);
fprintf(fileID,'%10s %8.3f %5.3f %5.3f %5.3f\n','c', sum(x(:,1)>0)/size(x,1), sum(x(:,2)>0)/size(x,1), sum(x(:,3)>0)/size(x,1), sum(x(:,4)>0)/size(x,1));
% Delta
x = comp_stats_Delta;
x = x(~isnan(x(:,1)),:);
fprintf(fileID,'%10s %8.3f %5.3f %5.3f %5.3f\n','Delta', sum(x(:,1)>0)/size(x,1), sum(x(:,2)>0)/size(x,1), sum(x(:,3)>0)/size(x,1), sum(x(:,4)>0)/size(x,1));

fclose(fileID);








