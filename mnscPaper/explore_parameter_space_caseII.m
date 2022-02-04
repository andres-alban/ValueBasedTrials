% ~15,000 instances and ~25'
rng(121); % for reproducibility. The maximization has random starting point
basic = MixInputConstructor({'c',1000,'pN',0,'mu0',0,'IS',0,'IN',0,...
    'zeta',100,'H',200,'Tmax',200,'rmax',200,'Delta',0.1*200,'Population',10000,'n0',1,'sigmaX',1000,...
    'cfix',0,'cr',1,'cT',0,'r',1,'horizon', 'MktExcl','CTdiscrate',log(1+0.01)/100,...
    'online',0,'relevantdifference',1,'alpha',0.05,'power',0.8,...
    'UnkVariance',false,'UnkVarianceShape',3});
basic.r = basic.rmax;
horizon = {'MktExcl','Patent'};
cr = [0];
theta = [1];
H = 200;
CTdiscrate = [0,0.5,0.1,0.01]./H;
sigmaX = 1000;
n0 = [1,0.1,10];
mu0 = [0,-1000,1000];
c = 1000;
IN = [0,0.01,0.1];
IS = IN;
pN = [0,0.2,0.4];
Delta = [0.1,0.01,0.3]*H;
zeta = [1e5,1e4,1e6];

combs = length(horizon)*length(cr)*length(theta)*length(H)*length(CTdiscrate)*length(sigmaX)...
    *length(n0)*length(mu0)*length(c)*length(IN)*length(IS)*length(pN)...
    *length(Delta)*length(zeta);

% results = array2table(zeros(combs,15),...
%     'VariableNames',{'FixedPool','cr','theta','CTdiscrate','n0','mu0','IN','IS',...
%     'pN','Delta','zeta','ENGstar','Tstar','rstar','Qstar'});
results = zeros(combs,15);


count = 1;
tic()
for horizon_i = 1:length(horizon)
for cr_i = 1:length(cr)
for theta_i = 1:length(theta)
for CTdiscrate_i = 1:length(CTdiscrate)
for n0_i = 1:length(n0)
for mu0_i = 1:length(mu0)
for IN_i = 1:length(IN)
for IS_i = 1:length(IS)
for pN_i = 1:length(pN)
for Delta_i = 1:length(Delta)
for zeta_i = 1:length(zeta)
    if strcmp(horizon{horizon_i},'MktExcl') && CTdiscrate(CTdiscrate_i) == 0
        combs = combs - 1;
        results = results(1:end-1,:);
        continue
    end
    results(count,1:11) = [strcmp(horizon{horizon_i},'MktExcl'),cr(cr_i),theta(theta_i),...
        CTdiscrate(CTdiscrate_i),n0(n0_i),mu0(mu0_i),IN(IN_i),IS(IS_i),pN(pN_i),...
        Delta(Delta_i),zeta(zeta_i)];
    count = count + 1;
end
end
end
end
end
end
end
end
end
end
end
toc()

tic()
for i = 1:combs
    zeta = results(i,11);
    basic.zeta = zeta;
    Delta = results(i,10);
    basic.Delta = Delta;
    basic.Tmax = basic.H-Delta;
    if results(i,1)
        basic.P = @(T) zeta*H;
    else
        basic.P = @(T) zeta.*(H-Delta.*(T>0)-T);
    end
    cr = results(i,2);
    theta = results(i,3);
    basic.ccap = @(r) cr*r.^(theta);
    CTdiscrate = results(i,4);
    basic.CTdiscrate = CTdiscrate;
    n0 = results(i,5);
    basic.n0 = n0;
    mu0 = results(i,6);
    basic.mu0 = mu0;
    IN = results(i,7);
    basic.IN = IN*zeta*H*sigmaX/sqrt(n0);
    IS = results(i,8);
    basic.IS = IS*zeta*H*sigmaX/sqrt(n0);
    pN = results(i,9);
    basic.pN = pN;
    [Design,ENGstar] = OneShotMaxExpectedNetGain(basic,true,'T',10,false);
    results(i,12) = ENGstar;
    results(i,13) = Design;
    results(i,14) = basic.r;
    results(i,15) = results(i,13)*results(i,14)/2;
    if mod(i,1e3) == 0
        fprintf('%d/%d\n',i,combs);
    end
end
toc()
T = array2table(results,...
    'VariableNames',{'FixedPool','cr','theta','CTdiscrate','n0','mu0','IN','IS',...
    'pN','Delta','zeta','ENGstar','Tstar','rstar','Qstar'});

save('explore_caseII.mat','results','T')



