load('explore_caseIV.mat')
% Make sure to update basic definition from
% explore_parameter_space_caseIV.m
basic = MixInputConstructor({'c',1000,'pN',0,'mu0',0,'IS',0,'IN',0,...
    'zeta',100,'H',200,'Tmax',200,'rmax',200,'Delta',0.1*200,'Population',10000,'n0',1,'sigmaX',1000,...
    'cfix',0,'cr',1,'cT',0,'r',1,'horizon', 'MktExcl','CTdiscrate',log(1+0.01)/100,...
    'online',0,'relevantdifference',1,'alpha',0.05,'power',0.8,...
    'UnkVariance',false,'UnkVarianceShape',3});
H = 200;
c = 1000;
sigmaX = 1000;
results_pos = results(results(:,15)~=0 & results(:,14)<basic.rmax-0.01 & results(:,13)<H-results(:,10)-0.01,:);
combs = size(results_pos,1);

results_pos = [results_pos zeros(combs,4)];

rng(121)
randomsample = datasample(1:combs,10000,'Replace',false);
results_pos = results_pos(randomsample,:);

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
    
    
    
    [Design,ENGstar] = OneShotMaxExpectedNetGain(basic,true,'Tr',100,false);
    results_pos(i,16) = ENGstar;
    results_pos(i,17) = Design(1);
    results_pos(i,18) = Design(2);
    results_pos(i,19) = Design(1)*Design(2)/2;


    if mod(i,1e3) == 0
        fprintf('%d/%d\n',i,size(results_pos,1));
    end
end
toc()

save('explore_caseIV_validation.mat','results_pos')

sum(results_pos(:,16) >  results_pos(:,12) + 1e3)
sum(results_pos(:,12) >  results_pos(:,16) + 1e3)

