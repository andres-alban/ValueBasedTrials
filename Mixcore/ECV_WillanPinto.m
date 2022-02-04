%% Willan Pinto
basic = MixInputConstructor({'mu0',68.97,'sigmaX',sqrt(2*217227),'n0',2*217227/3724.78,...
                    'Population',1e6,'cfix',5e5,'c',1600,'pN',0,'IS',0,...
                    'IN',0,'zeta',1,'H',1e6,'Tmax',1e6,'Delta',0,'r',1,...
                    'horizon','Patent','online',1});
                
ENG = @(x) OneShotExpectedNetGain(basic,x);
x = 0:1000;

figure(1)
plot(x,ENG(x) - ENG(0))
[m,i] = max(ENG(x) - ENG(0))

%% Eckerman Willan (Time dependence)
basic = MixInputConstructor({'mu0',68.97,'sigmaX',sqrt(2*217227),'n0',2*217227/3724.78,...
                    'Population',1e6,'cfix',5e5,'c',1600,'pN',0,'IS',0,...
                    'IN',0,'zeta',50000,'H',1e6/5e4,'Tmax',19.5,'Delta',0.5,'r',500    ,...
                    'horizon','MktExcl','online',1});
ENG = @(x) OneShotExpectedNetGain(basic,x);
x = 0:0.05:19;

figure(2)
plot(x,ENG(x) - ENG(0))
[m,i] = max(ENG(x) - ENG(0))

