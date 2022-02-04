function [Delaybasic,Delayadvanced] = Mix2DelayTransformer(basic)
% This functions transforms the problem specified by basic (output from
% MixInputConstructor) into the sequential problem as explained in the
% appendix of the Profher theory paper. When the transformation is not
% possible, an error is displayed.

% To-do:
% > When ICost < 0, DelayInputValidator.m will stop the code unnecessarily
% > Testing

% Take parameters from the structure basic
c = basic.c;
pN = basic.pN;
mu0 = basic.mu0;
IS = basic.IS;
IN = basic.IN;
H = basic.H;
Tmax = basic.Tmax;
Delta = basic.Delta;
n0 = basic.n0;
sigmaX = basic.sigmaX;
cfix = basic.cfix;
r = basic.r;
Population = basic.Population;
online = basic.online;
CTdiscrate = basic.CTdiscrate;
UnkVariance = basic.UnkVariance;
xi0 = basic.UnkVarianceShape;
horizon = basic.horizon;

N = FrequentistSampleSize(basic);

if ~strcmp(horizon,'MktExcl')
    error('The problem is not transformable for patent horizon.')
end

DTdiscrate = exp(CTdiscrate) - 1;
theta = (1 + DTdiscrate)^(-1/r);

if pN == 0
    [Delaybasic,Delayadvanced] = DelayInputConstructor({'c',2*c,'theta',theta,'PPatients',Population,...  % check definition of theta
        'ICost',IN,'tau',ceil(r*Delta),'TMax',floor(r*Tmax/2),'numpairs',N,...
        'online',online,'mu0',mu0,'t0',n0,'sigma',sigmaX},...
        {'UnkVariance',UnkVariance},'UnkVarianceShape',xi0);
elseif IN == 0 && IS == 0 && mu0 == 0
    if pN < 1/2
        [Delaybasic,Delayadvanced] = DelayInputConstructor({'c',2*c/(1 - 2*pN),'theta',theta,'PPatients',Population/(1 - 2*pN),...  % check definition of theta
            'ICost',0,'tau',ceil(r*Delta),'TMax',floor(r*Tmax/2),'numpairs',N,...
            'online',online,'mu0',mu0,'t0',n0,'sigma',sigmaX},...
            {'UnkVariance',UnkVariance,'UnkVarianceShape',xi0});
    elseif pN == 1/2
        [Delaybasic,Delayadvanced] = DelayInputConstructor({'c',2*c,'theta',theta,'PPatients',Population,...  % check definition of theta
            'ICost',0,'tau',ceil(r*Delta),'TMax',floor(r*Tmax/2),'numpairs',N,...
            'online',online,'mu0',mu0,'t0',n0,'sigma',sigmaX},...
            {'UnkVariance',UnkVariance,'UnkVarianceShape',xi0});
    end
else %Force decision to standard or new treatment
    warning('The problem cannot be transformed exactly but it will be transformed by forcing the decision to the standard or new treatments')
    if pN < 1/2
        [Delaybasic,Delayadvanced] = DelayInputConstructor({'c',2*c/(1 - 2*pN),'theta',theta,'PPatients',Population/(1 - 2*pN),...  % check definition of theta
            'ICost',(IN - IS)/(1 - 2*pN),'tau',ceil(r*Delta),'TMax',floor(r*Tmax/2),'numpairs',N,...
            'online',online,'mu0',mu0,'t0',n0,'sigma',sigmaX},...
            {'UnkVariance',UnkVariance,'UnkVarianceShape',xi0});
    elseif pN == 1/2
        [Delaybasic,Delayadvanced] = DelayInputConstructor({'c',2*c,'theta',theta,'PPatients',Population,...  % check definition of theta
            'ICost',IN - IS,'tau',ceil(r*Delta),'TMax',floor(r*Tmax/2),'numpairs',N,...
            'online',online,'mu0',mu0,'t0',n0,'sigma',sigmaX},...
            {'UnkVariance',UnkVariance,'UnkVarianceShape',xi0});
    end
end
