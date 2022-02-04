function Profherpn_Tr(dosave)
close all

basic = ProfherSetParameters();
I = [10e6 100e6];
pN = 0:0.03:0.5;
rng(122);

basic = MixInputModifier(basic,{'IN',I(1),'IS',I(1)});
for j = 1:length(pN)
    basic = MixInputModifier(basic,{'pN',pN(j)});
    [Varstarmkt(:,j),ENGstarmkt1(j)] = OneShotMaxExpectedNetGain(basic,true,true);
end
Tstarmkt1 = Varstarmkt(1,:);
rstarmkt1 = Varstarmkt(2,:);
Qstarmkt1 = Tstarmkt1 .* rstarmkt1./2;

sigmaZ = @(Q) sqrt(basic.sigmaX^2 .* 2*Q ./ (basic.n0 .* (2.*basic.n0 + 2*Q)));
Pdisc = @(T) basic.zeta .* (1 - exp(-basic.CTdiscrate .* basic.P(T)./basic.zeta))./basic.CTdiscrate;
alphaN = @(T) basic.IN ./ ((1-pN(14)) .* Pdisc(T));
ratiomkt1 = alphaN(Tstarmkt1(14))./sigmaZ(Qstarmkt1(14));


basic = MixInputModifier(basic,{'horizon','Patent'});
for j = 1:length(pN)
    basic = MixInputModifier(basic,{'pN',pN(j)});
    [Varstarpat(:,j),ENGstarpat1(j)] = OneShotMaxExpectedNetGain(basic,true,true);
end
Tstarpat1 = Varstarpat(1,:);
rstarpat1 = Varstarpat(2,:);
Qstarpat1 = Tstarpat1 .* rstarpat1./2;

sigmaZ = @(Q) sqrt(basic.sigmaX^2 .* 2*Q ./ (basic.n0 .* (2.*basic.n0 + 2*Q)));
Pdisc = @(T) basic.zeta .* (1 - exp(-basic.CTdiscrate .* basic.P(T)./basic.zeta))./basic.CTdiscrate;
alphaN = @(T) basic.IN ./ ((1-pN(14)) .* Pdisc(T));
ratiopat1 = alphaN(Tstarpat1(14))./sigmaZ(Qstarpat1(14));

basic = MixInputModifier(basic,{'horizon','MktExcl','IN',I(2),'IS',I(2)});
for j = 1:length(pN)
    basic = MixInputModifier(basic,{'pN',pN(j)});
    [Varstarmkt(:,j),ENGstarmkt2(j)] = OneShotMaxExpectedNetGain(basic,true,true);
end
Tstarmkt2 = Varstarmkt(1,:);
rstarmkt2 = Varstarmkt(2,:);
Qstarmkt2 = Tstarmkt2 .* rstarmkt2./2;

sigmaZ = @(Q) sqrt(basic.sigmaX^2 .* 2*Q ./ (basic.n0 .* (2.*basic.n0 + 2*Q)));
Pdisc = @(T) basic.zeta .* (1 - exp(-basic.CTdiscrate .* basic.P(T)./basic.zeta))./basic.CTdiscrate;
alphaN = @(T) basic.IN ./ ((1-pN(14)) .* Pdisc(T));
ratiomkt2 = alphaN(Tstarmkt2(14))./sigmaZ(Qstarmkt2(14));

basic = MixInputModifier(basic,{'horizon','Patent'});
for j = 1:length(pN)
    basic = MixInputModifier(basic,{'pN',pN(j)});
    [Varstarpat(:,j),ENGstarpat2(j)] = OneShotMaxExpectedNetGain(basic,true,true);
end
Tstarpat2 = Varstarpat(1,:);
rstarpat2 = Varstarpat(2,:);
Qstarpat2 = Tstarpat2 .* rstarpat2./2;

sigmaZ = @(Q) sqrt(basic.sigmaX^2 .* 2*Q ./ (basic.n0 .* (2.*basic.n0 + 2*Q)));
Pdisc = @(T) basic.zeta .* (1 - exp(-basic.CTdiscrate .* basic.P(T)./basic.zeta))./basic.CTdiscrate;
alphaN = @(T) basic.IN ./ ((1-pN(14)) .* Pdisc(T));
ratiopat2 = alphaN(Tstarpat2(14))./sigmaZ(Qstarpat2(14));

figure(1)
plot(pN,Qstarpat1,'b-',pN,Qstarpat2,'b:',pN,Qstarmkt1,'k--',pN,Qstarmkt2,'k-.','LineWidth',2)
ax = axis;
ax(3) = 0;
axis(ax)
xlabel('Fraction in new treatment (p_N)');ylabel('Optimal sample size Q^*')
lgd = legend(['Q^*_{H}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['Q^*_{H}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],...
    ['Q^*_{P}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['Q^*_{P}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],'Location','SouthEast');
% set(lgd,'Position',[0.56 0.28 0.29 0.25]);
UtilStdizeFigureAbs(1);

figure(2)
plot(pN,Tstarpat1,'b-',pN,Tstarpat2,'b:',pN,Tstarmkt1,'k--',pN,Tstarmkt2,'k-.','LineWidth',2)
ax = axis;
ax(3) = 0;
axis(ax)
xlabel('Fraction in new treatment (p_N)');ylabel('Optimal recruitment duration T^*')
lgd = legend(['T^*_{H}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['T^*_{H}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],...
    ['T^*_{P}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['T^*_{P}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],'Location','SouthEast');
% set(lgd,'Position',[0.56 0.52 0.29 0.25]);
UtilStdizeFigureAbs(2);

figure(3)
plot(pN,rstarpat1,'b-',pN,rstarpat2,'b:',pN,rstarmkt1,'k--',pN,rstarmkt2,'k-.','LineWidth',2)
ax = axis;
ax(3) = 0;
ax(4) = ax(4)+1;
axis(ax)
xlabel('Fraction in new treatment (p_N)');ylabel('Optimal recruitment rate r^*')
lgd = legend(['r^*_{H}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['r^*_{H}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],...
    ['r^*_{P}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['r^*_{P}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],'Location','SouthEast');
% set(lgd,'Position',[0.57 0.27 0.29 0.25]);
UtilStdizeFigureAbs(3);

figure(4)
plot(pN,ENGstarpat1./1e6,'b-',pN,ENGstarpat2./1e6,'b:',pN,ENGstarmkt1./1e6,'k--',pN,ENGstarmkt2./1e6,'k-.','LineWidth',2)
ax = axis;
ax(3) = 0;
axis(ax)
xlabel('Fraction in new treatment (p_N)');ylabel('Optimal expected net gain V^* (millions of £)')
lgd = legend(['V^*_{H}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['V^*_{H}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],...
    ['V^*_{P}, I_N=I_S=',num2str(I(1)/1e6),'Mill.'],['V^*_{P}, I_N=I_S=',num2str(I(2)/1e6),'Mill.'],'Location','East');
set(lgd,'Position',[0.55 0.45 0.29 0.25]);
UtilStdizeFigureAbs(4);

ratiopat1
ratiopat2
ratiomkt1
ratiomkt2

% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'Profherpn_Q',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'Profherpn_T',[],'eps');
    UtilSaveFigFile(3, 'Profher', 'Profherpn_r',[],'eps');
    UtilSaveFigFile(4, 'Profher', 'Profherpn_ENG',[],'eps');
end








