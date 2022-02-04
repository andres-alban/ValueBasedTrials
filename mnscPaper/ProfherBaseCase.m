function ProfherBaseCase(dosave)
close all
basic = ProfherSetParameters();
h = 0.3;
r = 0:h:20;
T = 0:h:33;

[rr,TT] = meshgrid(r,T);
ENG = zeros(size(rr));

for i = 1:length(T)
    for j = 1:length(r)
        ENG(i,j) = OneShotExpectedNetGain(basic,TT(i,j),rr(i,j));
    end
end
[astar,ENGstar] = OneShotMaxExpectedNetGain(basic,true,true);

figure(1)
Levels = [92,91,90,88,80];
[C,h] = contour(TT,rr,ENG./1e6,Levels,'Linewidth',2);
clabel(C,h,'FontSize',12)
hold on
p1 = plot(astar(1),astar(2),'xr','LineWidth',2,'MarkerSize',10);
text(astar(1),astar(2)-0.5,'92.3','FontSize',12)
p2 = plot(32,7.833,'sr','LineWidth',5,'MarkerSize',5);
xlabel('Recruitment duration T (mo)')
ylabel('recruitment rate r (mo^{-1})')
title('Expected net gain V (millions of £)')
legend([p1,p2],'Optimal design','Actual design')
UtilStdizeFigureAbs(1);

% Number of sites required
(basic.ccap(astar(2))-basic.cfix)/13700

% For fixed horizon
basic = MixInputModifier(basic,{'horizon','Patent'});

[rr,TT] = meshgrid(r,T);
ENG = zeros(size(rr));

for i = 1:length(T)
    for j = 1:length(r)
        ENG(i,j) = OneShotExpectedNetGain(basic,TT(i,j),rr(i,j));
    end
end
[astar,ENGstar] = OneShotMaxExpectedNetGain(basic,true,true);

figure(2)
Levels = [84,82,80,75,70];
[C,h] = contour(TT,rr,ENG./1e6,Levels,'Linewidth',2);
clabel(C,h,'FontSize',12)
hold on
p1 = plot(astar(1),astar(2),'xr','LineWidth',2,'MarkerSize',10);
text(astar(1),astar(2)-0.5,'85.0','FontSize',12)
p2 = plot(32,7.833,'sr','LineWidth',6,'MarkerSize',6);
xlabel('Recruitment duration T (mo)')
ylabel('recruitment rate r (mo^{-1})')
title('Expected net gain V (millions of £)')
legend([p1,p2],'Optimal design','Actual design')
UtilStdizeFigureAbs(2);

% Number of sites required
(basic.ccap(astar(2))-basic.cfix)/13700

% Save the figures if dosave=true
if nargin < 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher', 'ProfherBase_mkt',[],'eps');
    UtilSaveFigFile(2, 'Profher', 'ProfherBase_pat',[],'eps');
end

