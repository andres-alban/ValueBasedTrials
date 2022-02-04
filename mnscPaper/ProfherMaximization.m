function [Tstar,ENGstar] = ProfherMaximization(dosave)
% This function plots the figure of the maximization of the ENG that
% appears in section 7 of the paper
%
% Revision: AA 6-11-2017

close all
% Define the model with the appropriate parameters for Profher
basic = ProfherSetParameters();
% Plot the figures
% modifier = {{'online',1,'horizon','Patent'},{'online',1,'horizon','MktExcl'},...
%     {'online',0,'horizon','Patent'},{'online',0,'horizon','MktExcl'}};
% titles = {'Online learning - Patent Protection','Online learning - Market Exclusivity',...
%     'Offline learning - Patent Protection','Offline learning - Market exclusivity'};
% legendloc = {'NorthEast','East','NorthEast','East'};

CTdiscrate = basic.CTdiscrate;
modifier = {{'online',0,'horizon','MktExcl'},{'online',0,'horizon','Patent'}};
line = {'b-','b--'};
marker = {'ro','rs'};
legendloc = {'SouthEast','NorthEast'};
Tstar = zeros(2,length(modifier));
ENGstar = Tstar;
for i = 1:length(modifier)
    basic = MixInputModifier(basic,modifier{i});
    figure(i)
    for j = 1:2
        if j==1
            basic = MixInputModifier(basic,{'CTdiscrate',0});
        elseif j==2
            basic = MixInputModifier(basic,{'CTdiscrate',CTdiscrate});
        end
        % Solve the maximization problem
        ENG = @(T) OneShotExpectedNetGain(basic,T);
        [Tstar(j,i),ENGstar(j,i)] = OneShotMaxExpectedNetGain(basic);
        % Calculate the sample size that would be prescribed by the standard power
        % calculations
        N = FrequentistSampleSize(basic);
        Tfreq = 2*N/basic.r;
        ENGfreq = ENG(Tfreq);
        %Plot the graphs
        T = 0:(1/basic.r):(Tfreq*2.5);
        plot(T,ENG(T),line{j},Tstar(j,i),ENGstar(j,i),marker{j},'LineWidth',2,'MarkerSize',6,'MarkerFaceColor','r');
        hold on
    end
    if i == 1
        hax=axis;
    else
        axis(hax)
    end
    plot([Tfreq Tfreq],hax(3:4),'k:','LineWidth',2);
    hold off
    legend('V(T) undiscounted',...
        'Max. V(T) undiscounted',...
        'V(T) discounted',...
        'Max. V(T) discounted',...
        'Frequentist sample size','Location',legendloc{i})
    xlabel('Recruitment period duration (mo)')
    ylabel('Expected Net Gain (£)')
    % title(titles{i})
    UtilStdizeFigureAbs(i);
end
hold off

% Save the figure if dosave = true
if nargin ~= 1
    dosave = false;
end
if dosave
    UtilSaveFigFile(1, 'Profher','ProfherENGmaximisationMktExcl', [],'eps');
    UtilSaveFigFile(2, 'Profher','ProfherENGmaximisationPatent', [],'eps');
end