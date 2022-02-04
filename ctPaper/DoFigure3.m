function DoFigure3(basic, advanced, mat125, mat250, tProf, profPath, tboot, samplepathsboot, fignum)

set(gca,'FontSize',10)
plot_upper = 10000 ; % set upper limit of plot region ( + / 
plot_lower = -10000  ; % set lower limit of plot region
ylim( [ plot_lower,plot_upper] ) ;
figure(fignum); % always increment fignum and create a new figure
hold on ; 

p1 = plot( mat125.tvec, mat125.bndupper, 'r', 'Linewidth', basic.linewid ) ;
p2 = plot( mat125.tvec, mat125.bndlower, 'r', 'Linewidth', basic.linewid ) ;
p3 = plot( mat250.tvec,  mat250.bndupper, 'b--', 'Linewidth', basic.linewid ) ;
p4 = plot( mat250.tvec,  mat250.bndlower, 'b--', 'Linewidth', basic.linewid ) ;
xlim( [ 0, 300] ) ;
ylim( [ plot_lower,plot_upper] ) ;

% Plot the ProFHER path on the boundaries diagram and set up the bootstrap
p5 = plot( tProf, profPath,  'k-o', 'Linewidth', basic.linewid ) ; 
p6 = line( [ basic.t0, basic.t0 ], [ plot_lower, plot_upper ], 'Linestyle', ':', 'Color', 'k', 'Linewidth', basic.linewid ) ;
p7 = line( [ basic.t0 + basic.tau, basic.t0 + basic.tau ], [ plot_lower, plot_upper ], 'Linestyle', ':', 'Color', 'k', 'Linewidth', basic.linewid ) ;
legend( [ p4, p1, p5], {'Stopping boundary, Q_{max} = 250', 'Stopping boundary, Q_{max} = 125',  'ProFHER posterior mean'}, 'Fontsize', 12 ) ;

tProf2aa = tboot(1, tboot(1,:) > 0);
tProf2bb = tboot(2, tboot(2,:) > 0);
tProf2cc = tboot(3, tboot(3,:) > 0);
samplePathProf2aa = samplepathsboot(1,1:length(tProf2aa));
samplePathProf2bb = samplepathsboot(2,1:length(tProf2bb));
samplePathProf2cc = samplepathsboot(3,1:length(tProf2cc));
plot( tProf2aa, samplePathProf2aa, 'k--+', 'DisplayName', 'Resampled path 1', 'Color', 'magenta', 'LineWidth', 2  ) ;
plot( tProf2bb, samplePathProf2bb, 'k--s', 'DisplayName', 'Resampled path 2', 'Color', 'green', 'LineWidth', 2  ) ;
plot( tProf2cc, samplePathProf2cc, 'k--o', 'DisplayName', 'Resampled path 3', 'Color', 'cyan', 'LineWidth', 2   ) ;
ylabel('Prior/posterior mean for E[INMB]','FontSize',16,'FontName','Helvetica');
xlabel('n_0 + pairwise allocations','FontSize',16,'FontName','Helvetica'); 
UtilStdizeFigure_Sec4LargeFont(fignum,advanced,true);
UtilSaveFigEpsPdf(fignum,'Figure',strcat('pathsCT',''),'-r600');
hold off; 

end