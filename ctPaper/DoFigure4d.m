function DoFigure4d(basic,advanced,inmbStop125, inmbStop250, fignum)

figure(fignum)
hold on
h = histogram( inmbStop125, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'r', 'Linewidth', basic.linewid ) ;
h.LineWidth = basic.linewid ;

h = histogram( inmbStop250, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'b', 'Linewidth', basic.linewid, 'LineStyle', '--' ) ;
h.LineWidth = basic.linewid ;
legend( 'Q_{max} = 125', 'Q_{max} = 250' ) ;

xlabel( 'Posterior mean for E[INMB] at stopping', 'Fontsize', advanced.bigfontsize ) ; 
ylabel( 'Relative frequency', 'Fontsize', advanced.bigfontsize ) ; 
xlim( [ -6500, 5000] ) ;
UtilStdizeFigure_Sec4LargeFont(fignum,advanced,true);

UtilSaveFigEpsPdf(fignum,'Figure',strcat('histINMBStopCT',''),'-r600') ;
