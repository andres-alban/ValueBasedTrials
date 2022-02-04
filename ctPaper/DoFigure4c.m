function DoFigure4c(basic,advanced,trialSampleSize125,trialSampleSize250,fignum)

figure(fignum)
hold on
[ f1, ~ ] = cdfplot( trialSampleSize125 ) ;
set( f1, 'LineStyle', '-', 'Linewidth', basic.linewid, 'Color', 'r' ) ;

[ f2, ~ ] = cdfplot( trialSampleSize250 ) ;
set( f2, 'LineStyle', '--', 'Linewidth', basic.linewid, 'Color', 'b' ) ;

legend('Q_{max} = 125', 'Q_{max} = 250' ) ;
title('')
xlabel( 'Number of pairwise allocations made before stopping', 'Fontsize', advanced.bigfontsize ) ; 
ylabel( 'Cumulative relative frequency', 'Fontsize', advanced.bigfontsize ) ; 
ylim( [ 0, 1 ] ) ; 
UtilStdizeFigure_Sec4LargeFont(fignum,advanced,true);

UtilSaveFigEpsPdf(fignum,'Figure',strcat('cdfTrialSampleSizeCT',''),'-r600') ;