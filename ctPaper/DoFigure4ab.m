function DoFigure4ab(basic,advanced,trialSampleSize,fignum)
figure( fignum ); % always increment fignum and create a new figure
if (basic.TMax == 125)
    edges = [57 62:5:122 127] ;
    xlim( [ 50, 130 ] ) ;
else
    edges = [57 62:5:252 257] ;
    xlim( [ 50, 260 ] ) ;
end 
h = histogram( trialSampleSize,  edges ) ; 
h.Normalization = 'probability' ; 
h.LineWidth = basic.linewid ;
xlabel( 'Number of pairwise allocations made before stopping', 'Fontsize', advanced.bigfontsize ) ; 
ylabel( 'Relative frequency', 'Fontsize', advanced.bigfontsize ) ; 
if (basic.TMax == 125)
    xlim( [ 50, 130 ] ) ;
else
    xlim( [ 50, 260 ] ) ;
end 
UtilStdizeFigure_Sec4LargeFont(fignum,advanced,true);
if( basic.TMax == 125 ) 
    UtilSaveFigEpsPdf(fignum,'Figure',strcat('histTrialSampleSize125CT',''),'-r600');
else
    UtilSaveFigEpsPdf(fignum,'Figure',strcat('histTrialSampleSize250CT',''),'-r600');
end 
hold off ; 