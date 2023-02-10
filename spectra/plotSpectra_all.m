function plotSpectra_all(Period,amax,vmax, umax, fpath,fname)
% Plot response spectra
figure(2) ; clf ;
fsz = 16;
lw = 2;
set(gcf,'position', [50 200 1300 450]) ;


subplot_special(3,1, fsz ) ; hold;
loglog(Period, amax ,'-k','LineWidth',lw); hold;
ylabel('Pseudo spectral acceleration (cm/s^2)','FontSize',fsz,'FontWeight','normal','FontAngle','normal','fontname','times');
title( sprintf( 'Spectal acceleration')) ;
axis ("square") ; xlim([min(Period) max(Period)*1.25]); ylim([min(amax)/10 max(amax)*10]);
box on; grid on;


subplot_special(3,2, fsz ) ; hold ;
loglog(Period, vmax ,'-k','LineWidth',lw); hold;
ylabel('Spectral velocity (cm/s)','FontSize',fsz,'FontWeight','normal','FontAngle','normal','fontname','times');
title( sprintf( 'Spectal velocity')) ;
axis square; xlim([min(Period) max(Period)*1.25]);  ylim([min(vmax)/10 max(vmax)*10]);
box on; grid on;


subplot_special(3 , 3 , fsz ) ; hold ;
loglog(Period, umax ,'-k','LineWidth',lw); hold;
ylabel('Spectral displacement (cm)','FontSize',fsz,'FontWeight','normal','FontAngle','normal','fontname','times');
title( sprintf( 'Spectal displacement')) ;
axis square; xlim([min(Period) max(Period)*1.25]);  ylim([min(umax)/10 max(umax)*10]);
box on; grid on;


axes( 'Position', [0, 0.90 , 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
title_text = text( 0.5, 0.1 , 'SPECTRAL AMPLITUDES', 'FontSize', fsz +2 , 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;


hFig = gcf() ;
set ( hFig, "PaperType", "a4" ) ;
print( hFig, strcat(fpath, fname) , "-landscape" , "-fillpage", "-color", "-dpdf") ;
print( hFig, strcat(fpath, fname) , "-landscape" , "-color", "-djpg") ;

