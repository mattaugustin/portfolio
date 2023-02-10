function plotSpectra_all(Period,amax,vmax, umax, fpath,fname)
% Plot response spectra
close all ;
figure() ; clf ;
% clf ;
fsz = 16;
lw = 2;
set(gcf,'position', [50 200 1300 450]) ;

%Period = period ;
%amax = sp_acc ;
%vmax = sp_vel ;
% umax = sp_dis ;

subplot_special(3,1, fsz ) ; hold;
% set(gca,'FontSize',fsz);
loglog(Period, amax ,'-k','LineWidth',lw); hold;
% xlabel('Period (s)','FontSize',fsz,'FontWeight','normal','FontAngle','normal','fontname','times');
ylabel('Pseudo spectral acceleration (cm/s^2)','FontSize',fsz,'FontWeight','normal','FontAngle','normal','fontname','times');
title( sprintf( 'Spectal acceleration')) ;
axis ("square") ; xlim([min(Period) max(Period)*1.25]); ylim([min(amax)/10 max(amax)*10]);
box on; grid on;
set(gca,'FontSize',fsz);
set(gca,'Xtick',[0.01 0.1 1 10 100]);
set(gca,'XTickLabel',{'0.01s','0.1s','1s','10s','100s'},'fontsize',fsz);
set(gca,'fontsize',[fsz-1],'fontname','times');


% subplot(1,3,2) ; hold;
subplot_special(3,2, fsz ) ; hold ;
loglog(Period, vmax ,'-k','LineWidth',lw); hold;
ylabel('Spectral velocity (cm/s)','FontSize',fsz,'FontWeight','normal','FontAngle','normal','fontname','times');
title( sprintf( 'Spectal velocity')) ;
axis square; xlim([min(Period) max(Period)*1.25]);  ylim([min(vmax)/10 max(vmax)*10]);
box on; grid on;
set(gca,'FontSize',fsz);
set(gca,'Xtick',[0.01 0.1 1 10 100]);
set(gca,'XTickLabel',{'0.01s','0.1s','1s','10s','100s'},'fontsize',fsz);
set(gca,'fontsize',[fsz-1],'fontname','times');


% subplot(1,3,3);  hold;
subplot_special(3 , 3 , fsz ) ; hold ;
loglog(Period, umax ,'-k','LineWidth',lw); hold;
ylabel('Spectral displacement (cm)','FontSize',fsz,'FontWeight','normal','FontAngle','normal','fontname','times');
title( sprintf( 'Spectal displacement')) ;
axis square; xlim([min(Period) max(Period)*1.25]);  ylim([min(umax)/10 max(umax)*10]);
box on; grid on;
set(gca,'FontSize',fsz);
set(gca,'Xtick',[0.01 0.1 1 10 100]);
set(gca,'XTickLabel',{'0.01s','0.1s','1s','10s','100s'},'fontsize',fsz);
set(gca,'fontsize',[fsz-1],'fontname','times');


axes( 'Position', [0, 0.90 , 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
title_text = text( 0.5, 0.1 , 'SPECTRAL AMPLITUDES', 'FontSize', fsz +2 , 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

% delete(title_text)

% fig = gcf;
% fig = gca;
% mkdir fpath
% saveas ( fig , "./plots/figure1.png") ;
% saveas ( fig , "./plots/figure1.pdf") ;
% saveas ( fig , strcat( fpath , fname , ".pdf" )) ;
% saveas ( fig , "martina/martina_again.pdf" ) ;

% print ( fig , "martina/test.pdf" , "-dpdfcairo",  '-S1000,600');
% print( fig, "./martina/figure1.png" , "png") ;

hFig = gcf() ;
set ( hFig, "PaperType", "a4" ) ;
% set ( hFig, "PaperOrientation", "landscape" ) ;
% set ( hFig, "PaperPosition", "bestfit" ) ;
% print( hFig, "./martina/martina_again2.pdf" , "-landscape" , "-fillpage", "-color", "-dpdf") ;
print( hFig, strcat(fpath, fname) , "-landscape" , "-fillpage", "-color", "-dpdf") ;
print( hFig, strcat(fpath, fname) , "-landscape" , "-color", "-djpg") ;
