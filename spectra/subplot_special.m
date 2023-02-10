function subplot_special(nb_plots , id_plot , fsz )

% - Compute #rows/cols, dimensions, and positions of lower-left corners.
nCol = nb_plots ;  nRow = 1 ;
 rowH = 0.58 / nRow ;  colW = 0.7 / nCol ;
 colX = 0.09 + linspace( 0, 0.96, nCol+1 )  ;  colX = colX(1:end-1) ;
 rowY = 0.15  + linspace( 0.9, 0, nRow+1 )  ;  rowY = rowY(2:end) ;


% - Build subplots axes and add custom theme
rowId = 1   ;
colId = id_plot ;

axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
grid on ;
xlabel('Period (s)','FontSize',fsz,'FontWeight','normal','FontAngle','normal','fontname','times');
set(gca,'FontSize',fsz);
set(gca,'Xtick',[0.01 0.1 1 10 100]);
set(gca,'XTickLabel',{'0.01s','0.1s','1s','10s','100s'},'fontsize',fsz);
set(gca,'fontsize',[fsz-1],'fontname','times');


