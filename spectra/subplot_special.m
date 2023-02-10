function subplot_special(nb_plots , id_plot , fsz )

% - Build figure.
 % clf ;
 %set( gcf, 'Color', 'White', 'Unit', 'Normalized' , ...
 %   'Position', [0.1,0.1,0.6,0.6] ) ;
 % set(gca,'FontSize',fsz);
 % set(gca,'Xtick',[0.01 0.1 1 10 100]);
 % set(gca,'XTickLabel',{'0.01s','0.1s','1s','10s','100s'},'fontsize',fsz);
 % set(gca,'fontsize',[fsz-1],'fontname','times');
 
% nb_plots = 3 ;
% id_plot = 1 ;
    
% - Compute #rows/cols, dimensions, and positions of lower-left corners.
nCol = nb_plots ;  nRow = 1 ;
 rowH = 0.58 / nRow ;  colW = 0.7 / nCol ;
 colX = 0.09 + linspace( 0, 0.96, nCol+1 )  ;  colX = colX(1:end-1) ;
 rowY = 0.15  + linspace( 0.9, 0, nRow+1 )  ;  rowY = rowY(2:end) ;

 
% - Build subplots axes and plot data.
% dId = id_plot ;
rowId = 1   ;
colId = id_plot ;

axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
% plot( data(:,1), data(:,colId +2), 'b' ) ;
grid on ;
% xlabel( 'Period (s)' ) ;  
xlabel('Period (s)','FontSize',fsz,'FontWeight','normal','FontAngle','normal','fontname','times');
% set(gca,'FontSize',fsz);

% set(gca,'Xtick',[0.01 0.1 1 10 100]);
% set(gca,'XTickLabel',{'0.01s','0.1s','1s','10s','100s'},'fontsize',fsz);
% set(gca,'fontsize',[fsz-1],'fontname','times');
