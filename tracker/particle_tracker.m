close all ; clear all ; clc ;
pkg load statistics;

% READ MATRIX
function [m,n,i]=scanning(fid)
  i=0;
  while (1)
  [v,fn]=fscanf(fid,"%14f%14f%14f%14f%14f %1i%1i%1i%1i%1i%1i",11);
  if (feof (fid) )
    m=[]; n=[]; i=-1;
    break
  endif
  if (fn != 11)
    break
  endif
  i++;
  m(i,:)=v(1:5)';
  n(i,:)=v(6:11)';
  endwhile
  fscanf(fid,"%s",1);
  % printf("Read transformation matrix \n" );
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GENERATE ORIGINAL RANDOM-DISTRIBUTED PARTICLE BEAM
function p=origine(em0,ps);
  t0 = time ;
  for i=1:5
    p(:,i)=em0(i)*(betarnd(2,2,ps,1)*2-1);
   endfor
  p(:,6)=em0(6);%*round(betarnd(2,2,ps,1)*2-1);
  % printf("Generated %i particles in %fs\n", ps, time-t0);
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TRANSFORM PARTICLE BEAM USING TRANSFER MATRIX
function p=transform(p0,ps,m,n,l);
%transforms rays accoridng to transfer matrix
  t0 = time ;
  p=zeros(ps,5);
  v0=ones(ps,l);
  for j=1:l
    for i=1:6
      v0(:,j)=v0(:,j).*(p0(:,i).^n(j,i));
    endfor
  endfor
    for j=1:l
    for i=1:5
      p(:,i)+=m(j,i)*v0(:,j);
    endfor
  endfor
  % printf("Transformed %i particles from original beam in %fs\n", ps, time-t0);
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE FULL WIDTH HALF MAXIMUM FOR 1 MASS BEAM
function fwhm=histo_origin(vec,bin,pn);
%generate histogram from 'vec' into 'bin' bins
  x1=min(vec);
  x2=max(vec);

% set binsizes
  dx=(x2-x1)/bin;
  dxexp=floor(log10(dx)); % exponent    (b in n=a*10^b)
  dxsig=dx*10^-dxexp;     % significant (a in n=a*10^b)
  sig=[1,2,5]; % allowed factors for binsizes
  dx=max(sig(find(dxsig>=sig)))*10^dxexp;

% convert vec into bin indices
  veci=(vec-(x1+x2)/2)/dx;
  veci=floor(veci);
  veci0=1-min(veci);
  veci+=veci0;
  binn=max(veci);

% generate bin coordinates (LH side of bin)
  binx=1:binn;
  binx-=veci0;
  binx*=dx;
  binx+=(x1+x2)/2;

% fill bins
  for i=1:binn
    biny(i)=sum(veci==i);
  endfor

% find the FWHM
  hmax=max(biny)/2;
  tmp=find(biny>1.2*hmax);
  xl1=binx(tmp(1)); xr2=binx(tmp(length(tmp)));
  tmp=find(biny>0.8*hmax);
  xl2=binx(tmp(1)); xr1=binx(tmp(length(tmp)));
  xl=(xl1+xl2+dx)/2; xr=(xr1+xr2+dx)/2;
  fwhm=xl-xr;

% plot profile
  if (pn!=-1)
    subplot (2,2,2+pn); hold on ;
  	grid on;
	set(gca, 'GridLineStyle', '-');
    grid(gca,'minor') ;
	[xs,xys]= stairs( [binx,binx(binn)+dx], [biny,biny(binn)] );
  % [zs,zys]=stairs([binz,binz(binzn)+dz], [binzy,binzy(binzn)]);
  plot(xs,xys,"k");
	axis("labelx","ticx");  % axis([-5 5]) ;
  axis("labely","ticy"); axis([-5 5 0 1e4]) ;
	xlabel("particle horizontal position x"); ylabel("particle count") ;
  endif
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE FULL WIDTH HALF MAXIMUM FOR 3 NEAR-MASS BEAMS
function [fwhm,fwhm2,fwhm3,sp1,sp2]=histo_ens(v1,v2,v3,bin,pn);
%generate histogram from 'vec' into 'bin' bins
  x1=min(v1); x2=max(v1);
  y1=min(v2); y2=max(v2);
  z1=min(v3); z2=max(v3);

% set binsizes
  dx=(x2-x1)/bin;
  dxexp=floor(log10(dx)); % exponent    (b in n=a*10^b)
  dxsig=dx*10^-dxexp;     % significant (a in n=a*10^b)
  sig=[1,2,5]; % allowed factors for binsizes
  dx=max(sig(find(dxsig>=sig)))*10^dxexp;
% meme ope pour y et z
  dy=(y2-y1)/bin;
  dyexp=floor(log10(dy));
  dysig=dy*10^-dyexp;
  sig=[1,2,5];
  dy=max(sig(find(dysig>=sig)))*10^dyexp;
%
  dz=(z2-z1)/bin;
  dzexp=floor(log10(dz));
  dzsig=dz*10^-dzexp;
  sig=[1,2,5];
  dz=max(sig(find(dzsig>=sig)))*10^dzexp;

% convert vec into bin indices
  veci=(v1-(x1+x2)/2)/dx;
  veci=floor(veci);
  veci0=1-min(veci);
  veci+=veci0;
  binxn=max(veci);
%
  vecj=(v2-(y1+y2)/2)/dy;
  vecj=floor(vecj);
  vecj0=1-min(vecj);
  vecj+=vecj0;
  binyn=max(vecj);
 %
  veck=(v3-(z1+z2)/2)/dz;
  veck=floor(veck);
  veck0=1-min(veck);
  veck+=veck0;
  binzn=max(veck);
%
% generate bin coordinates (LH side of bin)
  binx=1:binxn;
  binx-=veci0;
  binx*=dx;
  binx+=(x1+x2)/2;
%
  biny=1:binyn;
  biny-=vecj0;
  biny*=dy;
  biny+=(y1+y2)/2;
 %
  binz=1:binzn;
  binz-=veck0;
  binz*=dz;
  binz+=(z1+z2)/2;
%
  % fill bins
  for i=1:binxn
    binxy(i)=sum(veci==i);
  endfor

  for j=1:binyn
    binyy(j)=sum(vecj==j);
  endfor

  for k=1:binzn
    binzy(k)=sum(veck==k);
  endfor

% find the FWHM
  hmax=max(binxy)/2;
  tmp=find(binxy>1.2*hmax);
  xl1=binx(tmp(1)); xr2=binx(tmp(length(tmp)));
  tmp=find(binxy>0.8*hmax);
  xl2=binx(tmp(1)); xr1=binx(tmp(length(tmp)));
  xl=(xl1+xl2+dx)/2; xr=(xr1+xr2+dx)/2;
  fwhm=xl-xr;

  hmax2=max(binyy)/2;
  tmp2=find(binyy>1.2*hmax2);
  yl1=biny(tmp2(1)); yr2=biny(tmp2(length(tmp2)));
  tmp2=find(binyy>0.8*hmax2);
  yl2=biny(tmp2(1)); yr1=biny(tmp2(length(tmp2)));
  yl=(yl1+yl2+dy)/2; yr=(yr1+yr2+dy)/2;
  fwhm2=yl-yr;

  hmax3=max(binzy)/2;
  tmp3=find(binzy>1.2*hmax3);
  zl1=binz(tmp3(1)); zr2=binz(tmp3(length(tmp3)));
  tmp3=find(binzy>0.8*hmax3);
  zl2=binz(tmp3(1)); zr1=binz(tmp3(length(tmp3)));
  zl=(zl1+zl2+dz)/2; zr=(zr1+zr2+dz)/2;
  fwhm3=zl-zr;

  % spread between curves regarding fwhm
  sp1=abs(xl-yr) %spread between M- and Mo
  sp2=abs(zl-xr) %spread between M- and Mo

% plot profile
  if (pn!=-1)
    subplot (2,2,2+pn); hold on ;
	grid on;
	set(gca, 'GridLineStyle', '-');
    grid(gca,'minor') ;
	%stairs( [binx,binx(binxn)+dx], [binxy,binxy(binxn)],"k",[biny,biny(binyn)+dy], [binyy,binyy(binyn)],"b",[binz,binz(binzn)+dz], [binzy,binzy(binzn)],"r" );
	[xs,xys]=stairs([binx,binx(binxn)+dx], [binxy,binxy(binxn)]);
	[ys,yys]=stairs([biny,biny(binyn)+dy], [binyy,binyy(binyn)]);
	[zs,zys]=stairs([binz,binz(binzn)+dz], [binzy,binzy(binzn)]);
	plot(xs,xys,"k",ys,yys,"b",zs,zys,"r");
	axis("labelx","ticx"); % axis([-5 5]) ;
  axis("labely","ticy"); axis([-5 5 0 1e4]) ;
	xlabel("particle horizontal position x"); ylabel("particle count") ;
	% hold off ;
  endif
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAME_AXES
function b=sameaxes
  for i=1:2
    subplot (2,2,i);
    a(i,:)=axis;
  endfor

  b(2)=max(abs(vec( a(:, 1:2) )));
  b(1)=-b(2);
  b(4)=max(abs(vec( a(:, 3:4) )));
  b(3)=-b(4);
  % b=[-30,30,-30,30];
  subplot (2,2,1); axis(b);
  subplot (2,2,2); axis(b);
  subplot (2,2,3); axis([b(1),b(2)]);
  subplot (2,2,4); axis([b(1),b(2)]);
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOTTING_X
function plottingx(v0,v1,v2,v3,pn);
 % recuperer colonnes avec x et a
	l0_x = v0(:,1)  ; l0_a = v0(:,2) ;
  l1_x = v1(:,1)  ; l1_a = v1(:,2) ;
	l2_x = v2(:,1) ;  l2_a = v2(:,2) ;
  l3_x = v3(:,1) ;  l3_a = v3(:,2) ;
  %
	nb_part = rows(v0);
	subplot(2,2,pn);  hold on ; % grid on ;
  axis([-5 5 -20 20]) ;
	plot( l(1:floor(nb_part/3)) , l2(1:floor(nb_part/3)),".k;beam source;");                     % box off ;
  plot( l( ceil(nb_part/3):floor(2*nb_part/3) ) , l2( ceil(nb_part/3):floor(2*nb_part/3) ),".r;beam source;"); % box off ;
  plot( l( ceil(2*nb_part/3):nb_part     ) , l2( ceil(2*nb_part/3):nb_part     ),".b;beam source;"); % box off ;
  % box off ;
  xlabel ("particle horizontal position x"); ylabel("particle horizontal angle a") ;
	histo_origin(l,100,pn);

  nb_bin = 100 ;
  histogram_beam( v0, v1 , nb_bin , pn ) ;

	subplot(2,2,pn+1); hold on ;
  % axis([-5 5 -20 20]) ;
	% grid on ;
	hold on ; grid on ; axis([-5 5 -75 75]) ;
	% plot(l3(1:n),l4(1:n),".k;mass 100.05;",l5(1:n),l6(1:n),".b;mass 99.95;",l7(1:n),l8(1:n),".r;mass 100.05;");
	plot(l3(1:n),l4(1:n),".k",l5(1:n),l6(1:n),".b",l7(1:n),l8(1:n),".r");
	plot([-50,50],[0 0], 'k-') ;
	plot([0 0],[-100 100], 'k-') ;
	hold off ;
	hold on ; grid on ;
	[fwhm,fwhm2,fwhm3,sp1,sp2]=histo_ens(l3,l5,l7,100,(pn+1));
	printf("fwhm for Mo: %f \n" , abs(fwhm))  ;
	printf("fwhm2 for M-: %f \n", abs(fwhm2)) ;
	printf("fwhm3 for M+: %f \n", abs(fwhm3)) ;
	% printf("spread between M- and Mo: %f \n",abs(sp1)) ;
	% printf("spread between M+ and Mo: %f \n",abs(sp2)) ;
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOTTING_Y
function plottingy(v1,v2,v3,v4,pn);
 % recuperer colonnes avec y et b
	l=v1(:,3) ; l2=v1(:,4) ; l3=v2(:,3) ; l4=v2(:,4) ;
	l5=v3(:,3) ; l6=v3(:,4) ; l7=v4(:,3) ; l8=v4(:,4) ;
	n=rows(v1);
	subplot(2,2,pn);
	plot(l(1:n),l2(1:n),".k;beam source;");
	histo_origin(l,100,pn);
	subplot(2,2,pn+1);
	hold on ; grid on ; axis([-10 10]) ;
	plot(l3(1:n),l4(1:n),".k",l5(1:n),l6(1:n),".b",l7(1:n),l8(1:n),".r");
	plot([-30,30],[0 0], 'k-') ;
	plot([0 0],[-30 30], 'k-') ;
	hold off ;
	% grid on ;
	histo_ens(l3,l5,l7,100,(pn+1));
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CREATE HISTOGRAMS (BIN OPTIMIZED) FROM ORIGINAL AND TRANSFORMED BEAM
function [ binx,binxy , biny,binyy ] = histogram_beam(v0,v1,nb_bin);
%
%generate histogram from 'vec' into 'bin' bins
    %  v0 = p0 ; v1 = p ;   % (example)
    %  nb_bin = 100 ;       % (example)
  x1=min(v0(:,1)) ; x2=max(v0(:,1));  % original    beam
  y1=min(v1(:,1)) ; y2=max(v1(:,1));  % transformed beam
  x_min = min ( x1 , y1 ) ;
  x_max = max ( x2 , y2 ) ;

% set binsizes
  dx=(x_max - x_min )/ nb_bin;
  dxexp=floor(log10(dx)); % exponent    (b in n=a*10^b)
  dxsig=dx*10^-dxexp;     % significant (a in n=a*10^b)
  sig=[1,2,5];            % allowed factors for binsizes
  dx=max(sig(find(dxsig>=sig)))*10^dxexp;

% convert vec into bin indices for both original and transformed beams
  veci=(v0(:,1) - (x1+x2)/2)/dx;
  veci=floor(veci);
  veci0=1-min(veci);
  veci+=veci0;
  binxn=max(veci)+1;  % add 1 to ensure symetry later of bin coordinates at the right-side

  vecj=(v1(:,1)-(y1+y2)/2)/dx ;
  vecj=floor(vecj);
  vecj0=1-min(vecj);
  vecj+=vecj0;
  binyn=max(vecj) + 1;

% generate bin coordinates (LH side of bin)
  binx=1:(binxn);
  binx-=veci0;
  binx*=dx;
  binx+=(x1+x2)/2;              % adjust beam indices if beam not centered at x=0
  binx = round(binx*100)/100 ;  % round at 2 digits precision to ensure "round" values for bin coordinates <> -0.03
  ## instead of -0.030008
%
  biny=1:binyn;
  biny-=vecj0;
  biny*=dx;
  biny+=(y1+y2)/2;
  biny = round(biny*100)/100 ;

% fill bins
  for i=1:(binxn)
    binxy(i)=sum(veci==i);
  endfor

  for j=1:binyn
    binyy(j)=sum(vecj==j);
  endfor

endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOT HISTOGRAM (using histogram results)
function plot_single_beam( beam_start , beam_main  , component , nb_bin , bin_step , colour_name ) ;  % red, black, blue, green, magenta, cyan

 %  colour_name = "red" ;  % = 1 ;
 %  beam_start  = beam_start_ter(:,1:2) ;    beam_main  = beam_delta_plus(:,1:2) ;

  n_part = rows(beam_start);
  min_beams_x  = min( min(beam_start(:,1)) , min(beam_main(:,1)) ) ;
  max_beams_x  = max( max(beam_start(:,1)) , max(beam_main(:,1)) ) ;
  x_boundary   = max( abs(min_beams_x) , max_beams_x  ) *1.1  ;

  min_beams_y  = min( min(beam_start(:,2)) , min(beam_main(:,2)) ) ;
  max_beams_y  = max( max(beam_start(:,2)) , max(beam_main(:,2)) ) ;
  y_boundary   = max( abs(min_beams_y) , max_beams_y  ) ;

% Prepare histogram values/settings prior to plotting
  % [ beam_start_bin , beam_start_count , beam_main_bin , beam_main_count ] = histogram_beam( beam_start , beam_main , nb_bin) ;
  [ beam_start_bin , beam_start_count , beam_main_bin , beam_main_count ] = set_histogram_beam( beam_start , beam_main , nb_bin , bin_step ) ;

  %beam_start_bin = binx ;
  %beam_start_count = binxy ;
  %beam_main_bin = biny ;
  % beam_main_count = binyy ;

  % dx = beam_start_bin(2) - beam_start_bin(1) ;
  dx = bin_step ;
  count_boundary = max( max(beam_start_count) , max(beam_main_count) )*1.2 ;
  [histo_beam_start_bin, histo_beam_start_count] = stairs( [beam_start_bin(1)-dx , beam_start_bin , beam_start_bin(end)+dx ], [ 0 , beam_start_count , 0 ] );
  [histo_beam_main_bin, histo_beam_main_count]   = stairs( [beam_main_bin(1)-dx  , beam_main_bin , beam_main_bin(end)+dx ]  , [ 0 , beam_main_count , 0 ] );

  if ( strncmpi(component , "horizontal" , 3 )    )
      label_x_axis = "particle horizontal position x (mm)" ;   label_y_axis = ("particle horizontal angle a (mrad)") ;
  else ( strncmpi( component , "vertical" , 3)   )
      label_x_axis = "particle vertical position y (mm)"   ;   label_y_axis = ("particle vertical angle b (mrad)")   ;
  endif

% Plot beam at start and after transformation
  subplot(2,2, 1 )  ;  hold on ;
  grid on ; box on ;
  set(gca, 'GridLineStyle', ':' ) ;
  set(gca,'fontsize',11,'fontname','times');
  %
  % rdm_rows = randi(rows(beam_start),5000,1) ;
  plot([-x_boundary x_boundary ],[0 0] , '--', "color",[0 0 0] + 0.5 ) ;  % Zero horizontal dashed grey (0.7)line
  plot([0 0],[-y_boundary y_boundary ] , '--', "color",[0 0 0] + 0.5 ) ;  %      vertical
    plot( beam_start( : , 1  ) , beam_start( : , 2 ) , ".","color", colour_name );
    % plot( beam_start( rdm_rows , 1  ) , beam_start( rdm_rows , 2 ) , ".","color", colour_name );
    % xlabel ("particle horizontal position x (mm)"); ylabel("particle horizontal angle a (mrad)") ;
    xlabel ( label_x_axis ); ylabel( label_y_axis ) ;
    axis("label","tic") ; % axis("labely","ticy") ;
    axis([-x_boundary x_boundary -y_boundary y_boundary]) ;


  subplot(2,2, 2 )  ;  hold on ;
  grid on ; box on ;
  set(gca, 'GridLineStyle', ':' ) ;
  set(gca,'fontsize',11,'fontname','times');
  %
  % rdm_rows = randi(rows(beam_start),5000,1) ;
  plot([-x_boundary x_boundary ],[0 0] , '--', "color",[0 0 0] + 0.5 ) ;  % Zero horizontal dashed grey (0.7)line
  plot([0 0],[-y_boundary y_boundary ] , '--', "color",[0 0 0] + 0.5 ) ;  %      vertical
    plot( beam_main( : , 1  ) , beam_main( : , 2 ) , ".","color", colour_name );
    % plot( beam_start( rdm_rows , 1  ) , beam_start( rdm_rows , 2 ) , ".","color", colour_name );
    % xlabel ("particle horizontal position x (mm)"); ylabel("particle horizontal angle a (mrad)") ;
    xlabel ( label_x_axis ); ylabel( label_y_axis ) ;
    axis("label","tic") ; % axis("labely","ticy") ;
    axis([-x_boundary x_boundary -y_boundary y_boundary]) ;


% Plot histograms of beams at start and after transformation
  subplot(2,2, 3  )  ;  hold on ;   % grid(gca,'minor') ;
  grid on ; box on ;
  set(gca, 'GridLineStyle', ':' ) ;
  set(gca,'fontsize',11,'fontname','times');
  %
  plot([-x_boundary x_boundary ],[0 0] , '--', "color",[0 0 0] + 0.7 ) ;  % Zero horizontal dashed grey (0.7)line
  plot([0 0],[-y_boundary y_boundary ] , '--', "color",[0 0 0] + 0.7 ) ;  %      vertical
  plot( histo_beam_start_bin , histo_beam_start_count ,"color", colour_name );
    % xlabel ("particle horizontal position x (mm)"); ylabel("particle count") ;
    xlabel ( label_x_axis ); ylabel("particle count") ;
    axis("label","tic") ;
    axis([-x_boundary x_boundary 0 count_boundary]) ;

  subplot(2,2, 4  )  ;  hold on ;
  grid on ; box on ;
  set(gca, 'GridLineStyle', ':' ) ;
  set(gca, 'fontsize',11,'fontname','times');
  %
  plot([-x_boundary x_boundary ],[0 0] , '--', "color",[0 0 0] + 0.7 ) ;  % Zero horizontal dashed grey (0.7)line
  plot([0 0],[-y_boundary y_boundary ] , '--', "color",[0 0 0] + 0.7 ) ;  %      vertical
  plot( histo_beam_main_bin , histo_beam_main_count ,"color", colour_name );
    % xlabel ("particle horizontal position x (mm)"); ylabel("particle count") ;
    xlabel ( label_x_axis ) ; ylabel("particle count") ;
    axis("label","tic") ;
    axis([-x_boundary x_boundary 0 count_boundary]) ;
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CALCULATE HISTOGRAM STEP
function [ dx ] = histogram_step(v1,v2,nb_bin);
%
%generate histogram from 'vec' into 'bin' bins
    %  v0 = p0 ; v1 = p ;   % (example)
    %  nb_bin = 100 ;       % (example)
  x1=min(v1(:,1)) ; x2=max(v1(:,1));  % original    beam
  y1=min(v2(:,1)) ; y2=max(v2(:,1));  % transformed beam
  x_min = min ( x1 , y1 ) ;
  x_max = max ( x2 , y2 ) ;

% set binsizes
  dx=(x_max - x_min )/ nb_bin;
  dxexp=floor(log10(dx)); % exponent    (b in n=a*10^b)
  dxsig=dx*10^-dxexp;     % significant (a in n=a*10^b)
  % sig=[1,2,5];            % allowed factors for binsizes
  sig=[1,2,3,4,5];            % allowed factors for binsizes
  dx=max(sig(find(dxsig>=sig)))*10^dxexp;
%
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CREATE HISTOGRAMS (BIN OPTIMIZED) FROM ORIGINAL AND TRANSFORMED BEAM
function [ binx,binxy , biny,binyy ] = set_histogram_beam(v1,v2,nb_bin , dx );
%
%generate histogram from 'vec' into 'bin' bins
    %  v0 = p0 ; v1 = p ;   % (example)
    %  nb_bin = 100 ;       % (example)
  x1 = min(v1(:,1)) ; x2 = max(v1(:,1)) ;  % original    beam
  y1 = min(v2(:,1)) ; y2 = max(v2(:,1)) ;  % transformed beam
  x_min = min ( x1 , y1 ) ;
  x_max = max ( x2 , y2 ) ;

% set binsizes
  % dx=(x_max - x_min )/ nb_bin;
  % dxexp=floor(log10(dx)); % exponent    (b in n=a*10^b)
  % dxsig=dx*10^-dxexp;     % significant (a in n=a*10^b)
  % sig=[1,2,5];            % allowed factors for binsizes
  % dx=max(sig(find(dxsig>=sig)))*10^dxexp;

% convert vec into bin indices for both original and transformed beams
  veci=(v1(:,1) - (x1+x2)/2)/dx;
  veci=floor(veci);
  veci0=1-min(veci);
  veci+=veci0;
  binxn=max(veci)+1;  % add 1 to ensure symetry later of bin coordinates at the right-side

  vecj=(v2(:,1)-(y1+y2)/2)/dx ;
  vecj=floor(vecj);
  vecj0=1-min(vecj);
  vecj+=vecj0;
  binyn=max(vecj) + 1;

% generate bin coordinates (LH side of bin)
  binx=1:(binxn);
  binx-=veci0;
  binx*=dx;
  binx+=(x1+x2)/2;              % adjust beam indices if beam not centered at x=0
  binx = round(binx*100)/100 ;  % round at 2 digits precision to ensure "round" values for bin coordinates <> -0.03
  ## instead of -0.030008
%
  biny=1:binyn;
  biny-=vecj0;
  biny*=dx;
  biny+=(y1+y2)/2;
  biny = round(biny*100)/100 ;

% fill bins
  for i=1:(binxn)
    binxy(i)=sum(veci==i);
  endfor

  for j=1:binyn
    binyy(j)=sum(vecj==j);
  endfor

endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CALCULATE BEAM STATISTICS
function [v_min,v_max,mean_v,v_sigma]=stats(v,npart)
	mean_v=0 ;
	for i=1:npart
		mean_v+=v(i,1)/npart;
	endfor
	v_min=min(v(:,1)) ;
	v_max=max(v(:,1)) ;
	v_mean=v(:,1)-mean_v*ones(rows(v),1) ;
	v_mean_sq=v_mean.*v_mean ;
	v_sigma =sqrt(sum(v_mean_sq)/npart) ;
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% READ MATRIX, GENERATE & TRANSFORM BEAM
function [p1,p0] = particles(fn,npart,nplane,delta);

  % printf( "Analysing output from <%s>\n", fn );

  fh=fopen (fn,"r");
  [m,n,l]=scanning(fh);
  em0=[abs(m(1,1)),abs(m(2,2)),abs(m(3,3)),abs(m(4,4)),0,delta];
  %
  t0 = time ;
  p0=origine(em0,npart);
  % printf("Generated %i particles in %fs\n", npart, time-t0);
  %
  % p0(:,5)=zeros(npart,1);
  % p0(:,6)=0.005*ones(npart,1);
  [m,n,l]=scanning(fh);
    pt=0;
    while (pt<nplane)
      pt++;
      [m,n,l]=scanning(fh);
      [m,n,l]=scanning(fh);
	    if (l==-1) break; endif
      %
      p1=transform(p0,npart,m,n,l);
      %
   endwhile
  fclose (fh);
  printf("Generated and transformed %i particles in %fs\n", npart, time-t0);

endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      MAIN  ROUTINE    -  MAIN  ROUTINE    - MAIN  ROUTINE    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dn='C:/Octave-7.3.0/cosy/mine/';               % directory
fn='output_ms120_mc90_wienfil.txt';

printf( "***************************\n");
printf( "** PARTICLE BEAM TRACKER **\n");
printf( "***************************\n\n\n");

printf( "Analysing output from <%s>\n\n", fn );

npart  = 5e4 ;
nplane = 7 ;

delta = 0 ;
delta1= -0.0001 ;
delta2=  0.0001 ;

time_start = time ;
[beam_delta, beam_origin]           =  particles( [dn,fn],npart,nplane,delta  ) ;
[beam_delta_minus , beam_start_bis] =  particles( [dn,fn],npart,nplane,delta1 ) ;
[beam_delta_plus  , beam_start_ter] =  particles( [dn,fn],npart,nplane,delta2 ) ;
printf("\n => All beams of %i particles ready in %fs\n", npart, time-time_start) ;

% CONVERT FROM M to MM and FROM RAD to MRAD
beam_origin*=1e3 ; beam_start_bis*=1e3   ; beam_start_ter*=1e3  ;
beam_delta*=1e3  ; beam_delta_minus*=1e3 ; beam_delta_plus*=1e3 ;


printf("Now preparing beam plots ... \n\n") ;
printf("horizontal plane on figure (1) ... \n") ;
printf("vertical plane on figure (2) ... \n")   ;
% close all ;
nb_bin = 100 ;


[ bin_step ] = histogram_step( beam_origin(:,1:2)    , beam_delta(:,1:2)       , nb_bin );
[ temp ]     = histogram_step( beam_start_bis(:,1:2) , beam_delta_minus(:,1:2) , nb_bin );
bin_step     = max( bin_step , temp ) ;
[ temp ]     = histogram_step( beam_start_ter(:,1:2) , beam_delta_plus(:,1:2) , nb_bin );
bin_step     = max( bin_step , temp ) ;


figure(1);  % Select and plot on horizontal plane
plot_single_beam( beam_start_bis(:,1:2), beam_delta_minus(:,1:2) , "horizontal" ,nb_bin , bin_step , "blue" );
plot_single_beam( beam_start_ter(:,1:2), beam_delta_plus(:,1:2)  , "horizontal" ,nb_bin , bin_step , "red" );
plot_single_beam( beam_origin(:,1:2)   , beam_delta(:,1:2)       , "horizontal" ,nb_bin , bin_step , "black" );

% print('mine/particle_tracker_horiplane.png','-dpng','-S1000,750');

[ bin_step ] = histogram_step( beam_origin(:,3:4)    , beam_delta(:,3:4)       , nb_bin );
[ temp ]     = histogram_step( beam_start_bis(:,3:4) , beam_delta_minus(:,3:4) , nb_bin );
bin_step     = max( bin_step , temp ) ;
[ temp ]     = histogram_step( beam_start_ter(:,3:4) , beam_delta_plus(:,3:4) , nb_bin );
bin_step     = max( bin_step , temp ) ;

figure(2);  % Select and plot on vertical plane
plot_single_beam( beam_start_bis(:,3:4), beam_delta_minus(:,3:4) , "vertical" ,nb_bin , bin_step , "blue" );
plot_single_beam( beam_start_ter(:,3:4), beam_delta_plus(:,3:4)  , "vertical" ,nb_bin , bin_step , "red" );
plot_single_beam( beam_origin(:,3:4)   , beam_delta(:,3:4)       , "vertical" ,nb_bin , bin_step , "black" );



% print('mine/particle_tracker_vertiplane.png','-dpng','-S1000,750');


