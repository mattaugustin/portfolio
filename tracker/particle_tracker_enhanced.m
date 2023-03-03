close all ; clear all ; clc ;
pkg load statistics;


function [nb_plane] = plane_scanning( fn )

  fid = fopen (fn,"r");
  % nb_plane = plane_scanning(fh);
  % nb_plane = nb_plane / 2 ;

  nb_plane = 0 ;
% Read the file line by line
while ~feof(fid)
    line = fgetl(fid);

    % Check if the line contains only "-"
    if strfind(line, " ----"  )
        nb_plane = nb_plane + 1;
    endif
endwhile
%
fclose(fid);
nb_plane = nb_plane / 2 ;

endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% CALCULATE HISTOGRAM STEP
function [ dx ] = histogram_step_enhanced( beam , component , nb_bin) ;

%  Prepare settings based on "vertical" or "horizontal"
  if ( strncmpi(component , "horizontal" , 3 )    )
      comp_col = 1 ;
  elseif ( strncmpi( component , "vertical" , 3)   )
      comp_col = 3 ;
  endif

%generate histogram from 'vec' into 'bin' bins

  nb_plane = length(beam) ;

  x1 = min( beam{1}(:,comp_col) )        ; x2 = max( beam{1}(:,comp_col) );         % original    beam
  y1 = min( beam{nb_plane}(:,comp_col) ) ; y2 = max( beam{nb_plane}(:,comp_col) );  % transformed beam
  x_min = min ( x1 , y1 ) ;
  x_max = max ( x2 , y2 ) ;

% set binsizes
  dx=(x_max - x_min )/ nb_bin;
  dxexp=floor(log10(dx));   % exponent    (b in n=a*10^b)
  dxsig=dx*10^-dxexp;       % significant (a in n=a*10^b)
  sig=[1,2,3,4,5];          % allowed factors for binsizes
  dx=max(sig(find(dxsig>=sig)))*10^dxexp;
%
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CREATE HISTOGRAMS (BIN OPTIMIZED) FROM ORIGINAL AND TRANSFORMED BEAM
function [stairs_histo_beam] = set_histogram_beam_enhanced( beam , component , nb_plane , nb_bin , dx );
%
%generate histogram from 'vec' into 'bin' bins
  % x1 = min(beam{1}(:,1))        ; x2 = max(beam{1}(:,1)) ;         % original    beam
  % y1 = min(beam{nb_plane}(:,1)) ; y2 = max(beam{nb_plane}(:,1)) ;  % final transformed beam
  % x_min = min ( x1 , y1 ) ;
  % x_max = max ( x2 , y2 ) ;

  if ( strncmpi(component , "horizontal" , 3 )    )
      comp_col = 1 ;
  else ( strncmpi( component , "vertical" , 3)   )
      comp_col = 3 ;
  endif

for beam_plane=1:nb_plane

  x1 = min(beam{beam_plane}(:,comp_col))  ;  x2 = max(beam{beam_plane}(:,comp_col)) ;         % original    beam

  % convert vec into bin indices for for beams in all planes
  veci  = (beam{beam_plane}(:,comp_col) - (x1+x2)/2)/dx;
  veci  = floor(veci)  ;
  veci0 = 1 - min(veci);
  veci+=veci0 ;
  binxn = max(veci) + 1  ;     % add 1 to ensure symetry later of bin coordinates at the right-side

  % generate bin coordinates (LH side of bin)
  binx = 1:(binxn);
  binx-=veci0;      % binx = binx - veci0 ;
  binx = binx*dx ;
  binx+=(x1+x2)/2;  % binx = binx + (x1+x2)/2;    % adjust beam indices if beam not centered at x=0
  binx = round(binx*100)/100 ;                    % round at 2 digits precision to ensure "round" values for bin coordinates <> -0.03
  ## instead of -0.030008

  % fill bins
  for i=1:(binxn)
    biny(i) = sum( veci==i );
  endfor
  %
  histo_beam{beam_plane} = [binx ; biny] ;
  %
  clear binx ; clear biny ;
endfor
%

% Create the histogram staircase
for beam_plane=1:nb_plane
      x = [histo_beam{beam_plane}(1,1)-dx , histo_beam{beam_plane}(1,:) , histo_beam{beam_plane}(1,end)+dx ] ;
      y = [ 0 , histo_beam{beam_plane}(2,:) , 0 ] ;
      %xy = stairs( x , y ) ;
      %stairs_histo_beam{beam_plane} = xy ;
      [x_stairs , y_stairs]         =  stairs ( x , y ) ;
      stairs_histo_beam{beam_plane} = [x_stairs , y_stairs] ;
     % clear x ; clear y ; clear x_stairs ; clear y_stairs ;
endfor
%

endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PLOT BEAMS & HISTOGRAMS FROM ORIGINAL AND TRANSFORMED BEAM
function plot_selected_beam_enhanced( beam , component , beam_plane , nb_bin , colour_name ) ;  % red, black, blue, green, magenta, cyan

% Prepare axis labels
  if ( strncmpi(component , "horizontal" , 3 )    )
      label_x_axis = "particle horizontal position x (mm)" ;   label_y_axis = ("particle horizontal angle a (mrad)") ;
      comp_col = 1 ;
  elseif ( strncmpi( component , "vertical" , 3)   )
      label_x_axis = "particle vertical position y (mm)"   ;   label_y_axis = ("particle vertical angle b (mrad)")   ;
       comp_col = 3 ;
  endif
%

% Prepare axis boundaries
  %beam = beam_delta  ;
  %component = "horizontal" ;
  % beam_plane = 5 ;
  %beam_plane = [2 4] ;
  % nb_part  = 50000 ;
  % dx = bin_step ;
  %colour_name = "black" ;
  % n_part = rows( beam{1} );

  % index_plane = 2 ;


  % Calculate bin_step for each plane
  % calculate boundaries on x , y and count number on each plane
  % recalculate histogram for each plane
  % plot for each plane

for index_plane=1:length(beam_plane)

  x_min  = min(  beam{beam_plane(index_plane)}(:,comp_col)  ) ;
  x_max  = max(  beam{beam_plane(index_plane)}(:,comp_col)  ) ;
  x_boundary(index_plane)   = max( abs(x_min)  , x_max  ) *1.1  ;

  min_beams_y  = min(  beam{beam_plane(index_plane)}(:,comp_col+1)  ) ;
  max_beams_y  = max(  beam{beam_plane(index_plane)}(:,comp_col+1)  ) ;
  y_boundary(index_plane)   = max( abs(min_beams_y)  , max_beams_y  ) *1.1  ;
%

% Calculate bin_step for plot horizontal axis
  % x_min = min( beam{beam_plane(index_plane)}(:,comp_col) )  ;  %  selected beam
  % x_max = max( beam{beam_plane(index_plane)}(:,comp_col) )  ;  %  selected beam

  dx=(x_max - x_min )/ nb_bin;
  dxexp=floor(log10(dx));   % exponent    (b in n=a*10^b)
  dxsig=dx*10^-dxexp;       % significant (a in n=a*10^b)
  sig=[1,2,3,4,5];          % allowed factors for binsizes
  dx=max(sig(find(dxsig>=sig)))*10^dxexp;


% Prepare histogram values/settings prior to plotting
  %  x1 = min_beams_x  ;  % min(beam{beam_plane(index_plane)}(:,comp_col))  ;
  % x2 = max_beams_x  ;  % max(beam{beam_plane(index_plane)}(:,comp_col))  ;         % original    beam

% convert vec into bin indices for for beams in all planes
  veci  = (beam{beam_plane(index_plane)}(:,comp_col) - (x_min+x_max)/2)/dx;
  veci  = floor(veci)  ;
  veci0 = 1 - min(veci);
  veci+=veci0 ;
  binxn = max(veci) + 1  ;     % add 1 to ensure symetry later of bin coordinates at the right-side

  % generate bin coordinates (LH side of bin)
  binx = 1:(binxn);
  binx-=veci0;      % binx = binx - veci0 ;
  binx = binx*dx ;
  binx+=(x_min+x_max)/2;  % binx = binx + (x1+x2)/2;    % adjust beam indices if beam not centered at x=0
  binx = round(binx*100)/100 ;                    % round at 2 digits precision to ensure "round" values for bin coordinates <> -0.03
  ## instead of -0.030008

% fill bins
  for i=1:(binxn)
    biny(i) = sum( veci==i );
  endfor

  %histo_beam{index_plane} = [binx ; biny] ;
  histo_beam  = [binx ; biny] ;

  %x = [histo_beam{index_plane}(1,1)-dx , histo_beam{index_plane}(1,:) , histo_beam{index_plane}(1,end)+dx ] ;
  %y = [ 0 , histo_beam{index_plane}(2,:) , 0 ] ;
  x = [histo_beam(1,1)-dx , histo_beam(1,:) , histo_beam(1,end)+dx ] ;
  y = [ 0 , histo_beam(2,:) , 0 ] ;

  [x_stairs , y_stairs]  =   stairs ( x , y )     ;
  stairs_histo_beam{index_plane} = [x_stairs , y_stairs] ;

  clear x ; clear y ; clear x_stairs ; clear y_stairs ; clear binx ; clear biny ; clear x_min ; clear x_max ; clear histo_beam ;
  % clear count_boundary ; clear histo_beam ; clear x_boundary ; clear y_boundary ;

  count_boundary(index_plane) =  max( stairs_histo_beam{index_plane}(:,2) )*1.2  ;

endfor
%%%

% x_boundary
% y_boundary
% count_boundary
% stairs_histo_beam{1}


%%  Plot beams in all plane !

  %figure(1) ; pause (0.1) ;
  %clf ; set(gcf,'position', [20 100 300*col_nb 450]) ;
  %fsz = 12 ;
  %lw = 2;
  col_nb = length( beam_plane ) ;


for index_plane=1:col_nb

  %  subplot_special( nb_plane, 1 , beam_plane , fsz ) ; hold on ;                 % subplot(nb_plane, rowId (1 if beam, 2 if histo) , plotId , fontsize )
  subplot( 2 , col_nb , index_plane ) ; hold on ;
  plot([-x_boundary(index_plane) x_boundary(index_plane) ],[0 0] , '--', "color",[0 0 0] + 0.5 ) ;        % Zero horizontal dashed grey (0.7)line
  plot([0 0],[-y_boundary(index_plane) y_boundary(index_plane) ] , '--', "color",[0 0 0] + 0.5 ) ;        %      vertical
  plot( beam{beam_plane(index_plane)}(:,comp_col) , beam{beam_plane(index_plane)}(:,comp_col+1) , ".","color", colour_name );

  xlabel( label_x_axis ) ; ylabel( label_y_axis ) ;
  title( sprintf( 'Beam in plane %i', beam_plane(index_plane) ) ) ;
  axis("label","tic") ; % axis("labely","ticy") ;
  axis( [-x_boundary(index_plane) x_boundary(index_plane) -y_boundary(index_plane) y_boundary(index_plane)] ) ;

% Plot corresponding histogram

  % subplot_special( nb_plane, 2 , beam_plane , fsz ) ; hold on ;                 %    subplot(nb_plane, rowId (1 if beam, 2 if histo) , plotId , fontsize )
  subplot( 2 , col_nb , index_plane + col_nb ) ; hold on ;
  plot([-x_boundary(index_plane) x_boundary(index_plane) ],[0 0] , '--', "color",[0 0 0] + 0.5 ) ;        %    Zero horizontal dashed grey (0.7)line
  plot([0 0],[-y_boundary(index_plane) y_boundary(index_plane) ] , '--', "color",[0 0 0] + 0.5 ) ;        %         vertical
  plot( stairs_histo_beam{index_plane}(:,1) , stairs_histo_beam{index_plane}(:,2) , "color" , colour_name );

  xlabel ( label_x_axis ) ; ylabel( "Particle count" ) ;
  title( ' {\fontsize{11}\it associated histogram} ') ;
  axis("label","tic") ; % axis("labely","ticy") ;
  axis([-x_boundary(index_plane) x_boundary(index_plane) 0 count_boundary(index_plane)]) ;

endfor
%

endfunction
%%%%%%%%%%%%%


% PLOT BEAMS & HISTOGRAMS FROM ORIGINAL AND TRANSFORMED BEAM
function plot_beam_enhanced( beam , component , nb_plane , nb_bin , bin_step , colour_name ) ;  % red, black, blue, green, magenta, cyan

% Prepare axis labels
  if ( strncmpi(component , "horizontal" , 3 )    )
      label_x_axis = "particle horizontal position x (mm)" ;   label_y_axis = ("particle horizontal angle a (mrad)") ;
      comp_col = 1 ;
  elseif( strncmpi( component , "vertical" , 3)   )
      label_x_axis = "particle vertical position y (mm)"   ;   label_y_axis = ("particle vertical angle b (mrad)")   ;
      comp_col = 3 ;
  endif
%

% Prepare axis boundaries
  % beam = beam_delta ;
  % nb_plane = 5 ;
  % n_part  = 5000 ;
  dx = bin_step ;
  % colour_name = "red" ;

  % n_part = rows( beam{1} );
  min_beams_x  = min( min(beam{1}(:,comp_col)) , min(beam{nb_plane}(:,comp_col)) ) ;
  max_beams_x  = max( max(beam{1}(:,comp_col)) , max(beam{nb_plane}(:,comp_col)) ) ;
  x_boundary   = max( abs(min_beams_x)  , max_beams_x  ) *1.1  ;

  min_beams_y  = min( min(beam{1}(:,comp_col+1)) , min(beam{nb_plane}(:,comp_col+1)) ) ;
  max_beams_y  = max( max(beam{1}(:,comp_col+1)) , max(beam{nb_plane}(:,comp_col+1)) ) ;
  y_boundary   = max( abs(min_beams_y)  , max_beams_y  ) *1.1  ;
%

% Prepare histogram values/settings prior to plotting

  [ stairs_histo_beam ] = set_histogram_beam_enhanced( beam , component , nb_plane , nb_bin , dx ) ;  % bin_step = dx
  count_boundary = max(  max(stairs_histo_beam{1}(:,2) ) , max( stairs_histo_beam{nb_plane}(:,2) )  )*1.2  ;


%%  Plot beams in all plane

  % figure(1) ;  % pause (0.1) ; % clf ;
  % set(gcf,'position', [20 100 1350 450]) ;
  fsz = 12 ;
  lw = 2;

for beam_plane=1:nb_plane

  subplot( 2 , nb_plane , beam_plane ) ; hold on ;                              %  subplot ( rows ,col , index )  read L to R , T to B
  plot([-x_boundary x_boundary ],[0 0] , '--', "color",[0 0 0] + 0.5 ) ;        %  Zero horizontal dashed grey (0.7)line
  plot([0 0],[-y_boundary y_boundary ] , '--', "color",[0 0 0] + 0.5 ) ;        %      vertical
  plot( beam{beam_plane}(:,comp_col) , beam{beam_plane}(:,comp_col+1) , ".","color", colour_name );

  xlabel( label_x_axis ) ; ylabel( label_y_axis ) ;
  title( sprintf( 'Beam in plane %i', beam_plane )) ;
  axis("label","tic") ; % axis("labely","ticy") ;
  axis( [-x_boundary x_boundary -y_boundary y_boundary] ) ;

% Plot corresponding histogram
%
  subplot( 2 , nb_plane , beam_plane + nb_plane ) ; hold on ;
  plot([-x_boundary x_boundary ],[0 0] , '--', "color",[0 0 0] + 0.5 ) ;        %    Zero horizontal dashed grey (0.7)line
  plot([0 0],[-y_boundary y_boundary ] , '--', "color",[0 0 0] + 0.5 ) ;        %         vertical
  plot( stairs_histo_beam{beam_plane}(:,1) , stairs_histo_beam{beam_plane}(:,2) , "color" , colour_name );

  xlabel ( label_x_axis ) ; ylabel( "Particle count" ) ;
  title( ' {\fontsize{11}\it associated histogram} ') ;
  axis("label","tic") ; % axis("labely","ticy") ;
  axis([-x_boundary x_boundary 0 count_boundary]) ;

endfor
%
endfunction
%%%%%%%%%%%%%


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
function [beam] = plane_particles( fn ,npart, nb_plane ,delta );

  % printf( "Analysing output from <%s>\n", fn );
  % fn = [dn,file_output]
  % fh=fopen (fn,"r");
  % nb_plane = plane_scanning(fh);
  % nb_plane = nb_plane / 2 ;
  % fclose(fh);
  %
  fh=fopen (fn,"r");
  [m,n,l]=scanning(fh);
  em0=[abs(m(1,1)),abs(m(2,2)),abs(m(3,3)),abs(m(4,4)),0,delta];
  %
  t0 = time ;
  p0=origine(em0,npart);
  % beam{pt} = origine(em0, npart) ;  %
  pt = 1 ;
  beam{pt} = p0 ;
  % printf("Generated %i particles in %fs\n", npart, time-t0);
  %
  [m,n,l]=scanning(fh);
  %
  pt++ ; % pt = 1 ;
    while ( pt <= nb_plane )
      % pt++;
      [m,n,l] = scanning(fh) ;
      [m,n,l] = scanning(fh) ;
	    if ( l==-1 ) break; endif
      %
      beam{pt} = transform( p0,npart,m,n,l );
      %
      pt++ ;
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

printf( "*********************************\n");
printf( "** PARTICLE BEAM TRACKER       **\n");
printf( "*********************************\n\n\n");

printf( "Analysing output from <%s>\n\n", fn );


% 1) DEFINE NB OF PARTICLES TO SIMULATE AND MASS DIFFERENCE W.R.T REFERENCE MASS (EX. M=100)
nb_part  = 50000 ;
delta = 0 ;
delta1= -0.0001 ;
delta2=  0.0001 ;

% 2) SIMULATE AND CALCULATE PROJECTED PARTICLES ACROSS ALL OBSERVATION PLANES
time_start = time ;
nb_plane = plane_scanning( [dn,fn] ) ;
printf("\n => The file contains %i observations planes \n\n", nb_plane ) ;
[beam_delta]       = plane_particles( [dn,fn],nb_part, nb_plane , delta  ) ;
[beam_delta_minus] = plane_particles( [dn,fn],nb_part, nb_plane , delta1 ) ;
[beam_delta_plus ] = plane_particles( [dn,fn],nb_part, nb_plane , delta2 ) ;
printf("\n => All beams of %i particles ready in %fs\n\n", nb_part, time-time_start) ;


% CONVERT FROM M to MM and FROM RAD to MRAD
for i=1:length(beam_delta)
    beam_delta{i} = beam_delta{i} * 1000 ;
endfor
for i=1:length(beam_delta_minus)
    beam_delta_minus{i} = beam_delta_minus{i} * 1000 ;
endfor
for i=1:length(beam_delta_plus)
    beam_delta_plus{i} = beam_delta_plus{i} * 1000 ;
endfor
%%%%%%%%%%%%%%%

printf("Now preparing beam plots ... \n\n") ;


% close all ;

nb_bin = 100 ;



% 3) HORIZONTAL PLANE PLOTS

printf("horizontal plane on figure (1) ... \n") ;

[ temp_delta ]       = histogram_step_enhanced( beam_delta       , "horizontal" , nb_bin )  ;
[ temp_delta_minus ] = histogram_step_enhanced( beam_delta_minus , "horizontal" , nb_bin ) ;
[ temp_delta_plus ]  = histogram_step_enhanced( beam_delta_plus  , "horizontal" , nb_bin )  ;
bin_step     = max( [temp_delta , temp_delta_minus , temp_delta_plus ] ) ;

figure(1);  % Select and plot on horizontal plane
pause (0.1) ; set(gcf,'position', [20 250 1350 400]) ;
plot_beam_enhanced( beam_delta       , "horizontal" , nb_plane , nb_bin , bin_step , "black" ); pause(0.3) ;
plot_beam_enhanced( beam_delta_plus  , "horizontal" , nb_plane , nb_bin , bin_step , "red"   ); pause(0.3) ;
plot_beam_enhanced( beam_delta_minus , "horizontal" , nb_plane , nb_bin , bin_step , "blue"  );

% print('mine/partracker_horiplane_start.png','-dpng','-S1000,750');



% 4)  VERTICAL PLANE PLOTS

printf("vertical plane on figure (2) ... \n")   ;

[ temp_delta ]       = histogram_step_enhanced( beam_delta       , "vertical" , nb_bin )  ;
[ temp_delta_minus ] = histogram_step_enhanced( beam_delta_minus , "vertical" , nb_bin ) ;
[ temp_delta_plus ]  = histogram_step_enhanced( beam_delta_plus  , "vertical" , nb_bin )  ;
bin_step     = max( [temp_delta , temp_delta_minus , temp_delta_plus ] ) ;


figure(2);  % Select and plot on vertical plane
pause (0.1) ; set(gcf,'position', [20 50 1350 400]) ;
plot_beam_enhanced( beam_delta       , "vertical" , nb_plane , nb_bin , bin_step , "black" ); pause(0.3) ;
plot_beam_enhanced( beam_delta_plus  , "vertical" , nb_plane , nb_bin , bin_step , "red"   ); pause(0.3) ;
plot_beam_enhanced( beam_delta_minus , "vertical" , nb_plane , nb_bin , bin_step , "blue"  );

% print('mine/partracker_vertiplane_start.png','-dpng','-S1000,750');



% 5)  INSPECT SPECIFIC PLANES, USING PLANE INDIVIDUAL AXIS SCALE (instead of axis scale common to all planes)

% Select beam of interest (ex: beam_delta) and select index # of specific planes to look at (ex. #2 and #4 )
beam_plane = [2 4] ;

printf("horizontal plane for user-selected beams on figure (3) ... \n" )   ;

figure(3); clf ; % Select and plot on vertical plane
pause (0.1) ; set(gcf,'position', [20 100 300*length(beam_plane) 400]) ;
plot_selected_beam_enhanced( beam_delta , "horizontal" , beam_plane , nb_bin  , "black" ) ;







#{
figure(1);  % Select and plot on horizontal plane
pause (0.1) ; set(gcf,'position', [20 100 1350 450]) ;
plot_single_beam_enhanced( beam_delta         , "horizontal" , nb_plane , nb_bin , bin_step , "black" ); hold ; % pause(0.3) ;
plot_single_beam_enhanced( beam_delta_plus  , "horizontal" , nb_plane , nb_bin , bin_step , "red"   ); hold on ; pause(0.3) ;
plot_single_beam_enhanced( beam_delta_minus , "horizontal" , nb_plane , nb_bin , bin_step , "blue"  );
#}



% beam_delta_selected = {beam_delta{2}, beam_delta{4}}   ;




