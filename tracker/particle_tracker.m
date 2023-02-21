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

  [ beam_start_bin , beam_start_count , beam_main_bin , beam_main_count ] = set_histogram_beam( beam_start , beam_main , nb_bin , bin_step ) ;

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
    xlabel ( label_x_axis ); ylabel( label_y_axis ) ;
    axis("label","tic") ; 
    axis([-x_boundary x_boundary -y_boundary y_boundary]) ;


  subplot(2,2, 2 )  ;  hold on ;
  grid on ; box on ;
  set(gca, 'GridLineStyle', ':' ) ;
  set(gca,'fontsize',11,'fontname','times');
  %
  plot([-x_boundary x_boundary ],[0 0] , '--', "color",[0 0 0] + 0.5 ) ;  % Zero horizontal dashed grey (0.7)line
  plot([0 0],[-y_boundary y_boundary ] , '--', "color",[0 0 0] + 0.5 ) ;  %      vertical
    plot( beam_main( : , 1  ) , beam_main( : , 2 ) , ".","color", colour_name );
    xlabel ( label_x_axis ); ylabel( label_y_axis ) ;
    axis("label","tic") ; 
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
    xlabel ( label_x_axis ) ; ylabel("particle count") ;
    axis("label","tic") ;
    axis([-x_boundary x_boundary 0 count_boundary]) ;
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CALCULATE HISTOGRAM STEP
function [ dx ] = histogram_step(v1,v2,nb_bin);
%
%generate histogram from 'vec' into 'bin' bins
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
  x1 = min(v1(:,1)) ; x2 = max(v1(:,1)) ;  % original    beam
  y1 = min(v2(:,1)) ; y2 = max(v2(:,1)) ;  % transformed beam
  x_min = min ( x1 , y1 ) ;
  x_max = max ( x2 , y2 ) ;

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

  fh=fopen (fn,"r");
  [m,n,l]=scanning(fh);
  em0=[abs(m(1,1)),abs(m(2,2)),abs(m(3,3)),abs(m(4,4)),0,delta];
  %
  t0 = time ;
  p0=origine(em0,npart);
  % printf("Generated %i particles in %fs\n", npart, time-t0);
  %
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


