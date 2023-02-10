% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                                                 % %
% %   CALCULATION OF PSEUDO SPECTRAL ACCELERATION   % %
% %                  PSEUDO SPECTRAL VELOCITY       % %
% %                  PSEUDO SPECTRAL DISPLACEMENT   % %
% %                                                 % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % %


close all
clear
clc

% Load acceleration file and prepare settings
    filename = sprintf('javak_ew_accel.txt');
    savename = sprintf('results/spectra_javak_ew.txt');
    # data_a=load ("-ascii", filename);
    printf("\n\n Reading text file ... \n\n");
    [a, b] =  textread( filename , "%f %f" );
    printf("\n\n Displaying first rows of text file and 'cleaning' ... \n");
    a(1:5)
    b(1:5)
    a(1,:)=[] ;
    b(1,:)=[] ;
    a(1:5)
    b(1:5)


% Prepare accelerogram data and associated parameters (time step, nb of datapoints...)
    acceleration = b  ;        %  convert here to m/s2 from (cm/mm/?) WHEN required
    dt=a(2,:)         ;        %  time step
    printf( "\n The timestep is %f s", dt ) ;
    nStep = numel(a);
    printf( "\n Number of datapoints is %d", nStep ) ;
    time  = a         ;
    % time  = dt:dt:endt;
    printf( "\n Signal duration is %f s", time(end) ) ;  %time(end)                 % signal duration

% Define set of periods for computation of response spectra (choose a desired set below)
    printf( "\n Set of structural periods (in s) for calculation of spectral values is \n") ;
    %period = [ 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00 2.00 3.00]
    period = logspace( -2 , 1 , 100 )      % p = logspace(minPower,maxPower,nPeriod);
    % period = 0.01:0.01:10 ;

% build response spectra (input: signal, signal duration and number of signal datapoints)
    printf( "\n Computation of spectral values starting ... ...  " ) ;
   [period , f, umax, vmax, amax] = respSpectra_period(acceleration , time(end), nStep , period );
    period=transpose(period);
    f=transpose(f);
    amax=100*transpose(amax);
    vmax=100*transpose(vmax);
    umax=100*transpose(umax);
    % tamax= [period f amax];       % save in csv format as period, frequency, value rsp spectra
    tamax= [period f amax vmax umax];

%% Save results in subdirectory
    [status, message, message_id] = mkdir("results");
    if status ~= 1
    disp("Error creating directory: " + message);
    endif

   printf("\n\n          Computation done. Now saving spectral values in %s\n" , savename);
   dlmwrite( savename , tamax);
   if exist(savename, "file") == 2
      disp("File successfully written!");
    else
      disp("Error writing file!");
    endif


%% Prepare for plotting (basic)
printf( "\n\n Prepare for plotting ... " ) ;
data = dlmread(  savename );
period=data(:,1);
frequency=data(:,2);
sp_acc=data(:,3);
sp_vel=data(:,4);
sp_dis=data(:,5);


%% Plot the response spectra (basic)
printf( "\n\n Prepare for basic plotting ... " ) ;
figure(1)
plot(period,sp_acc)
% semilogx( period, sp_acc)   % alternative with x-axis in logscale
xlabel('Period [s]')
ylabel('Pseudo Spectral acceleration [cm/s2]')
xlim([min(period) max(period)])


%% Plot the response spectra (enhanced)
printf( "\n\n Prepare for enhanced plotting ... " ) ;
mkdir plots
plotSpectra_all(  period,sp_acc, sp_vel , sp_dis , 'plots/','response_spectra'  )


