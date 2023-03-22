
***In a nutshell :*** *I created this code to obtain engineering parameters from SAC earthquake signals (i.e. time histories) previously processed using a SAC (Seismic Analysis Code) code sequence shown in portfolio/SAC. This code reads and further processes SAC signals (clipping, metadata-matching...) before calculating the desired parameters using a set of functions while keeping track of the range of settings associated to the SAC signals (stations, geocoordinates...). I included several plotting options to clearly visualize the data and results. The final output is a dataframe containing seismic parameters derived from earthquake observations that can be compared to earthquake simulation data in a companion project section.*

A table containing information related to the earthquake events is provided to connect with the engineering parameters.
With this code is provided a set of SAC earthquake signals and associated station metadata referred to as PoleZero files (PZ).
 https://drive.google.com/drive/folders/1FRe2bp1G7aQ-0hjDYUV57vSaD_2b9zix?usp=sharing


### Background 
&nbsp; &nbsp; &nbsp; &nbsp; Previously (see portfolio/SAC), SAC signals were processed and saved in SAC format, and a set of files was available, in connection with horizontal component type (East-West or North-South), instrumentation channel (broadband or high broadband), filtering frequency (upper corner frequency), earthquake motion type (acceleration ...). Further processing of SAC files needs to be achieved in R code in order to ultimately obtain earthquake parameters of relevant engineering interest (PGA etc ...).    


### Task
Following the processing of RAW SAC signals in portfolio/SAC, a set of SAC files and associated PoleZero files are now available in the folder /observations, with the  sample below as an example of file organisation.
- folkestone/
- rasen/
- swansea/
   - Base_Bchan/
   - Base_Hchan/
   - Bchan_10Hz/
   - Bchan_25Hz/
   - Hchan_10Hz/
   - Hchan_35Hz/
     - accel
     - vel/
     - disp/
         - CWF.HHE
         - CWF.HHN
         - HPK.HHE
         - HPK.HHN
         - MCH1.HHE
         - MCH1.HHN
         - STNC.HHE
         - STNC.HHN
         - SWN1.HHE
         - SWN1.HHN
          
In this section, the goal is to further process the SAC signals and obtain several engineering parameters. These include Peak amplitudes (Peak Ground Velocity, Acceleration and Displacement), spectral accelerations at structural periods of interest, parameters related to Arias intensity. The final outcome is a table containing all desired parameters as well as information on the SAC signal parameters were derived from.


### Action

&nbsp; &nbsp; &nbsp; &nbsp; Calculating engineering parameters from the SAC files firstly requires to "prepare" said SAC files, a.k.a. ensuring these SAC files contain metadata (start time, coordinates, station name ...) consistent with PoleZero files and have a suitably optimized file size for later calculations. 

Thanks to the table "earthquake_info", I created and verified corresponding file paths, loaded the PZ and SAC files, and assembled a new table listing useful characteristics associated to each SAC file. These notably include earthquake name, the monitoring channel (B-chan or H-chan), the filtering corner frequency, seismicity recording station name, the horizontal component, the ground motion type (acceleration/velocity/displacement).     
Using the "eseis" package, it is easy to plot the SAC signal with certain associated metadata in caption (Figure 1, left).

> Figure 1. Clipping the signals to their main physical contents. 
<img src="https://user-images.githubusercontent.com/61290423/210902751-f54c43f5-3a4a-4eca-b091-af74226b83ae.png" alt="SAC signals plot" width="750"/>
     
As shown, earthquake-related signal (visible via the amplitude) lasting a few minutes may be included as part of 1-hour-long SAC file. The long file duration has a direct impact on file size as well as computer-resources in later calculations, so it is necessary to reduce the signal to its significant part within a time window bracketing the earthquake.        
Signal clipping may however fail if metadata on signal start date are inconsistent with earthquake occurrence date. I therefore identified and corrected such inconsistencies where they existed in few SAC files, which allowed for signal clipping and file size reducing (Figure 1, right), and used PZ files and the "earthquake info" table to rectify other spotted mismatches (station names, geocoordinates ... ).

With the SAC files now ready for subsequent calculations, I re-organized the files using a Tree Structure which, firstly reflected the file management on the computer, secondly allowed to track via the Tree branches the different settings associated to the SAC files (monitoring channel, ground motion type ...).

I then created a code sequence to calculate the several engineering parameters requested while keeping track of the SAC file characteristics they are associated to (station name, earthquake name etc ...). The parameters related to spectral accelerations and Arias intensity were obtained using functions defined at the early beginning of the script.

### Results

&nbsp; &nbsp; &nbsp; &nbsp; I finally obtained a table containing the different engineering parameters. Based on this table, I defined different plotting options to provide insights on the earthquake ground motions, notably one showing selected engineering parameters (Figure 2), and another one describing the impact of different filter corner frequencies as well as seismic attenuation across distance (Figure 3).





> Figure 2. Selecting and plotting PGA PGV and PGD
<img src="https://user-images.githubusercontent.com/61290423/211605116-49ff43f6-3e3c-4c2e-b13c-48f7057c4438.png" width="750"/>



> Figure 3. Plotting selected parameters to show impact of filtering and epicentral distance.
<img src="https://user-images.githubusercontent.com/61290423/210989406-c875d05f-8d55-4798-97ef-c0e399dbab50.png" height="700"/>
     

      





I manipulated the SAC files thanks to the "eseis" package, accessing not only the earthquake ground motion datapoints but also all metadata associated to the earthquake record such as seismic network, station instrumentation details, signal sampling time step ... .      
Before calculating the engineering parameters, few operations are required on the SAC signals. Notably, seismograms including the earthquake signal may cover a duration (~1 hour) much longer than the earthquake signal (few minutes) owing to the seismic station settings. To avoid calculations whose time and computer resources consumption grow with the signal duration, signal clipping (ie reducing) is one of these required operations. Other operations include assigning SAC signals to stations metadata, detecting and correcting inconsistencies in SAC signals metadata ... .    
To take into account the different types of SAC files (earthquake event, monitoring stations, instrumentation channel, motion types etc...) populating the different folders, a Tree Structure was set to classify the SAC files and keep track of all types and characteristics when later calculating the engineering parameters. In other words, which earthquake amplitudes were registered for which station/earthquake/filtering frequencies etc ... . 
When all SAC signals are "ready" and suitably sorted in the Tree, they can be sent through a routine to extract the desired engineering parameters.



Steps in the code are as follows :
+ defining functions used for later engineering parameter calculations (Arias intensity, spectral accelerations),
+ building the file architecture in order to read SAC files (containing the earthquake signals) and the PZ files containing recording stations metadata,
+ collecting earthquake event details and corresponding recording stations features (coordinates, distance ...),
+ retrieving and assigning stations metadata to SAC files,
+ clipping the signals to their main physical contents in order to reduce file size and to ease data manipulation and calculation,
+ creating a tree architecture to keep track of identification characteristics such as earthquake name, filtering frequency, station recording channel etc ... ,
+ extracting earthquake engineering parameters and derive other parameters from signal processing,
+ plotting selected parameters.

