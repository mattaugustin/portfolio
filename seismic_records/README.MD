This section focuses on acquiring the metrics related to earthquake observation, i.e. the engineering parameters such as Peak Ground Acceleration (PGA) and spectral accelerations (SA) at selected structural periods.

Previously (see portfolio/SAC), SAC signals were processed and saved in SAC format, and a set of files was available, in connection with horizontal component type (East-West or North-South), instrumentation channel (broadband or high broadband), filtering frequency (upper corner frequency), earthquake motion type (acceleration...). Further processing of SAC files is achieved in R code in order to ultimately obtain earthquake parameters of relevant engineering interest (PGA etc...).      
In this section, the code shows how to read and manipulate SAC files thanks to the "eseis" package, accessing not only the earthquake ground motion datapoints but also all metadata associated to the earthquake record such as seismic network, station instrumentation details, signal sampling time step ... . Before calculating the engineering parameters, few operations are required on the SAC signals. Notably, seismograms including the earthquake signal may cover a duration (~1 hour) much longer than the earthquake signal (few minutes) owing to the seismic station settings. To avoid calculations whose time and computer resources consumption grow with the signal duration, signal clipping (ie reducing) is one of these required operations. Other operations include assigning SAC signals to stations metadata, detecting and correcting inconsistencies in SAC signals metadata ... .    
To take into account the different types of SAC files (earthquake event, monitoring stations, instrumentation channel, motion types etc...) populating the different folders, a Tree Structure was set to classify the SAC files and keep track of all types and characteristics when later calculating the engineering parameters. In other words, which earthquake amplitudes were registered for which station/earthquake/filtering frequencies etc ... . 
When all SAC signals are "ready" and suitably sorted in the Tree, they can be sent through a routine to extract the desired engineering parameters.


With this code is provided a set of SAC earthquake signals and associated station metadata referred to as PoleZero files (PZ).
 https://drive.google.com/drive/folders/1FRe2bp1G7aQ-0hjDYUV57vSaD_2b9zix?usp=sharing

Steps in the code are as follows :
+ defining functions used for later engineering parameter calculations (Arias intensity, spectral accelerations),
+ building the file architecture in order to read SAC files (containing the earthquake signals) and the PZ files containing recording stations metadata,
+ collecting earthquake event details and corresponding recording stations features (coordinates, distance ...),
+ retrieving and assigning stations metadata to SAC files,
+ clipping the signals to their main physical contents in order to reduce file size and to ease data manipulation and calculation,
+ creating a tree architecture to keep track of identification characteristics such as earthquake name, filtering frequency, station recording channel etc ... ,
+ extracting earthquake engineering parameters and derive other parameters from signal processing,
+ plotting selected parameters.
     
      
      
> clipping the signals to their main physical contents 
<img src="https://user-images.githubusercontent.com/61290423/210902751-f54c43f5-3a4a-4eca-b091-af74226b83ae.png" alt="SAC signals plot" width="750"/>
     
       
> plotting selected parameters to show impact of filtering and epicentral distance
<img src="https://user-images.githubusercontent.com/61290423/210989406-c875d05f-8d55-4798-97ef-c0e399dbab50.png" width="750"/>
     
      
> Selecting and plotting PGA PGV and PGD
<img src="https://user-images.githubusercontent.com/61290423/211605116-49ff43f6-3e3c-4c2e-b13c-48f7057c4438.png" width="750"/>
