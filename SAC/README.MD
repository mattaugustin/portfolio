This code reads earthquake signals saved in SAC format (Seismic Analysis Code) from local directory folders and following several manipulations extracts and calculate common engineering parameters such as Peak Ground Acceleration (PGA) and spectral accelerations (SA) at selected structural periods.

Earthquake signals as recorded by seismic stations are usually available for both orthogonal horizontal components (e.g. East-West and North-South) in SEED or SAC format on dedicated websites. 
  Upon retrieving station signals associated to an earthquake, a SAC code routine (not provided here) was used to deconvolve the RAW signals from instrumental electronic amplifying using station metadata, filter the signals using bandpass filter, and finally derive accelerations, velocity and displacement signals. The resulting SAC earthquake signals can then be further utilized in software code such as R to obtain earthquake parameters of relevant engineering interest.
  In this example, signals are registered on broadband (B channel) or higher broadband (H channel) stations, depending on the regional monitoring network. The bandpass filtering frequency is set according to station instrumental recording rate and users decision, hence in this example at 35 Hz (or 25 Hz owing to station instrumental limitations) and 10 Hz.

With this code is provided a set of SAC earthquake signals and associated station metadata referred to as PoleZero files (PZ).
 https://drive.google.com/drive/folders/1FRe2bp1G7aQ-0hjDYUV57vSaD_2b9zix?usp=sharing

Steps are as follows :
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
