***In a nutshell :*** *I retrieved time histories from UK seismic monitoring network recorded in SAC (Seismic Analysis Coode) format and created a sequence in SAC code to process these RAW time histories. The purpose is to obtain processed, station metadata-convolved and filtered SAC files providing displacement, velocity and acceleration ground motions, ready to be futher processed in an external software.*


### Background

&nbsp; &nbsp; &nbsp; &nbsp; Earthquake signals (i.e. time histories, seismograms) are available as SAC files (Seismic Analysis Code) on seismology database platforms such as IRIS Wilber 3 or the British Geological Survey earthquake database. These signals "contain" the physics of earthquake; they can provide relevant engineering earthquake parameters commonly used as inputs in anti-seismic design or be directly used as dynamic load in modern structural analysis. Prior to serving their engineering purpose though, it is paramount for these signals to be adequately processed in order to be meaningful. 


### Task

&nbsp; &nbsp; &nbsp; &nbsp; The signals are provided for each station on separate horiozontal components (East-West EW and North-South NS), on two different monitoring channels (broadband B-channel and high-broadband H-channel). The signals are also recorded by seismic monitoring stations (including amplification gain and sensitivity), and only show the device response to the earthquake ground motion excitation (Figure 1).       
These signals do not show meaningful information, as they come with apparently inconsistent amplitude units and it is unclear which type of signal they refer to (velocity, acceleration ?), therefore they cannot be used directly. Instead, they should be combined with seismic station metadata (so-called PoleZero PZ files) in order to obtain relevant units and identify the ground motion nature of the signals. Besides, signal noise removal and signal filtering should also be applied.


> Figure 1. Plotting the raw SAC signal for the 2018 Swansea earthquake recorded atstation MCH1.

<img src="https://user-images.githubusercontent.com/61290423/214091098-6c2dfc43-76f2-4a15-bebe-90fa983a8db1.PNG" alt="raw signal" width="600"/>


### Action

&nbsp; &nbsp; &nbsp; &nbsp; I created a code routine "sample.sac" in the SAC directory which goes through the steps of SAC processing, including reading, deconvolving earthquake signal with station metadata (Figure 2), filtering (Figure 3), deriving other signals (e.g. acceleration) and plotting of seismic signals. The SAC signal for station MCH1 on H-channel is processed as an example. Connecting the station signal to corresponding metadata indicates the RAW signals are velocity time histories (Figure 4).       
This SAC routine is then incorporated in several macro sequences applicable to all stations and monitoring channels.  


> Figure 2. Plotting the SAC signal deconvolved from the station metadata (visual differences from previous plot highlighted in red).

 <img src="https://user-images.githubusercontent.com/61290423/214099639-337418d3-0b50-4ada-98d2-30a7f2b97a87.PNG" alt="mch1.hhe deconvolved " width="600"/>


> Figure 3. Plotting the SAC signal following filtering using a Butterworth bandpass (visual differences from previous plot highlighted in red).

<img src="https://user-images.githubusercontent.com/61290423/214099684-0f6d3fdb-a89a-4533-b95d-574f035b6f70.PNG" alt="mch1 hhe unfiltered" width="600"/>





> Figure 4. Plotting velocity (post plot cosmetic improvements).

<img src="https://user-images.githubusercontent.com/61290423/214099728-414a91e8-f56d-4163-92a2-b412d45bd5a9.PNG" alt="filtered velocity signal" width="600"/>



&nbsp; &nbsp; In the directories /Bchan_25Hz and /Hchan_35Hz, I created other SAC sequences to create the required directories (acceleration/velocity/displacement), clean unnecessary files produced by the routine, as well as managing the signal plots. Since RAW signal are velocity files, the set of commands in "sample.sac" was then adapted to derive or integrate the earthquake signals to obtain acceleration (signal_vel_acc.mac, Figure 5) or displacement (signal_dis.mac) files as required.

All sequences are finally incorporated into a larger macro sequence "loopXXX.mac" (e.g. loopH35.mac), which handles signals for all stations, monitoring channels and horizontal components.    


> Figure 5. Plotting acceleration (after plot cosmetic improvements on axes and captions).

<img src="https://user-images.githubusercontent.com/61290423/214099750-83829289-465d-434b-bb10-ad9dcf2b72d7.PNG" alt="filtered acc signal" width="600"/>


### Results

&nbsp; &nbsp; &nbsp; &nbsp; The macro sequence supplies processed SAC files, representative of velocity, acceleration and displacement time histories, as well as PDF as corresponding plots. The full routine is using an example set of SAC files associated to the 2018 Mw 4.3 Swansea earthquake.     
SAC files produced by the macro routines serve as input in the routine shown in portfolio/seismic_records. Directories are structured based on the station recording channels (B-channel or H-channel) and upper frequency filtering (25Hz or 35Hz).








