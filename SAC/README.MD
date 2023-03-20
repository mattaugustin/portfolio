***In a nutshell :*** *I retrieved time histories from UK seismic monitoring network recorded in SAC (Seismic Analysis Coode) format and created a sequence in SAC code to process these RAW time histories. The purpose is to obtain processed, station metadata-convolved and filtered SAC files providing displacement, velocity and acceleration ground motions, ready to be futher processed in an external software.*


### Background

&nbsp; &nbsp; &nbsp; &nbsp; Earthquake signals (i.e. time histories, seismograms) are available as SAC files (Seismic Analysis Code) on seismology database platforms such as IRIS Wilber 3 or the British Geological Survey earthquake database. These signals "contain" the physics of earthquake; they can provide relevant engineering earthquake parameters commonly used as inputs in anti-seismic design or be directly used as dynamic load in modern structural analysis. Prior to serving their engineering purpose though, it is paramount for these signals to be adequately processed in order to be meaningful. 


### Task

&nbsp; &nbsp; &nbsp; &nbsp; The signals are provided for each station on separate horiozontal components (East-West and North-South), as recorded by seismic monitoring stations (including amplification gain and sensitivity) , and only show the device response to the earthquake ground motion excitation (Figure 1). These signals do not show meaningful information - with apparently inconsistent amplitude units and confusing about which type of signal they refer to (velocity, acceleration ?) -, and cannot therefore be used directly. Instead, they should be combined with seismic station metadata (so-called PoleZero PZ files) in order to obtain relevant units and identify the ground motion nature of the signals. Besides, signal noise removal and signal filtering should also be applied.


> Figure 1. Plotting the raw SAC signal for the 2018 Swansea earthquake recorded atstation MCH1.

<img src="https://user-images.githubusercontent.com/61290423/214091098-6c2dfc43-76f2-4a15-bebe-90fa983a8db1.PNG" alt="raw signal" width="600"/>


Upon adequate processing of the SAC files, the user can then obtain relevant time-series data (acceleration, velocity...) that can be post-processed using different code solution (notably using R code on this Github, in portfolio/seismic_records).


### Action

&nbsp; &nbsp; &nbsp; &nbsp; The code routine "sample.sac" in the SAC directory provides a glimpse in the steps of SAC processing, including reading, deconvolving earthquake signal with station metadata, filtering, deriving other signals (e.g. acceleration) and plotting of seismic signals.
In the directories /Bchan_25Hz and /Hchan_35Hz, macro SAC sequences "loopX.mac" using smaller sequences and similar commands as in sample.sac are available. The macro routines supply acceleration, velocity and displacement signals, as well as plots in pdf format, for an example set of SAC files associated to the 2018 Mw 4.3 Swansea earthquake used as reference.


### Results

SAC files produced by the macro routines serve as input in the routine shown in portfolio/seismic_records. Directories are structured based on the station recording channels (B or H, for broadband and high broadband) and upper frequency filtering (25Hz or 35Hz), which determined the same structure observed in portfolio/seismic_records. 
 


> Figure 2. Plotting the SAC signal deconvolved from the station metadata (visual differences from previous plot highlighted in red).

![mch1_hhe_deconvolved](https://user-images.githubusercontent.com/61290423/214099639-337418d3-0b50-4ada-98d2-30a7f2b97a87.PNG)
 <img src="https://user-images.githubusercontent.com/61290423/214091098-6c2dfc43-76f2-4a15-bebe-90fa983a8db1.PNG" alt="raw signal" width="600"/>


> Figure 3. Plotting the SAC signal following filtering using a Butterworth bandpass (visual differences from previous plot highlighted in red).

![mch1_hhe_decon_filtered](https://user-images.githubusercontent.com/61290423/214099684-0f6d3fdb-a89a-4533-b95d-574f035b6f70.PNG)
<img src="https://user-images.githubusercontent.com/61290423/214091098-6c2dfc43-76f2-4a15-bebe-90fa983a8db1.PNG" alt="raw signal" width="600"/>


> Figure 4. Plotting velocity (post plot cosmetic improvements).

![mch1_hhe_decon_filt_vel](https://user-images.githubusercontent.com/61290423/214099728-414a91e8-f56d-4163-92a2-b412d45bd5a9.PNG)
<img src="https://user-images.githubusercontent.com/61290423/214091098-6c2dfc43-76f2-4a15-bebe-90fa983a8db1.PNG" alt="raw signal" width="600"/>


> Figure 5. Plotting acceleration (post plot cosmetic improvements).

![mch1_hhe_decon_filt_acc](https://user-images.githubusercontent.com/61290423/214099750-83829289-465d-434b-bb10-ad9dcf2b72d7.PNG)
<img src="https://user-images.githubusercontent.com/61290423/214091098-6c2dfc43-76f2-4a15-bebe-90fa983a8db1.PNG" alt="raw signal" width="600"/>


