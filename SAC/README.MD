

  Earthquake signals (i.e. time histories) are available as SAC files (Seismic Analysis Code) on seismology database platforms such as IRIS Wilber 3.
These signals are the first step towards extracting and/or deriving relevant engineering earthquake parameters commonly used as inputs in modern structural analysis. Suitable treatment of these earthquake signals is therefore paramount.
The signals are provided as recorded by seismic monitoring stations (including amplification gain and sensitivity) and thus should be examined in association with seismic station metadata in order to obtain relevant units. Upon adequate processing of the SAC files, the user can then obtain relevant time-series data (acceleration, velocity...) that can be post-processed using different code solution (notably using R code on this Github, in portfolio/seismic_records).
  
  The code routine in this SAC directory shows the detailed steps of SAC processing, including reading, deconvolving earthquake signal with station metadata, filtering, deriving other signals (e.g. acceleration) and plotting of seismic signals.
  
  This code mainly provides the SAC files serving as input in the routine shown in portfolio/seismic_records. Directories are structured based on the station recording channels (B or H, for broadband and high broadband) and upper frequency filtering (25Hz or 35Hz), which determined the same structure observed in portfolio/seismic_records. An example set of SAC files associated to the 2018 Mw 4.3 Swansea earthquake is supplied for reference.
 

Plotting the raw SAC signal for station.

![mch1_hhe_raw](https://user-images.githubusercontent.com/61290423/214091098-6c2dfc43-76f2-4a15-bebe-90fa983a8db1.PNG)


Plotting the SAC signal deconvolved from the station metadata.

![mch1_hhe_deconvolved](https://user-images.githubusercontent.com/61290423/214099639-337418d3-0b50-4ada-98d2-30a7f2b97a87.PNG)
 

Plotting the SAC signal following filtering using a Butterworth bandpass.

![mch1_hhe_decon_filtered](https://user-images.githubusercontent.com/61290423/214099684-0f6d3fdb-a89a-4533-b95d-574f035b6f70.PNG)


Plotting velocity (post plot cosmetic improvements).

![mch1_hhe_decon_filt_vel](https://user-images.githubusercontent.com/61290423/214099728-414a91e8-f56d-4163-92a2-b412d45bd5a9.PNG)


Plotting acceleration (post plot cosmetic improvements).

![mch1_hhe_decon_filt_acc](https://user-images.githubusercontent.com/61290423/214099750-83829289-465d-434b-bb10-ad9dcf2b72d7.PNG)


