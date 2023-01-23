

  Earthquake signals (i.e. time histories) are available as SAC files (Seismic Analysis Code) on seismology database platforms such as IRIS Wilber 3.
These signals are the first step towards extracting and/or deriving relevant engineering earthquake parameters commonly used as inputs in modern structural analysis. Suitable treatment of these earthquake signals is therefore paramount.
The signals are provided as recorded by seismic monitoring stations (including amplification gain and sensitivity) and thus should be examined in association with seismic station metadata in order to obtain relevant units. Upon adequate processing of the SAC files, the user can then obtain relevant time-series data (acceleration, velocity...) that can be post-processed using different code solution (notably using R code on this Github, in portfolio/seismic_records).
  
  The code routine in this SAC directory shows the detailed steps of SAC processing, including reading, deconvolving earthquake signal with station metadata, filtering, deriving other signals (e.g. acceleration) and plotting of seismic signals.
  
  This code mainly provides the SAC files serving as input in the routine shown in portfolio/seismic_records. Directories are structured based on the station recording channels (B or H, for broadband and high broadband) and upper frequency filtering (25Hz or 35Hz), which determined the same structure observed in portfolio/seismic_records. An example set of SAC files associated to the 2018 Mw 4.3 Swansea earthquake is supplied for reference.
 

Plotting the raw SAC signal for station.
![mch1_hhe_raw](https://user-images.githubusercontent.com/61290423/214091098-6c2dfc43-76f2-4a15-bebe-90fa983a8db1.PNG)


Plotting the SAC signal deconvolved from the station metadata.
![mch1_hhe_deconvolved](https://user-images.githubusercontent.com/61290423/214094772-81e5defa-2754-4a68-b6cc-31ffd6f62dbc.PNG)


Plotting the SAC signal following filtering using a Butterworth bandpass.
![mch1_hhe_decon_filtered](https://user-images.githubusercontent.com/61290423/214097171-1d456fb3-b64a-4f4b-bfd7-7f25cde1590d.PNG)

Plotting velocity (post plot cosmetic improvements).
![mch1_hhe_decon_filt_vel](https://user-images.githubusercontent.com/61290423/214097325-1ee08cbd-ac4e-4fba-84a8-3e5c8dec395c.PNG)


Plotting acceleration (post plot cosmetic improvements).
![mch1_hhe_decon_filt_acc](https://user-images.githubusercontent.com/61290423/214097349-31aa4651-a6a7-4f2f-afe4-8ff091ad6042.PNG)

