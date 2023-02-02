This work may constitute the final part of analysing earthquake records, coming after I) Processing SAC signals, II) Deriving earthquake parameters from SAC signals, and III) Deriving parameters from NDSHA-based earthquake signals.

In this work, the dataset of predictions from NDSHA signals (300 values for each engineering parameter, associated to the 99 realisations of fault rupture on each of the three directivity angles 0, 90 and 180) is aggregated to obtain the median and mean value of each parameter (e.g. Peak Ground Acceleration).
The prediction dataset is then compared against the observation dataset (based on SAC signals) in order to calculate corresponding residuals.
Plotting said residuals provides the user a glimpse about how performant the predictive model (NDSHA) behaves, but also shows the impact of several settings such as directivity, geologies ...


> PGA observations against median of all 3000 predicted values, using signals filtered at 35 Hz.
![NDSHA_vs_SAC_median](https://user-images.githubusercontent.com/61290423/216394375-de324b97-fbf2-4c64-b71c-82f957c3bae4.png)



> PGA Observations compared against full range of predicted values
![NDSHA_vs_SAC_swarm](https://user-images.githubusercontent.com/61290423/216393809-783966d6-4c4c-45ce-8e4f-25359daeae7a.png)



> Impact of geologies on residual distribution
![Allgeo_vs_Rock](https://user-images.githubusercontent.com/61290423/216394729-7627fa26-29b0-4237-a357-bcc7aa192590.png)
