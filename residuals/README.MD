This work may constitute the final part of analysing earthquake records, coming after 
+ I) Processing SAC signals (portfolio/SAC/),
+ II) Deriving earthquake parameters from SAC signals (portfolio/seismic_records/), and 
+ III) Deriving parameters from NDSHA-based earthquake signals (portfolio/NDSHA/).

In this work, the dataset of predictions from NDSHA signals (300 values for each engineering parameter, associated to the 99 realisations of fault rupture on each of the three directivity angles 0, 90 and 180) is aggregated to obtain the median and mean value of each parameter (e.g. Peak Ground Acceleration).
The prediction dataset is then compared against the observation dataset (based on SAC signals) in order to calculate corresponding residuals.

Plotting said residuals provides the user a glimpse about how performant the predictive model (NDSHA) behaves, but also shows the impact of several settings such as directivity, geologies ...

The starting dataset prior to aggregation is located at the following link owing to its size:
https://drive.google.com/drive/folders/14DTBeyYuI1QqDD22NuBmGvperO-_bVv_?usp=sharing


> PGA observations against median of all 3000 predicted values, using signals filtered at 35 Hz.

![Allgeo_vs_Rock_small](https://user-images.githubusercontent.com/61290423/216395826-9f872cb6-6ad7-4892-9436-d39ce33c0843.png)



> PGA Observations compared against full range of predicted values

![NDSHA_vs_SAC_swarm_small](https://user-images.githubusercontent.com/61290423/216395890-a760919b-f907-4ec5-b0cd-57e6c432930a.png)



> Impact of geologies on residual distribution

![Allgeo_vs_Rock_small](https://user-images.githubusercontent.com/61290423/216396185-5f61cc88-0257-4a1a-a567-4463bf96e9c3.png)

