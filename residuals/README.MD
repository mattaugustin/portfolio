This work is the last step in assessing how well the NDSHA synthetic earthquake signals perform when compared to actual seismic observations recorded in the UK, by means of a residual analysis.

In this work, the dataset of predictions from NDSHA signals (300 values for each engineering parameter, associated to the 99 realisations of fault rupture on each of the three directivity angles 0, 90 and 180) is aggregated to obtain the median and mean value of each parameter (e.g. Peak Ground Acceleration).
The prediction dataset is then compared against the observation dataset (based on SAC signals) in order to calculate corresponding residuals.

This work constitutes the final part of analysing earthquake records, coming after: 
+ I) Processing SAC signals (*portfolio/SAC/*),
+ II) Deriving earthquake parameters from SAC signals (*portfolio/seismic_records/*), and 
+ III) Deriving parameters from NDSHA-based earthquake signals (*portfolio/NDSHA/*).

Plotting the residuals provides the user a glimpse about how performant the predictive model (NDSHA) behaves, but also shows the impact of several settings such as directivity, geologies ...

The starting dataset prior to aggregation is located at the following link owing to its size:
https://drive.google.com/drive/folders/14DTBeyYuI1QqDD22NuBmGvperO-_bVv_?usp=sharing


> PGA observations against median of all 3000 predicted values, using signals filtered at 35 Hz.

![NDSHA_vs_SAC_median_small](https://user-images.githubusercontent.com/61290423/216592445-2f4fe860-b52b-4903-af52-d7749dc24a5f.png)


> PGA Observations compared against full range of predicted values

![NDSHA_vs_SAC_swarm_small](https://user-images.githubusercontent.com/61290423/216591159-c656777d-5d7a-4933-ad10-209a602c3246.png)


> Impact of geologies on residual distribution

![Allgeo_vs_Rock_small](https://user-images.githubusercontent.com/61290423/216591242-c44b89c4-5569-4547-a456-fc8b9909fa04.png)



