***In a nutshell :*** *I created this code as a final step in analysing the performance of the Modal Summation technique featured in the NDSHA methodology. 
Using previous section results, I firstly aggregated the dataset obtained from the MS simulated earthquake ground motions, then joined them to the observation dataset, allowing the calculation and analysis of residuals between earthquake ground motion observations and simulations. I included several visualization options showing not only the comparison prediction-observation but also impact of several simulation or/and observation settings (rupture front directivity angles, corner frequencies ...).*


### Background
&nbsp; &nbsp; &nbsp; &nbsp; In relevant project companion sections, I firstly assembled an earthquake observation dataset made of several engineering parameters derived from the processing of SAC files associated to recent UK earthquakes. I then assembled an earthquake simulation dataset taking into account simulation specifics (rupture front directivity angles and randomness in fault-rupturing) which resulted in a large dataset of engineering parameters in correspondance with the observation dataset. Both datasets were stepping stones in comparing earthquake ground motion observations vs predictions.


### Dataset
The earthquake observation dataset is included in this directory.
The starting dataset of simulated ground motions prior to aggregation is located at the following link owing to its size:
https://drive.google.com/drive/folders/14DTBeyYuI1QqDD22NuBmGvperO-_bVv_?usp=sharing


### Task
&nbsp; &nbsp; &nbsp; &nbsp; In this work, the large simulation dataset of engineering parameters derived from MS time histories ("simulation dataset") is to be aggregated, prior to comparing with the earthquake observation dataset. The comparison is to be done by means of residual analysis.


### Action
&nbsp; &nbsp; &nbsp; &nbsp;The work starts with managing the simulation dataset; it contains, for each combination of the four characteristics earthquake/monitoring channel/station/horizontal component (E/MC/S/HC), several engineering parameters whose amplitudes vary based on the directivity angles and the fault-rupturing random realisation index (the meaning of these two earthquake properties is discussed in portfolio/NDSHA). In order to compare engineering parameters from the simulation dataset to the observation dataset, a reduced number of values for each combination E/MC/S/HC is desired.         

I therefore manipulated the dataset, first grouping amplitudes by dropping the column "real" and "directivity", and then grouping by dropping the column "real", determining the aggregated mean and median value of the entire amplitude range. This means that the former grouping will show aggregated values from all 300 amplitudes, and the latter grouping will show aggregated values from 100 amplitudes while still showing the impact of the parameter "directivity".
This provides a dataset in which the impact of "directivity" can still be observed in the parameter amplitudes, with aggregated values associated to a specific directivity angle (e.g. directivity="bil_dir000" for directivity angle 0) or reflecting all directivities (directivity="combined"). This allows me to assess the impact of directivity in the residual calculations later in the code.    
                    
I then joined the aggregated simulation dataset to the observation dataset and, upon selecting some variables of interest (engineering parameter, filtering corner frequency for the observations...), plotted observations and simulations together for all earthquakes (Figure 1). 

> Figure 1. PGA observations against median of all 300 predicted values, using observation signals filtered at 35 Hz.

<img src="https://user-images.githubusercontent.com/61290423/216592445-2f4fe860-b52b-4903-af52-d7749dc24a5f.png" alt="ndsha aggregated vs sac" width="600"/>


I then calculated residuals defined as log<sub>10</sub>( predictions / observations ) , using both aggregated mean and median. I additionally calculated residual averages across distance bins in order to obtain additional insight on the performance of the predictive model across distance, adding the soil geology as a differentiating setting in the calculation.

I finally coded several visualization options to assess the impact of various settings contained in the table of residuals between prediction and simulations.


### Results
&nbsp; &nbsp; &nbsp; &nbsp; Ultimately, I obtained and saved a table containing residuals between predictions and observations for a range of engineering parameters. The table columns and organisation help to keep track of various settings possibly affecting the residuals such as directivity, distance, soil geology ... .         
The final residual table and the visualization code sequences allow me to see the impact of several settings including : directivity in simulated motions, filtering corner frequency in observations, soil geology, distance ... (Figure 2). I also added a sequence to plot observations against the initial range of all 300 predicted amplitude values per engineering parameter (Figure 3).


> Figure 2. Impact of geologies on residual distribution

<img src="https://user-images.githubusercontent.com/61290423/216591242-c44b89c4-5569-4547-a456-fc8b9909fa04.png" alt="all geo vs rock" width="600"/>



> Figure 3. PGA Observations compared against full range of predicted values

<img src="https://user-images.githubusercontent.com/61290423/216591159-c656777d-5d7a-4933-ad10-209a602c3246.png" alt="ndsha swarm" width="600"/>



