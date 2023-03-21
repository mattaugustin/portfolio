
### Background

&nbsp; &nbsp; &nbsp; &nbsp; This work relates to a section of my PhD project and aims at developing a regional seismic hazard map for the United Kingdom (UK). The work started with calculating Peak Ground Acceleration (PGA) values at nodes of a regularly-spaced grid (0.2 °) mapped over the UK land, using the Neo-Deterministic Seismic Hazard Analysis methodology detailed in Panza *et al.* (2001, 2012). The gridded PGA distribution, plotted in Figure 1, provides useful insight on seismic hazard levels, but lacks the precision to provide the seismic hazard level at precise locations (e.g Edinburgh) not covered by the distribution. Such precision for any geocoordinates associated to the UK land can be achieved by smoothing the distribution available.  


> Figure 1. Map of Peak Ground Acceleration values gridded over the United Kingdom.

<img src="https://user-images.githubusercontent.com/61290423/208784049-dc958df0-a162-48df-82d0-ffe6f84c42b8.png" alt="krigged map" width="600"/>


### Task

&nbsp; &nbsp; &nbsp; &nbsp; Using the PGA distribution available, the goal is to smooth the hazard data over the entire UK land using a kriging technique, ensuring that potential high hazard individual nodes are not smoothed out in the final distribution prior to the plotting of hazard map. A corresponding plot will help visualize the results.


### Action

&nbsp; &nbsp; &nbsp; &nbsp; Any smoothing technique (krigging, nearest-neighbour ...) requires coordinates to be defined using a Cartesian definition, a.k.a. a distance measure, as well as a "surface" to smooth the distribution over. Defining these systems is therefore the first step.
The seismic PGA distribution is associated to GPC coordinates using latitude and longitude, hence I converted it to Cartesian coordinates using the British National Grid (BNG) system. In parallel, the smoothing surface is the UK land, loaded as a shape file, which I then converted to Cartesian coordinates based on BNG and rasterized to a user-defined resolution.
As implied by kriging interpolation technique, I created an empirical variogram to visualize the correlation PGA hazard-distance, and then defined a variogram model (between exponential, spherical or matern) which could fit the empirical variogram (Figure 2).     

> Figure 2. Empirical variogram and variogram model.

<img src="https://user-images.githubusercontent.com/126000617/226683834-6ec46f3f-2353-4216-8fd9-21fbc892cc7f.png" alt="variogram" width="500"/>


With the variogram model in hand, I could then spatially interpolate the existing PGA distribution over the UK land.
Interpolation values "cover" slightly beyond the borders of UK land so I intersected this now-interpolated (i.e. krigged) distribution with the UK raster map, finally obtaining the PGA prediction distribution shown in Figure 3. The default colour code on the predictive map obviously contrasts with the original gridded PGA distribution (Figure 1) so needs to be adapted. 

> Figure 3. Plot of krigged PGA values before and after intersection with UK map. On the right plots, maps show the range of predicted values (reflecting PGA) and the low variance of these predictions across the land.

<img src="https://user-images.githubusercontent.com/126000617/226683196-f961562f-a5b1-48e7-bed3-3a2aff8e59e5.png" alt="kriged distribution" width="1000"/>


A histogram plotting the range of available krigged PGA values helps to decide on possible new PGA amplitude bins for subsequent plotting (Figure 4). Eventually, the same bins used in Figure 1 were selected and applied to the predicted PGA dataset.

> Figure 4. Histogram showing bins of PGA values in the krigged distribution.

<img src="https://user-images.githubusercontent.com/126000617/226682477-dcac9663-f271-4113-a08b-77ba4623466c.png" alt="histogram" width="500"/>



Plotting the krigged map required importing additional items such as map of the UK borders and surrounding countries in the form of dataframes (for use with ggplot).
I added a customized plotting style to improve font and font size, colour palette and titles.


### Results

&nbsp; &nbsp; &nbsp; &nbsp; The krigged map is shown in Figure 5. The individual high hazard points located in Northwest Wales and near Dover are visible in the forms of small hazard pockets in the krigged map, proving that the smoothing technique did not smooth away these individual points.

> Figure 5. Map of Peak Ground Acceleration values krigged over the United Kingdom.

<img src="https://user-images.githubusercontent.com/61290423/208784062-dc5e7056-aa57-406f-ab36-2365e571a204.png" alt="krigged map" width="600"/>


     
### References           
Panza, G.F., Romanelli, F. and Vaccari, F., 2001. Seismic wave propagation in laterally heterogeneous anelastic media: theory and applications to seismic zonation. In Advances in geophysics (Vol. 43, pp. 1-95). Elsevier.

Panza, G.F., La Mura, C., Peresan, A., Romanelli, F. and Vaccari, F., 2012. Seismic hazard scenarios as preventive tools for a disaster resilient society. In Advances in geophysics (Vol. 53, pp. 93-165). Elsevier.





