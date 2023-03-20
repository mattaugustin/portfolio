
### Background

&nbsp; &nbsp; &nbsp; &nbsp; This work relates to a section of my PhD project and aims at developing a regional seismic hazard map for the United Kingdom (UK). The work started with calculating Peak Ground Acceleration (PGA) values at nodes of a regularly-spaced grid (0.2 °) mapped over the UK land, using the Neo-Deterministic Seismic Hazard Analysis methodology detailed in Panza *et al.* (2001, 2012). The gridded PGA distribution, plotted in Figure 1, provides useful insight on seismic hazard levels, but lacks the precision to provide the seismic hazard level at precise locations (e.g Edinburgh) not covered by the distribution. Such precision for any geocoordinates associated to the UK land can be achieved by smoothing the distribution available.  


> Figure 1. Map of Peak Ground Acceleration values gridded over the United Kingdom.

<img src="https://user-images.githubusercontent.com/61290423/208784049-dc958df0-a162-48df-82d0-ffe6f84c42b8.png" alt="krigged map" width="600"/>


### Task

&nbsp; &nbsp; &nbsp; &nbsp; Using the PGA distribution available, the goal is to smooth the hazard data over the entire UK land using a krigging technique, ensuring that potential high hazard individual nodes are not smoothed out in the final distribution prior to the plotting of hazard map.


### Action

&nbsp; &nbsp; &nbsp; &nbsp; Steps include preparing relevant spatial items (dataframes, grid, raster), creating and fitting a variogram for the kriging operation, assessing PGA value range, creating amplitude bins and plotting the krigged hazard map.     


### Results

&nbsp; &nbsp; &nbsp; &nbsp; The krigged map is shown in Figure 2. The individual high hazard points located in Northwest Wales and near Dover are visible in the forms of small hazard pockets in the krigged map, proving that the smoothing technique did not smooth away these individual points.

> Figure 2. Map of Peak Ground Acceleration values krigged over the United Kingdom.

<img src="https://user-images.githubusercontent.com/61290423/208784062-dc5e7056-aa57-406f-ab36-2365e571a204.png" alt="krigged map" width="600"/>


     
### References           
Panza, G.F., Romanelli, F. and Vaccari, F., 2001. Seismic wave propagation in laterally heterogeneous anelastic media: theory and applications to seismic zonation. In Advances in geophysics (Vol. 43, pp. 1-95). Elsevier.

Panza, G.F., La Mura, C., Peresan, A., Romanelli, F. and Vaccari, F., 2012. Seismic hazard scenarios as preventive tools for a disaster resilient society. In Advances in geophysics (Vol. 53, pp. 93-165). Elsevier.





