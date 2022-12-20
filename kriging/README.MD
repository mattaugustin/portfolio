This work relates to a section of my PhD project and aims at developing a regional seismic hazard map for the United Kingdom (UK). The work started with calculating Peak Ground Acceleration (PGA) values at nodes of a regularly-spaced grid (0.2 °) mapped over the UK land, using the Neo-Deterministic Seismic Hazard Analysis methodology detailed in Panza et al (2001, 2012). 
In the proposed code, the gridded PGA distribution was then smoothed over the entire UK land using a kriging technique, ensuring that potential high hazard individual nodes are not smoothed out in the final distribution prior to the plotting of hazard map.
Steps include preparing relevant spatial items (dataframes, grid, raster), creating and fitting a variogram for the kriging operation, assessing PGA value range, creating amplitude bins and plotting the kriged hazard map.

References:
Panza, G.F., Romanelli, F. and Vaccari, F., 2001. Seismic wave propagation in laterally heterogeneous anelastic media: theory and applications to seismic zonation. In Advances in geophysics (Vol. 43, pp. 1-95). Elsevier.

Panza, G.F., La Mura, C., Peresan, A., Romanelli, F. and Vaccari, F., 2012. Seismic hazard scenarios as preventive tools for a disaster resilient society. In Advances in geophysics (Vol. 53, pp. 93-165). Elsevier.

original PGA value distribution
![regiomap_pga_original_small](https://user-images.githubusercontent.com/61290423/208784049-dc958df0-a162-48df-82d0-ffe6f84c42b8.png)


kriged PGA distribution
![regiomap_pga_kriged - small](https://user-images.githubusercontent.com/61290423/208784062-dc5e7056-aa57-406f-ab36-2365e571a204.png)
