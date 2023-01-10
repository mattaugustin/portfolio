This code reads and processes earthquake acceleration signals produced by the Neo-Deterministic Seismic Hazard Analysis (NDSHA) methodology (Panza 2001,2012).

The steps are as follows : 
+ creating functions for later operations (Arias Intensity, spectral accelerations, dataframe and list re-organizing)
+ preparing the file architecture for subsequent reading, 
+ reading and processing dataframe, 
+ creating a Tree Structure for organizing the files, 
+ processing files to extract relevant earthquake parameters such as PGA, spectral acceleration values, Arias-intensity related parameters,
+ plotting parameters related to earthquake events.

An example dataset with earthquake signals corresponding to different events, recording channels and stations is provided at the following link:
https://drive.google.com/drive/folders/1jdH0llIJzKXL89UEYmAZGeMVjHX_uXgJ?usp=sharing

A Tree Structure was decided from the onset, owing to the fact that acceleration data for each earthquake simulation (*.plt files) do not contain any metadata such as earthquake realisation number, event, directivity or station ID. The signals were therefore created using a tree structure holding such info using the directory and sub-directory names (ex: bil000/quake_dudley_Bchan/real01/). Signals were then read into R and processed following the same tree pattern.
  In terms of tree structure and how it connects to simulation processes, for each earthquake event and corresponding recording stations, 99 realisations of earthquake fault rupture are simulated for each of the 3 directivity angles 0, 90 and 180 degrees. In this code, as an example, only 3 realisations (instead of 99) for each of directivity angle are considered, and earthquake signals correspond to 4 separate seismic events, each of which recorded on 3 different stations.

References: 
+ Panza, G.F., Romanelli, F. and Vaccari, F., 2001. Seismic wave propagation in laterally heterogeneous anelastic media: theory and applications to seismic zonation. In Advances in geophysics (Vol. 43, pp. 1-95). Elsevier.

+ Panza, G.F., La Mura, C., Peresan, A., Romanelli, F. and Vaccari, F., 2012. Seismic hazard scenarios as preventive tools for a disaster resilient society. In Advances in geophysics (Vol. 53, pp. 93-165). Elsevier.

> Geology distribution associated to the seismic monitoring stations

<img src="https://user-images.githubusercontent.com/61290423/211395831-d7d4aeed-bae8-4f52-8a7b-dfe53c04a143.png" alt="geo distr" width="750"/>

> PGA Amplitude distribution associated to the seismic wavefront directivity angle in case of the Swansea earthquake

<img src="https://user-images.githubusercontent.com/61290423/211396042-991788bc-cfcd-451a-ab35-22b58d07fc4e.png" alt="PGA amp distr" width="750"/>

> Map of earthquakes and associated monitoring stations

<img src="https://user-images.githubusercontent.com/61290423/211563354-b06ce175-33f4-4b98-9f8a-5b616f313949.png" alt="map quakes station" width="750"/>

