This code reads and processes earthquake acceleration signals produced by the Neo-Deterministic Seismic Hazard Analysis (NDSHA) methodology (Panza 2001,2012).

The steps are as follows : 
+ creating functions for later operations (Arias Intensity, spectral accelerations, dataframe and list re-organizing)
+ preparing the file architecture for subsequent reading, 
+ reading and processing dataframe, 
+ creating a Tree Structure for organizing the files, 
+ processing files to extract relevant earthquake parameters such as PGA, spectral acceleration values, Arias-intensity related parameters.


A Tree Structure was decided from the onset, owing to the fact that acceleration data for each earthquake simulation (*.plt files) do not contain any metadata such as earthquake realisation number, event, directivity or station ID. The signals were therefore created using a tree structure holding such info using the directory and sub-directory names (ex: bil000/quake_dudley_Bchan/real01/). Signals were then read into R and processed following the same tree pattern.
  In terms of tree structure and how it connects to simulation processes, for each earthquake event and corresponding recording stations, 99 realisations of earthquake fault rupture are simulated for each of the 3 directivity angles 0, 90 and 180 degrees. In this code, as an example, only 3 realisations for each directivity angle and earthquake event are shown and processed.

