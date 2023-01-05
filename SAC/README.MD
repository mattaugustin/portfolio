This code reads earthquake signals saved in SAC format (Seismic Analysis Code) from local directory folders and following several manipulations extracts and calculate common engineering parameters such as Peak Ground Acceleration (PGA) and spectral accelerations (SA) at selected structural periods.

Steps are as follows :
+ building the file architecture in order to read SAC files (containing the earthquake signals) and the PoleZero (PZ) files containing recording stations metadata
+ clipping the signals to reduce file size and to ease data manipulation and calculation
+ creating a tree architecture to keep track of identification characteristics such as earthquake name, filtering frequency, station recording channel etc ...
+ extracting earthquake engineering parameters and derive others from signal processing.
