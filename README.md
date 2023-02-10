# portfolio

Welcome !     
This directory provides select examples of code I have written in the last few years. This covers manipulating of small and large datasets, deriving diverse physics-related parameters, conducting data analysis, visualizing results in a clear and eye-appealing fashion, as well as producing maps and manipulating several spatial objects. Most of the code is written in R, but sequences in Octave (~Matlab) and SAC (Seismic Analysis Code) are also available.

+ **SAC**  
Earthquake signals recorded by seismic monitoring stations are frequently made available in so-called SAC format, requiring the use of Seismic Analysis Code (SAC) for their handling. Using RAW seismic data from stations operated during a recent earthquake, this work shows code routines for reading, identifying, treating, plotting and exporting suitably processed SAC signals and could be easily modified to handle additional settings (station number, filtering, different earthquakes).
  

+ **Seismic_records**   
This project is the next logical step upon pre-processing the SAC signals. Such SAC signals are read into R code in order to extract common engineering earthquake parameters such as PGA and PGV or to derive more complex parameter such as spectral accelerations over a range of structural periods and Arias intensity, with relevant plots illustrating the entire sequence.
This routine can directly follow from the SAC-based codes presented in directory /SAC, or be applied to other (already processed) SAC files.


+ **NDSHA**    
In this project, a dataset of synthetic earthquake seismograms generated by the NDSHA methology are read into and processed using R code, in order to extract and derive selected earthquake engineering parameters. The code shows how to incorporate the large dataset by means of a tree structure that allows to keep track of the several relevant simulation settings and parameters.


+ **Residuals**       
This section is the final stage upon obtaining engineering parameters from seismic records (=observations) and from predictions (=simulated synthetic signals). The large prediction dataset is validated against observations using a residual analysis and then relying on various ggplot-based visualisations, which shows how performing the predictive model behaves. Plots are also utilised to assess the impact of other parameters/settings involved in deriving predictions and observations.


+ **Kriging**    
Using earthquake acceleration data defined at the nodes of a regularly-spaced grid covering the United Kingdom, this project shows the steps leading to smoothing of said data over the land by means of kriging technique. This work combines dataframe reshaping and spatial objects (polygons, dataframe ...). The purpose is to obtain a higher resolution national hazard map covering the entire territory.
    
+ **Spectra**    
A code written in GNU Octave to compute, plot and save spectral values (acceleration, velocity and displacement) based on an externally-read earthquake signal. This code served as a basis when ported into R when treating earthquake signals at scale (see Seismic_records or NDSHA folders).

