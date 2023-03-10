# portfolio

Welcome !     
This directory provides select examples of code I have written in the last few years. This covers manipulating of small and large datasets, deriving diverse physics-related parameters, conducting data analysis, visualizing results in a clear and eye-appealing fashion, as well as producing maps and manipulating several spatial objects. Most of the code is written in R, but sequences in Octave (~Matlab) and SAC (Seismic Analysis Code) are also available.


### Kriging (R code)   
The purpose is to obtain a higher resolution national hazard map covering the entire territory. Using earthquake acceleration data defined at the nodes of a regularly-spaced grid covering the United Kingdom, this project shows the steps leading to smoothing of said data over the land by means of kriging technique. This work combines reshaping dataframes and manipulating spatial objects (polygons, dataframe ...).
    
### Spectra (Octave)          
A code written to compute, plot and save earthquake response spectral values (acceleration, velocity and displacement) associated to an earthquake signal. This code served as a basis when ported into R when treating earthquake signals at scale (see "Seismic_records" or "NDSHA" folders), but offers faster performance when working with many structural periods.

### Tracker (Octave)            
A tracker which calculates and plots the projected position of beam particles across several user-defined locations along a beamline, ultimately providing insights on the performances of a high resolution mass spectrometry setup.


### Predictions vs Observations   
_**Context:**_     
Seismologists and engineers are most often in need of earthquake ground motion data (the seismograms recorded when an earthquake occurs) for use in various fields (prediction equations, structural analysis, damage mitigation strategies...). Ideally, for the sake of specificity, these ground motion data are available for many geographical regions (ex: UK), tectonic environment (seismically active, moderate or low), type of soils (rock, hard and soft soils...), across a range of magnitude and distance from the earthquake locus, etc ... . Unfortunately for science, such amount and diversity of data is not always available, as it depends on earthquake occurrence, since strong earthquakes (magnitude above 6) rarely happen and only in limited parts of the globe.      
One field of earthquake engineering with growing interest and addressing this issue is so-called "earthquake ground motion synthesis", whereby artificial but realistic earthquake seismograms are generated by means of a physics model involving the knowledge of seismic fault rupture and seismic wave propagation. The rationale is: rather than "waiting" for strong but rare earthquakes to obtain ground motion data, how about artifically "generating" these signals using physics-based models? This would be ideal, for as long as the artificially generated signals can faithfully enough reproduce the signals recorded during actual earthquakes. In other words, for a given model, how well do artificial signals compare to real signals ? And how to compare ?         
This 4-part project aims at assessing the performance of the NDSHA, one model that generates artificial earthquake signals tailored to a region of interest (ex: UK), relying on the Modal Summation method. This technique was chosen among others owing to the little amount of geophysical & seismological data required in the model implementation. The earthquake signals were generated by yours truly, using Fortran codes proprietary of the Geophysics Dept of the University of Trieste, Italy.
The comparison between predictions and observations relied on several engineering parameters associated to earthquake-induced ground shaking, such as Peak Ground Acceleration or earthquake response spectral accelerations.        

+ **Predictions vs Observations I : SAC**  (SAC code)      
Earthquake signals recorded by seismic monitoring stations are frequently made available in so-called SAC format, requiring the use of Seismic Analysis Code (SAC) for their handling. Using RAW seismic data from stations operated during a recent earthquake, this work shows code routines for reading, identifying, treating, plotting and exporting suitably processed SAC signals and could be easily modified to handle additional settings (station number, filtering, different earthquakes).
  

+ **Predictions vs Observations II : Seismic_records** (R code)       
This work comes upon pre-processing the SAC signals and is the completion step in obtaining the engineering parameters associated to earthquake observations. The SAC signals are read into R code in order to extract parameters such as PGA and PGV or to derive more complex parameter such as spectral accelerations over a range of structural periods as well as Arias intensity that reflects the seismic energy. Relevant plots illustrate parts of the sequence.
This routine can directly follow from the SAC-based codes presented in directory /SAC, or be applied to other (already processed) SAC files.


+ **Predictions vs Observations III : NDSHA** (R code)           
In this section, a large dataset of synthetic earthquake seismograms generated by the NDSHA methology are read into and processed using R code, in order to extract and derive selected earthquake engineering parameters later used for comparing against seismic observations. The code shows how to incorporate the large dataset by means of a tree structure that allows to keep track of the several model simulation settings (front rupture directivity, fault rupturing realisation...) and parameters (filtering frequencies, earthquake motion types...).


+ **Predictions vs Observations IV : Residuals**  (R code)              
This section is the final stage upon obtaining engineering parameters from seismic records (=observations) and from predictions (=simulated synthetic signals). The large prediction dataset is validated against observations using a residual analysis and then relying on various ggplot-based visualisations, which shows how performing the predictive model behaves. Plots are also utilised to assess the impact of other parameters/settings involved in deriving predictions and observations.
