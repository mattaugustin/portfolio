In this project, the code computes spectral values (acceleration, velocity and displacement) based on an earthquake acceleration signal read from an external file, and plots said spectral values as well as saves them in external pdf and jpg files.

The code contains several options for the user to comment/uncomment (set of periods, plotting, saving...).
Two plotting options are suggested, a "basic" one where only spectral accelerations are displayed, and an "enhanced" version where spectral accelerations, velocities and displacement are all plotted along with cosmetic plot improvements and then saved into jpg and pdf format.

I "ported" the function for spectral value computation into R code, as can be seen in portfolio/NDSHA/XX.R  or portfolio/seismic_records/XX.R , in order to "keep all operations under the same roof"
The Octave code however costs less in time and computer-resources, owing to Octave performance in matrix computing, and is therefore a great tool when computing spectral values for a large number of periods (typically more than two dozens periods based on user experience).


> Spectral amplitudes when using "enhanced" plotting option

![response_spectra](https://user-images.githubusercontent.com/61290423/218084123-997fb4d3-ff8f-43f0-aa08-ead55717604c.jpg)
