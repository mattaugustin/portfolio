In high resolution spectrometry, one main challenge is to separate a beam of near-identical mass particles into several beams each composed of only one mass and hence one single particle type, thus ensuring "purity" of each beam for ensuing relevant experiments.
As an example, a optical setup (made of multipoles, beam drifting tubes, magnet separators) able to separate particles with mass M1=100 and M2=100.0001 (i.e. a $\delta$ = 0.0001) respectively would exhibit a corresponding resolution R = 1 / $\delta$ = 10 000.

Particle beam transport modelling software (i.e. SIMION, CosyScript ...) commonly show beam trajectories , track beam transformation coefficients at user-selected locations along the optical setup, can account higher order effects owing to real-world optics challenges, but do not provide a visually-clear insight on the separation ability of the optical setup for different particle mass beams. 

> beamline layout example (3D Autocad) and corresponding beam trajectory modelling in the horizontal plane (modelled in CosyScript).     
      
![layoutoptic2](https://user-images.githubusercontent.com/126000617/220473800-5e5050ef-093c-494b-91ca-4d4fc37dbc20.png)


This gap is filled by means of external software, such as with this Octave code. The code output shows whether, and if so, how well, a multimass particle beam can be separated into single mass smaller beams. Using the beam transformation data (i.e. transfer matrices) recorded across the optical setup, the Octave routine projects and displays the beams spatial distribution in all planes where beam characteristics were recorded (via transfer matrices).       
A beam is composed of N particles, each characterized by its size x (e.g mm) and angle orientation a (e.g. mrad) in the horizontal plane, and size y and angle b in the vertical plane respectively.


> beam initial definition plotted in both horizontal and vertical planes.    

![beamstart](https://user-images.githubusercontent.com/126000617/220470657-c73feaa1-bce4-4610-9f18-399615f7f661.png)

As the beam travels through the optical setup, particles behave differently based on their respective mass. By "recording" the beam transfer matrix at selected places (=planes), aka the coefficients whose linear combination allow to calculate the position of a given particle, the user can determine the projected beam of particles in
these planes.      
Each particle is indeed input to the location-relevant transfer matrix to calculate the parameters of the corresponding projected particle in both horizontal and vertical planes. Applying this operation to all particles of a beam means that projected beams can be plotted to assess how well the optical setup performs in terms of separation ability, i.e. the resolution of the optical setup.
       
          
The Octave routine includes :
+ generating a beam of N particles within the boundaries of the beam size,
+ calculating projected particles using the transfer matrix file provided by the beam transport modelling software,
+ creating and adequately parametrizing a histogram associated to the particle spatial distribution,
+ setting and plotting beam particles across **all** planes with common axes, using either beam horizontal or vertical component,
+ setting and plotting beam particles across **selected** planes with tailored axes, using either beam horizontal or vertical component.



> Three beams of nearly-identical mass particles projected across all observation planes, along with a corresponding particle count histogram. Axes are common to all planes.       
      
![partrack_enh_horiplane](https://user-images.githubusercontent.com/126000617/222761048-fbfbc6c2-f85b-40b5-a097-b44310533fd7.png)


> Beams plotted in selected planes worthy of interest to the user. Axes are suited to each individual plane.      
      
![partrack_enh_horiplane_selected](https://user-images.githubusercontent.com/126000617/222768657-dc6456c4-f714-4b37-9114-f8a1eca77f7f.png)


