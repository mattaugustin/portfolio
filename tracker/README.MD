In high resolution spectrometry, one main challenge is to separate a beam of near-identical mass particles into several beams each composed of only one mass and hence one single particle type, thus ensuring "purity" of each beam for ensuing relevant experiments.
As an example, an optical setup  (made of multipoles, beam drifting tubes, magnet separators) as pictured in Figure 1 is able to separate particles with mass M1=100 and M2=100.0001 (i.e. a $\delta$ = 0.0001) respectively. This means the setup has a resolution R = 1 / $\delta$ = 10 000.

In evaluating the performance of an optical setup, particle beam transport modelling software (i.e. SIMION, CosyScript ...) commonly show beam trajectories (Figure 1), calculate beam transformations (compared to the original beam profile) along the optical setup, model complex optical abberations inherent to the setup, but do not provide a visually-clear insight on the separation ability of the optical setup for different particle mass beams. 

> Figure 1. Left: beamline layout example (3D Autocad). Right: corresponding beam trajectory modelling in the horizontal plane (modelled in CosyScript).     
      
![layoutoptic2](https://user-images.githubusercontent.com/126000617/220473800-5e5050ef-093c-494b-91ca-4d4fc37dbc20.png)


This gap is filled by means of external software, such as with this Octave code. The code output shows whether, and if so, how well, a multimass particle beam can be separated into single mass smaller beams. Using the beam transformation data (i.e. transfer matrices) recorded across the optical setup, the Octave routine projects and displays the beams spatial distribution in all planes where beam characteristics were recorded (via transfer matrices). 
       
Each particle of the beam is characterized by its position x and y (usually expressed in mm) in the horizontal and vertical plane respectively, but also by its angular orientation (in mrad) in each plane that indicate "where" in the plane the particle is heading. Together, all particles form a beam whose boundaries are thus defined by the outermost particle position (Figure 2), and the beam boundaries are most often referred to as beam emittance. A beam with boundaries of 0.3 mm and 10 mrad thus has an emittance 3 $\pi.mm.mrad .
As the beam travels through the optical setup, particles behave differently based on their respective mass. By "recording" the beam transfer matrix at selected places (=planes), aka the coefficients whose linear combination allow to calculate the position of a given particle, the user can determine the projected beam of particles in these planes.      


> Figure 2. Beam initial definition plotted in both horizontal and vertical planes.  

![beamstart](https://user-images.githubusercontent.com/126000617/220470657-c73feaa1-bce4-4610-9f18-399615f7f661.png)

      
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


