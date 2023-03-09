In high resolution spectrometry, one main challenge is to separate a beam of near-identical mass particles into several beams each composed of only one mass and hence one single particle type, thus ensuring "purity" of each beam for ensuing relevant experiments.
As an example, an optical setup  (made of multipoles, beam drifting tubes, magnet separators) as pictured in Figure 1 is able to separate particles with mass M1=100 and M2=100.0001 (i.e. a $\delta$ = 0.0001) respectively. This means the setup has a resolution R = 1 / $\delta$ = 10 000.

In evaluating the performance of an optical setup, particle beam transport modelling software (i.e. SIMION, CosyScript ...) commonly show beam trajectories (Figure 1), calculate beam transformations (compared to the original beam profile) along the optical setup, model complex optical abberations inherent to the setup.

> Figure 1. Left: beamline layout example (3D Autocad). Right: corresponding beam trajectory modelling in the horizontal plane (modelled in CosyScript), where the beam is represented using a few rays and beam envelope (outermost ray).     
      
![layoutoptic2](https://user-images.githubusercontent.com/126000617/220473800-5e5050ef-093c-494b-91ca-4d4fc37dbc20.png)


Particle transport modelling software do not however provide a visually-clear insight on the separation ability of the optical setup for different particle mass beams, since only one beam can be simulated in a transport and optimization routine. This gap is filled by means of external software, such as with this project Octave code. The code output shows whether, and if so, how well, a multimass particle beam can be separated into single mass smaller beams.

       
Each particle of the beam is characterized by its position x and y (usually expressed in mm) in the horizontal and vertical plane respectively, but also by its angular orientation (in mrad) in each plane that indicate "where" in the plane the particle is heading. Together, all particles form a beam whose boundaries are thus defined by the outermost particle position (Figure 2), and the beam boundaries are most often referred to as beam emittance. A beam with boundaries of 0.3 mm and 10 mrad thus has an emittance 3 $\pi$.mm.mrad .    


> Figure 2. Beam initial definition plotted in both horizontal and vertical planes.  

![beamstart](https://user-images.githubusercontent.com/126000617/220470657-c73feaa1-bce4-4610-9f18-399615f7f661.png)


The COSY software simulated the beam envelope particle trajectory through the setup for a reference mass M=100, running optimization routines on few setup parameters (voltages, distances...). As the beam undergoes transformations resulting from the setup configuration, the user can select locations at which the software "records" such transformation in a matrix, the coefficients of which could be applied to any particle from the starting beam to calculate its position in the "recording" plane. Besides the transformation coefficients, the actual particle characteristics can also be recorded by the software (Figure 3).       
Beam characteristics in any recording plane "n" can be calculated using the transformation matrix M of the relevant plane as per Figure 4 formula. Coefficients of each column of M relate to one of the beam characteristics, and each line is associated to different order effects. As an example, for a particle defined in plane 0, the horizontal position x in plane n=1 is expressed using matrix coefficients in column 1 and beam order coefficients in column 6, as follows:

x(n) =    $x_{1}^{1}$ . $a^{0}_{0}.y^{0}_{0}.b^{0}_{0}.t^{0}_{0}.m^{0}_{0}$     ( 1st row)       
       +  0$x^{0}_{0}$.a^{1}_{0}.y^{0}_{0}.b^{0}_{0}.t^{0}_{0}.m^{0}_{0}$            
       +  $0.x^{0}_{0}.a^{0}_{0}.y^{1}_{0}.b^{0}_{0}.t^{0}_{0}.m^{0}_{0}$     ( 3rd row)       
       +  $0.x^{0}_{0}.a^{0}_{0}.y^{0}_{0}.b^{1}_{0}.t^{0}_{0}.m^{0}_{0}$     ( 4th row)      
       +  $0.x^{0}_{0}.a^{1}_{0}.y^{0}_{0}.b^{0}_{0}.t^{1}_{0}.m^{0}_{0}$     ( 5th row)      
       +  $0.x^{0}_{0}.a^{1}_{0}.y^{0}_{0}.b^{0}_{0}.t^{0}_{0}.m^{1}_{0}$     ( 6th row)     

$0.x^{0}_{0}.a^{1}_{0}.y^{0}_{0}.b^{0}_{0}.t^{0}_{0}.m^{0}_{0}$

> Figure 3. Example of beam characteristics and beam transformation matrices, as recorded by CosyScript.

![transfomatrix](https://user-images.githubusercontent.com/126000617/224098629-1b9227f3-21d4-43bd-9b0c-3e45ca2af1db.PNG)


As the beam travels through the optical setup, particles behave differently depending on their respective mass and hence see their flying path altered. Using the transfer matrices obtained when simulating the beamline for the reference mass M=100, a beam containing 50 000 particles of near-identical masses ( 99.9999, 100, 100.0001) was input into the transfer matrices defined in all recording planes in order to calculate the respective beam projections in both horizontal and vertical planes. In effect, components of each particle (positions x, y and plane orientation angles a, b, time of flight t and mass difference $\delta$) are multiplied by the matrix coefficients according to the relation (1) to obtain the coordinates of the projected particle in the selected plane.         
Plotting the projected particles in all planes, in particular the last one reflecting the impact of the entire beamline, reveals how well the optical setup performs in terms of separation ability, i.e. the resolution of the optical setup.

          
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


