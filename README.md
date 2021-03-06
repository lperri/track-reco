# track-reco
Project related to CMS Particle Track Reconstruction

**Experiments the Large Hadron Collider (LHC) at CERN:**
The particle collider consists of various experiments, including the Compact Muon Solenoid (CMS), A Large Ion Collider Experiment (ALICE), A Toroidal LHC Apparatus (ATLAS). For a visual representation of the layout of experiments with respect to the LHC, click [here](http://cds.cern.ch/images/OPEN-PHO-ACCEL-2013-056-1)).

**Scientific goals of the LHC:**
One of the reasons why the LHC was built was to discover the [Higgs Boson](https://en.wikipedia.org/wiki/Higgs_boson). This goal was accomplished in 2012. However, there is more physics to be done. The Higgs has been found, but further experiments are required to confirm that the Higgs we found is the same one that was predicted by theorists. Additionally, measurements of its properties remain important in helping us uncover physics beyond the [Standard Model](https://home.cern/science/physics/standard-model) of particle physics. The Standard Model represents humanity's best, yet incomplete, understanding of the building blocks (elementary particles) making up our universe. Searching for physics Beyond the Standard Model (BSM) is another goal of the LHC.

**Physics knowledge required to understand my project:**
Charged particles leave "tracks" in magnetic fields. CMS contains a 4T magnet, which is responsible for bending the tracks of chargedparticles born in proton-proton collisions (the beams used at LHC are proton beams). A simple way to visualize a spary of particles resulting from collisions is a car crash: when two cars collide head-first, tons of pieces of shrapnel fly out. A useful physical quantity that tells us a lot about a particle is momentum. Luckily, there is a fundamental relationship between the amount of momentum a particle possesses and how much it bends in a magnetic field (the more momentum the particle has, the harder it is for the magnetic field to bend it). 

**Basics of muon reconstruction in CMS**
Collisions in the detector happen very, very quickly. Therefore, we need a way to slow things down in order to analyze what actually happened. "Reconstruction" refers to looking at the aftermath of an event, stored as digitized data, and draw a physical picture of what particles resulted from the proton-proton collision.

Protons are composite particles and, as a reuslt, their collisions produce sprays of all kinds of particles that we can detect. In addition, heavy particles decay into lighter particles. One such light particle that CMS was built to detect is the [muon](https://en.wikipedia.org/wiki/Muon). We have known of the existence of muons for almost a century, however, muons are important to the particle reconstruction process for technical reasons (in a nutshell, they provide us with a very clean signal, which aides in the discovery of more exotic particles). 

Muons don't interact with the matter placed on the inner layers of the detector, so they are found in the endcaps of the [detector](http://cms.web.cern.ch/news/muon-detectors).
The muon reconstruction chain starts with the "local reconstruction". First, hits in DTs, CSCs and RPCs are reconstructed from digitized electronics signals. Hits within each DT and CSC chamber are then matched to form "segments" (track stubs). More information can be found [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMuonAnalysis).

**Physics variables glossary:**
- 'beam-line': proton beams travel along the z-axis of {CMS coordinate system}
- 'interaction point': location of collision along beam-line
- 'eta': stands for [pseudorapidity](https://en.wikipedia.org/wiki/Pseudorapidity), a quantity often used in place of 'theta,' the polar angle in a spherical coordinate system, within the field of experimental high energy physics; eta is an angular measure of how parallel a particle's trajectory is to the beam-line (beam-line is at eta=0, completely orthogonal to beam-line is eta=infinity)
- 'phi': the azimuthal angle in a cylindrical coordinate system
- 'k': in particle physics, lowercase k or the greek letter kappa refers to a quantity called "curvature" which is a measure of how much something bends
- 'pt': in particle physics, refers to transverse momentum, a 

**CMS-specific glossary:**
- 'event': an event is [recorded data] from a collision
- 'the barrel': the main portion of the detector (does not include endcaps), containing DTs and RPCs (see definitions below)
- 'DT': Drift Tube: structures that are responsible for measuring muon's position when it is in the 'barrel'
- 'RPC': Resistive Plate Chamber: structures responsible for giving quick measurements of muon momentum  
- 'CSC': Cathode Strip Chambers: trapezoidal chambers that make up the CMS Endcap Muon system
- 'Geant': shorthand for [Geant4](http://geant4.web.cern.ch/), which is simulation software developed by CERN for the purpose of simulating; any variable or function containting the name 'gen' refers to simulated data generated with this software
- a 'hit'/'stub': trigger primitives, software objects that indicate charge has built up on one of the wires, indicating that a particle has crossed that location on its trajectory
- 'wire': wires run azimuthally (in the phi direction) along panels in the CSCs and define a hit's radial coordinate
- 'strip': strips run radially (in the r direction) along panels in the CSCs and define a hit's angular coordinate
- 'sector': 60 degree portions of the muon system 

**About This Repository**
Guide to using the software:
First build the [Work Area](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookSetComputerNode#Create_a_work_area_and_build_the).
Here is a brief description of the functionality of each directory and its contents:
- **bash_scripts**: running setup_env.sh in working directory will set up CMS environment
- **non_digital_units**: digital units here refers to whether or not the quantities being manipulated are raw detector output (bits) or non-integer quantities, for example, an angle measured in radians. This project was originally done in non-digital units (i.e. not with raw detector data), but later was put into digital units so that it could accomodate raw data.
- **old_scratch_not_used**: self explanatory... just keeping it for my records.
- **angleManipulations.py**: makes sure that angles are in the correct range (and if they are not, it moves them to the correct range); also converts angles measured in radians to raw detector ("digital") units
- **associatedStubsAlgorithm.py**: a muon object (representing a real muon) in this project is defined by a few attributes (angles of trajectory and curvature) as well as the number of associated stubs. If the object has more than 2 associated stubs (which are determined to be associated based on "distance" from each other in "phi-space"), then it is considered to be a real muon.
- **checkRange_thetaFP.py**: a script that outputs a 1d histogram of all theta_fp values, purely for debugging purposes; thetaFP stands for "fixed point theta", theta being the polar angle in a spherical coordinate system and "fixed point" referring to it being measured in bits -- raw detector data.
- **makeHistograms.py**: uses PyROOT to set up histgram objects, later to be filled by fillHistograms.py.
- **fillHistograms.py**: fills histograms for the purpose of calibrating angles phi & theta; takes in all information from retrieveDataAsEvents.py, makeHistograms.py and angleManipulations.py. 
- **fitHistograms.py**: if you want to fit the histgorams filled by fillHistograms.py, this script fits the 2D calibration histograms (projected into 1D profiles) to the proper functions (based on physics concepts)
- **generate3Dhistograms.py**: generates 3D histograms (x-axis=k,y-axis=difference in phi,z-axis=thetaFP); the point of these histograms is to check if our phi calibration has any dependence on theta, i.e. if the difference in phi depends on where we are in "theta-space"; it also has funcitonality to write projections of the 3D histogram to the root file.
- **retrieveDataAsEvents.py**: specifies which tag you want to use and then cleans that data prior to use in analysis.
<!---Calibration, Propagation, Algorithm, Efficiency // Digital & Not Digital-->
