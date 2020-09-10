# 3D-GEB-mixed-mode-model
This project develops a novel 3D mixed-mode model to predict kinetics of austenite-ferrite phase transformations in low-alloyed steels under various conditions of heat treatments. This model takes into account both nucleation and interface migration. The nucleation is calculated according to continuous nucleation theory. The interface migration is computed by deriving interface velocity that meets the Gibbs energy balance (GEB) beteween chemical driving force and energy dissipations due to interface friction and trans-diffusion of substitutional elements inside the interface. The latter is normally termed as solute drag. Since an analytical solution for calculating the solute drag is proposed, the computation speed has been increased significantly compared to the conventional approach. New features of this model are its efficient algorithm to compute energy dissipation by solute drag, its capabilities of predicting the microstructural state for spatially resolved grains and the minimal adjustment of modelling parameters.

# Features of the model
- austenite-ferrite phase transformations during continuous cooling, isothermal holding and thermal cycling in the two-phase region.
- Nucleation and interface migrates isotropically
- Gibbs energy balance between chemical driving force and energy dissipations due to interface friction and solute drag 

# Denpendencies of the code
1. Install the Multi-Parametric Toolbox 3 (mpt3): https://www.mpt3.org/ into the same folder, for generating Voronoi cells to represent austenite grains.
2. Make sure that 'Optimization toolbox' and 'Symbolic Math toolbox' are included in your own Matlab package. These toolboxes should be included by default.
All codes have been tested executable with Matlab 2014b or above.

# Running the simulation
Very simple. Just run the file 'ferrite_3d_model_voronoin_PBC_ND_CNT_GEB.m'.
But, remember to first set up the chemical compositions, heat treatment parameters and including key thermodynamic data in the file 'SimulCond.m'.
Set the 'CyclicFlag' for different thermal routes:
1-thermal cycling;
0-isothermal holding;
-1-continuous cooling.

# Run on a linux cluster
The code can run on a linux cluster by running the bash file 'FeCMn_3D_GEB.pbs'.

# License
This package is free to use, ditribute and adapt, but no warranty and liability to any kinds of simulations.
Citing our work (will be published and updated here) is strongly encouraged if you use or get inspired by our code.
See the __LICENSE__ file for license rights and limitations (GNU General Public License v3.0).

## Contact via hfang@tudelft.nl or haixingfang868@gmail.com

