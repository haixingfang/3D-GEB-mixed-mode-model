# [3D-GEB-mixed-mode-model](https://github.com/haixingfang/3D-GEB-mixed-mode-model)
## Austenite-ferrite TransModel 2.0
The code were developed for predicting kinetics of austenite-ferrite phase transformations in low-alloyed steels under various conditions of heat treatments. This model takes into account both nucleation and interface migration. The nucleation is calculated according to continuous nucleation theory. The interface migration is computed by deriving interface velocity that meets the Gibbs energy balance (GEB) beteween chemical driving force and energy dissipations due to interface friction and trans-diffusion of substitutional elements inside the interface. The latter is normally termed as solute drag. Since an analytical solution for calculating the solute drag is proposed, the computation speed has been increased significantly compared to the conventional approach. New features of this model are its efficient algorithm to compute energy dissipation by solute drag, its capabilities of predicting the microstructural state for spatially resolved grains and the minimal adjustment of modelling parameters. The codes were a result of a follow-up research to Haixing Fang's [PhD thesis](https://repository.tudelft.nl/islandora/object/uuid%3Aecd8e101-3164-4227-b47b-13a04bc4b8fb?collection=research) supervised by [Dr.ir. N.H. van Dijk](https://www.tudelft.nl/en/faculty-of-applied-sciences/about-faculty/departments/radiation-science-technology/research/research-groups/fundamental-aspects-of-materials-and-energy/people/niels-van-dijk/) and [Prof.dr.ir. S. van der Zwaag](https://www.tudelft.nl/lr/organisatie/afdelingen/aerospace-structures-and-materials/novel-aerospace-materials/people/personal-pages-novam/academic-staff/s-van-der-zwaag-sybrand/) at Delft University of Technology.

# Features of the model
- Austenite-ferrite phase transformations during continuous cooling, isothermal holding and thermal cycling in the two-phase region.
- Nucleation and interface migrates isotropically
- Gibbs energy balance between chemical driving force and energy dissipations due to interface friction and solute drag 

# Denpendencies of the code
- Install the Multi-Parametric Toolbox 3 (mpt3): https://www.mpt3.org/ into the same folder, for generating Voronoi cells to represent austenite grains.
- Make sure that 'Optimization toolbox' and 'Symbolic Math toolbox' are included in your own Matlab package. These toolboxes should be included by default.
All codes have been tested executable with Matlab 2014b or above.

# Run the code on PC
Very simple. Just run [ferrite_3d_model_voronoin_PBC_ND_CNT_GEB.m](https://github.com/haixingfang/3D-GEB-mixed-mode-model/blob/master/ferrite_3d_model_voronoin_PBC_ND_CNT_GEB.m). <br>
But, remember to first set up the chemical compositions, heat treatment parameters and include thermodynamic data (which can be parameterized first with [Thermo-Calc software](https://www.thermocalc.com/)) in [SimulCond.m](https://github.com/haixingfang/3D-GEB-mixed-mode-model/blob/master/SimulCond.m). <br>
<br>
Within this file, set the variable [CyclicFlag](https://github.com/haixingfang/3D-GEB-mixed-mode-model/blob/master/SimulCond.m) for different thermal routes:
- 1-thermal cycling; <br>
- 0-isothermal holding; <br>
- -1-continuous cooling. <br>

# Run the code on a linux cluster
The code can run on a linux cluster by running the bash file [FeCMn_3D_GEB.pbs](https://github.com/haixingfang/3D-GEB-mixed-mode-model/blob/master/FeCMn_3D_GEB.pbs).

# License
This package is free to use, ditribute and adapt, but no warranty and liability to any kinds of simulation results. <br>
Citing our work (will be published and updated here) is strongly encouraged if you use or get inspired by our code. <br>
See the [LICENSE](https://github.com/haixingfang/3D-GEB-mixed-mode-model/blob/master/LICENSE) for license rights and limitations (GNU General Public License v3.0).

## Contact via hfang@tudelft.nl or haixingfang868@gmail.com

