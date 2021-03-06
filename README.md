<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<h1 align="center">
    <a>Development of fracture critical material</a>
</h1>

<h4 align="center">
    <a href="#material-overview">Material Overview</a> |
     <a href="#material-compiling">Material Compiling</a> |
    <a href="#connection-model">Connection model</a> |
    <a href="#steel-frame-performance">Steel Frame Performance</a>
</h4>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->

This repo contains the descriptions of the fracture critical material on [OpenSees](https://opensees.berkeley.edu/wiki/index.php/Main_Page) and how to compile it onto Stanford Cluster System ([Sherlock](https://www.sherlock.stanford.edu/)) for parallel computing. The use of the material is further described in the [connection model](#connection-model) and [steel frame performance](#steel-frame-performance) sections.

## Material Overview
In the proposed method, the connection is modeled using fiber beam elements, in which the hinge is discretized into fracture critical flanges and web fibers. In this model, the flange fiber constitutive (i.e. stress-strain) relationship represents the series of the weld and beam-hinging region of the connection. The flange fiber loses tensile resistance when fracture occurs, while the behavior under compressive forces is unaffected considering the re-contacting of the beam flange and column face. After connection fracture, the compressive resistance acts like a gap material, resembling the open and closing of the gap between the beam flange and column face.

- Material Constitutive
<p align="center"><img style="max-width:500px" width="640" src="https://github.com/wenyiyen/Fracture-critical-material/blob/master/fig/constitutive.png" alt="constitutive"></p>

- The material described in this study also includes a stress-based damage index function, which combines linear S-N curve, modified cycle counting method, and the Miner's rule.
    
    - <img src="https://latex.codecogs.com/svg.latex?\Large&space;\sigma_{cr}=90,\;m=0.275" title="\Large \sigma_{cr}=90,\;m=0.275" />
    <p align="center"><img style="max-width:500px" width="640" src="https://github.com/wenyiyen/Fracture-critical-material/blob/master/fig/DI.gif" alt="DI"></p>


## Material Compiling
The source code of the material is added in `./OpenSeesMP_Sherlock/SRC/material/uniaxial/`: [SteelFractureDI.h](https://github.com/wenyiyen/Fracture-critical-material/blob/master/OpenSeesMP_Sherlock/SRC/material/uniaxial/SteelFractureDI.h) and [SteelFractureDI.cpp](https://github.com/wenyiyen/Fracture-critical-material/blob/master/OpenSeesMP_Sherlock/SRC/material/uniaxial/SteelFractureDI.cpp).

### Compiling OpenSeesMP on Sherlock (user local)
- **Set up**\
Copy the full directory [OpenSeesMP_Sherlock](https://github.com/wenyiyen/Fracture-critical-material/tree/master/OpenSeesMP_Sherlock) to <code>/home/users/**_user_**/</code> on [Sherlock](https://www.sherlock.stanford.edu/).\
Install [tcl](https://www.tcl.tk/software/tcltk/) to <code>/home/users/**_user_**/tcl/</code>.\
    Proceed with Sherlock terminal:
<pre>
mkdir /home/users/<b><i>user</b></i>/lib
mkdir /home/users/<b><i>user</b></i>/bin
cd /home/users/<b><i>user</b></i>/OpenSeesMP_Sherlock
</pre>

- **Delete previous make items**
```
make clean
make wipe
```
- **Load needed modules for parallel computing**
```
module load scalapack
module load mumps
module load metis
module load parmetis
module load petsc
module load python # (when building python)
```
- **Compile**
```
make
```
The OpenSeesMP binary will be built in <code>/home/users/<b><i>user</b></i>/bin</code>.
### Adding New Material into OpenSeesMP on Sherlock - *Optional*
If interested in adding newly developed material into OpenSees, please refer to instructions below.
- **Include material source code**
    - Add files *newMaterial*.h and *newMaterial*.cpp into [/OpenSeesMP_Sherlock/SRC/material/uniaxial/](https://github.com/wenyiyen/Fracture-critical-material/tree/master/OpenSeesMP_Sherlock/SRC/material/uniaxial).
- **Create internal links for OpenSees to recognize the new material**
    - Add lines in [OpenSeesUniaxialMaterialCommands.cpp](https://github.com/wenyiyen/Fracture-critical-material/blob/master/OpenSeesMP_Sherlock/SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp)
    - Add lines in [TclModelBuilderUniaxialMaterialCommand.cpp](https://github.com/wenyiyen/Fracture-critical-material/blob/master/OpenSeesMP_Sherlock/SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp)
    - Add lines in [classTags.h](https://github.com/wenyiyen/Fracture-critical-material/blob/master/OpenSeesMP_Sherlock/SRC/classTags.h)
- **Editing makefiles for OpenSeesMP**
    - Add line <code>*newMaterial*.o</code> in [OpenSeesMP_Sherlock/SRC/material/uniaxial/Makefile](er/OpenSeesMP_Sherlock/SRC/material/uniaxial/Makefile)
    - Add line <code>*newMaterial*.o</code> in [OpenSeesMP_Sherlock/SRC/Makefile](https://github.com/wenyiyen/Fracture-critical-material/blob/master/OpenSeesMP_Sherlock/SRC/Makefile)
- **Compile**
    - Follow section ["Compiling OpenSeesMP on Sherlock (user local)"](#compiling-openseesmp-on-sherlock-user-local)

## Connection Model
To resemble the lab test beam-column connection assembly, we model the essential components as the following:
- Beam: A force-based element with hinge sections discretized into multiple fibers. The two flange fibers (top and bottom flange fibers) are modeled using the fracture critical material described in the [first section](#material-overview).
- Column: Elastic elements.
- Panel zone: Krawinkler's model.
- Below figure shows the lab test data vs the simulation result:
<p align="center"><img style="max-width:500px" width="640" src="https://github.com/wenyiyen/Fracture-critical-material/blob/master/fig/SAC.gif" alt="SAC"></p>

## Steel Frame Performance
The study includes the performance analysis of a 31 story steel moment frame building.
- For multi-stripe analysis, please refer to [MRF31_Multi-stripe](https://github.com/wenyiyen/Fracture-critical-material/tree/master/MRF31_Multi-stripe).
- For incremental dynamic analysis, please refer to [MRF31_IDA](https://github.com/wenyiyen/Fracture-critical-material/tree/master/MRF31_IDA).

Plots below shows the basic plans and geometry for the building, followed by the incremental dynamic analysis results and the collapse fragility curve.

<p align="center"><img style="max-width:500px" width="600" src="https://github.com/wenyiyen/Fracture-critical-material/blob/master/fig/plan.png" alt="plan"></p>
<p align="center"><img style="max-width:500px" width="800" src="https://github.com/wenyiyen/Fracture-critical-material/blob/master/fig/IDA.png" alt="IDA"></p>


