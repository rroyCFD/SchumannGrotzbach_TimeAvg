## Application : SchumannGrotzbach_TimeAvg (Wall shear-stress boundary condition for inhomogeneous surface)

### Author:
- Rajib Roy
- University of Wyoming
- rroy@uwyo.edu, roy.rajib@live.com

### Description:

Wall boundary condition for total shear stress at a solid (impermiable) boundary. Stresses are computed using the Schumann/Grotzbach formulation. In an inhomogeneous surface, average velocity is obtained though a time average instead of a horizontal-average or local mode.

The Boundary Condition (BC) is developed based on the [SchumannGrotzbach](https://github.com/NREL/SOWFA/tree/master/src/finiteVolume/fields/fvPatchFields/derived/surfaceStressModels/SchumannGrotzbach) BC implemented in the [SOWFA](https://github.com/NREL/SOWFA) toolbox for OpenFoam developed at [NREL](https://nwtc.nrel.gov/SOWFA).

**Usage:**
```
an_inhomogeneous_boundary
    {
        type            SchumannGrotzbach_TimeAvg;
        // default parameter values
        kappa           0.4; 
        betaM           16;
        gammaM          5;
        z0              uniform 0.02; // surface roughness length
        averageType     "timeAveraged";
        // place-holder for surface-parallel velocity
        U1_Paral_local_timeAve     uniform (0 0 0); 
        restart         true;
        TimeRange       2500; // time averaging window
        value           uniform (0 0 0 0 0 0); // place-holder for the stress tensor
    }
```

### Disclaiimer:

This application is built based on [OpenFOAM version-2.1.x](https://openfoam.org/release/2-1-0/). Please read the _About OpenFOAM_ section to learn more on OpenFOAM.

The application is free to use. The author neither provide any warranty nor shall be liable for any damage incurred from this application.



#### About OpenFOAM

OpenFOAM is the leading free, open source software for computational fluid dynamics (CFD), owned by the OpenFOAM Foundation and distributed exclusively under the [General Public Licence (GPL)](http://www.gnu.org/copyleft/gpl.html). The GPL gives users the freedom to modify and redistribute the software and a guarantee of continued free use, within the terms of the licence. To learn more visit [https://openfoam.org/](https://openfoam.org/)
