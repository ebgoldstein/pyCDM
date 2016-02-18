"""
STEP 1 in the CDM
    Compute the shear stress field. The shear stress field is computed by
    using the u* defined for each time iteration and the surface topography.
    The surface topography sets the shear stress through a linearized RANS
    formulation (Weng et al. 1991). This is an analytical solution utilizing
    a fourier transform. The topographic surface is
    modified by the existence of a seperation bubble when there is a dune
    brink. The resultant shear stress, from the linearized RANS formulation,
    is then modified by the vegetation field

"""
import numpy as np
from CDMparams import *
import CDMfunctions as cdmfns

def shearfield(Topo, Veg, Tau)

    #if brink (slope of 20 deg or greater; Duran and Moore 2013)
        #separation bubble
        #make new combined surface
    #else
        # Topowind =Topo

    #calculate tau pert

    #calculate combined tau

    #calculate reduction factor for Vegetation

    return Tau
