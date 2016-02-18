"""
rewrite of CDM by EBG 2/16

This is the function which computes the actual change in the surface
profile. The CDM operates by iterating over 5 operations using 2 state variables:

The State Variables are:
 	-'Topo' (the topographic domain)
	-'Veg' (the vegetation field)

The 5 operations are:
	1) Compute the shear stress field. The shear stress field is computed by
		using the u* defined for each time iteration and the surface topography.
		The surface topography sets the shear stress through a linearized RANS
		formulation (Weng et al. 1991). This is an analytical solution utilizing
		a fourier transform. The topographic surface is
		modified by the existence of a seperation bubble when there is a dune
		brink. The resultant shear stress, from the linearized RANS formulation,
		is then modified by the vegetation field
	2) A sand flux field is computed from the shear stress field. Flux in the
		separation zone is set to 0.
	3) The flux divergence is computed for a given time interval (forward Euler)
		and this sets the change in topography
	4) The new topo field is polled for slopes greater than the angle of repose
		(34 degrees). For slopes greater than this, avalanches occur down the
		steepest descent gradient
	5) The vegetation field is modified based on growth equations, which
		encompass the accretion/erosion rate of the sand surface.

"""

import numpy as np
from CDMparams import *
import CDMfunctions as cdm


def step_implementation(Veg,Topo,Tau,ustar):

	#1. Compute the stress field

	#2. Compute the flux field

	#3. Compute flux divergence and update topo

	#4. Poll the grid for avalanches

    #5. Compute the vegetation field
    Veg = cdm.vegevol(Veg, Topo, dhdt);

	return Veg,Topo
