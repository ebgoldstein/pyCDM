"""
CDM Wind

rewrite of CDM by EBG 2/16

CDM has 5 types of wind:
	constant
	real
	flatrand
	sine
	bidirectional

	The most critical right now is constant and real.

"""

from CDMparams import *


#constant wind
m_ustar = constwind;
m_dir = 0; # onshore wind

#real wind file
