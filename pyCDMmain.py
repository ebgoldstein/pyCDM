#Main pyCDM script. runs the model
import CDMparams

#Make/get the initial conditions (should make these 3D for time)
#sand surface
CDMparams.Topo=CDMparams.TopoDomain(shore_HMWL,shore_watertable,beach_angle,dx,nx,ny)
#vegetation surface
CDMparams.VegDomain(nx,ny)


#for initial time: timestep: final time
  #DuneEvolution
  #save each increment
