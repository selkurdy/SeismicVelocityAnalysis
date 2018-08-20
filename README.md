# SeismicVelocityAnalysis
Fit quadratic functions to seismic velocities for depth conversion with probabilities  
**svelansegy.py** reads either a flat file of x,y,t[z],vrms[vi,vav] or segy file of  
stacking or migration velocities. It generates time-depth pairs and computes quadratic  
coefficients  and its range of probabilities.  

The resulting file can be used to depth or time convert any horizon. The approach spares  
the user from having to build velocity models. The whole process is based on generating just  
2 or 3 grids of quadratic coefficients.
