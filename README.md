# SeismicVelocityAnalysis
Fit quadratic functions to seismic velocities for depth conversion with probabilities  
**svelansegy.py** reads either a flat file of x,y,t[z],vrms[vi,vav] or segy file of  
stacking or migration velocities. It generates time-depth pairs and computes quadratic  
coefficients  and its range of probabilities.  

The resulting file can be used to depth or time convert any horizon. The approach spares  
the user from having to build velocity models. The whole process is based on generating just  
2 or 3 grids of quadratic coefficients.

**quadreg.py** generates time, depth, and average velocity slices to QC quadratic functions.  
It can also read in a flat x y z file for depth/time conversion. It needs the quadratic coefs  
file generated from **svelansegy.py**. The output is a text file that can be imported into Petrel  
or any mapping system.
