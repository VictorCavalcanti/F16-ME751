# F16-ME751
SimEngine3D is a modest 3D-Kinematics and Dynamics simulator.
Users must provide both "project_wrapper.mdl" and "bodies_info.adm" 
and could also provide a "project.acf" file.

"engine_wrapper.mdl" - Provides the bulk of the calls to simEngine3D.
Needed for any post processing.

"bodies_info.adm" - Provides information about the attributes
and constraints of the bodies involved in the simulation.

(optional but recommended) "project.adm" - Provides information 
about the type of simulation to be run, the duration of the 
simulation, step size in integrator, etc. A default setting will be used 
if no file is provided.
