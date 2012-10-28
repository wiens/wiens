This is a small tool to calculate the flow through phase space as 
encountered by an electron in a semiconductor. The semiconductor is
modelled using a single non parabolic band. Several scattering mechanisms
are included, which represent the electron's interaction with the crystal.
The simulation setup allows for the simulation of the free evolution of 
an initial condition for several time steps under the influence of 
scattering. In the case a quantum mechanical initial condition is 
supplied, it is possible to observe the process of decoherence due to the 
interactions with the environment.

Compilation is handled by the waf build tool. You can find detailed
documentation on its web site: http://code.google.com/p/waf/.

To just get you started you first need to configure the build
environment by running:
   ./waf configure
This should locate all the required libraries and tools. In particular
a compiler capable of C++0x is required along with the Boost libraries
(www.boost.org) and the Lua interpreted language (www.lua.org).
A subsequent call of:
   ./waf
Will initiate the build and should result in the executable in 'built'
subdirectory. After the build has finished, the example can be started 
  ./built/phase_space_evolution_demo config.lua parameters.lua

