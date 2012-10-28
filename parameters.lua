-- temperature of the crystal lattice
TL = 77

-- effective mass of the simulated particle (electron)
effmass = 0.067

-- non parabolicity factor for the energy dispersion
alpha = 0.61 

-- lattice constant
distance = 5.642e-10


-- ionized impurity scattering parameters

-- impurity density
NI= 1
-- ionization rate
Z = 0.5e-12
-- screening energy
energy_beta = 1


-- acoustic phonon intravalley scattering parameters:

-- mass density of the material
rho = 5.36e3
-- sound velocity in the material
us = 1./3. * ( 2. * 3.0e3  + 5.24e3 )
-- accoustic deformation potential
DA = 7.0


-- polar optical and optical phonon intravalley scattering parameters:

-- optical deformation potential
DO = 1
-- optical phonon energy
optical_phonon_energy = 0.0343
-- polar optical phonon energy
polar_optical_phonon_energy = 0.036
-- dielectric constant (hf)
optical_permitivity = 10.92;
-- dielectric constant (static)
static_permitivity = 12.9;


-- free flight coefficient
free_flight_coeff = 1. / 2.6e13 ;
