#ifndef PARAMETER_SETUP_PHYSICAL_HPP
#define PARAMETER_SETUP_PHYSICAL_HPP

#include "scattering_tags.hpp"

#include <boost/units/systems/si.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/cmath.hpp>

#include <boost/random.hpp>

#include <boost/fusion/adapted/struct/adapt_assoc_struct.hpp>

namespace boost {
   namespace units {
      typedef derived_dimension< length_base_dimension, -3 >
      ::type density_dimension;  // derived dimension for density : l^-3 
      typedef derived_dimension<length_base_dimension
                                , -2
                                , mass_base_dimension
                                , -1
                                , time_base_dimension
                                , 2 >
      ::type alpha_dimension;  // derived dimension for alpha : J^-1 
  }
}

namespace boost {
  namespace units {
    namespace si {
      typedef unit< density_dimension, si::system > density;
      typedef unit< alpha_dimension, si::system > alpha;
    }
  }
}

typedef boost::units::divide_typeof_helper<boost::units::si::energy
                                           ,boost::units::si::length>
::type energy_over_length;

namespace constantss = boost::units::si::constants::codata;

struct random_generator_wrapper
{
   typedef double numeric_type;
   typedef boost::mt19937 engine_type;
   typedef boost::uniform_real<numeric_type> distribution_type;
   typedef boost::variate_generator<engine_type&, distribution_type > 
   generator_type;
   
   engine_type engine ;
   distribution_type distribution ;
   generator_type generator ;
   
   random_generator_wrapper() : 
      engine(), distribution(), generator(engine, distribution)
   {
      std::cerr << "constructed random generator" << std::endl ;
   }
   typedef numeric_type result_type ;

   numeric_type operator()()
   {
      numeric_type value ;
      {
         value = (generator()) ;
         
      }
      return value ;
   }
};

struct params_physical
{
   boost::units::quantity<boost::units::si::density> NI;
   boost::units::quantity<boost::units::si::temperature> TL;
   double Z;
   boost::units::quantity<boost::units::si::wavenumber> beta;
   boost::units::quantity<boost::units::si::mass> effmass;
   boost::units::quantity<boost::units::si::alpha> alpha;
   boost::units::quantity<boost::units::si::length> distance;
   boost::units::quantity<boost::units::si::energy> energy_beta;
   boost::units::quantity<boost::units::si::electric_charge> e;
   boost::units::quantity<constantss::energy_time> hbar;
   boost::units::quantity<constantss::capacitance_over_length> epsilon_0;
   boost::units::quantity<constantss::energy_over_temperature> k_B;

   boost::units::quantity<boost::units::si::mass_density> rho;
   boost::units::quantity<boost::units::si::velocity> us;
   boost::units::quantity<boost::units::si::energy> DA;

   boost::units::quantity<energy_over_length> DO;
   boost::units::quantity<boost::units::si::energy> optical_phonon_energy;
   double NO;

   boost::units::quantity<constantss::capacitance_over_length> 
   optical_permitivity ;
   boost::units::quantity<constantss::capacitance_over_length>
   static_permitivity ;

   double NLO ;
   boost::units::quantity<boost::units::si::energy>
   polar_optical_phonon_energy ;
   boost::units::quantity<boost::units::si::time> free_flight_coeff;
   
   params_physical() : generator() 
   {
      NI= 1 / boost::units::si::cubic_meters ;
      TL = 300 * boost::units::si::kelvin ;
      Z=0.5e-12;

      effmass = constantss::m_e * 0.067;
      alpha = 0.61 / constantss::e / boost::units::si::volt;
      distance = 5.642e-10 * boost::units::si::meters;

      e = constantss::e;
      hbar = constantss::hbar;
      epsilon_0 = constantss::epsilon_0;
      k_B = constantss::k_B;

      beta =       
         sqrt(
            constantss::e * constantss::e * NI / 
            constantss::epsilon_0 / constantss::k_B / TL
         )
         ;
      energy_beta = 
         constantss::hbar * constantss::hbar * beta * beta / 2. / effmass
         ;
      rho = 5.36e3 * boost::units::si::kilogram_per_cubic_meter;
      us = 1./3. 
         * ( 2. * 3.0e3  + 5.24e3 ) 
         * boost::units::si::meters / boost::units::si::seconds ;
      DA = 7.0 * constantss::e * boost::units::si::volt ;

      beta =       
         sqrt(
            constantss::e * constantss::e * NI / 
            constantss::epsilon_0 / constantss::k_B / TL
         )
         ;
      energy_beta = 1 * boost::units::si::joule;

      DO = 1 * boost::units::si::joule / boost::units::si::meters ;
      optical_phonon_energy = 0.0343 * constantss::e * boost::units::si::volt ;
      NO =  1. / (exp(optical_phonon_energy / constantss::k_B / TL ) - 1.);
   
      polar_optical_phonon_energy = 
         0.036 * constantss::e * boost::units::si::volt ;
      optical_permitivity = constantss::epsilon_0 * 10.92;
      static_permitivity = constantss::epsilon_0 * 12.9;
      NLO =  1. 
         / (exp(polar_optical_phonon_energy / constantss::k_B / TL ) - 1.);

      free_flight_coeff = 1. / (2.6e13 * boost::units::si::hertz) ;
      alpha = 0.61 
         / constantss::e 
         / boost::units::si::volt ;
   }

   template<typename T>
   void reset_temperature(T new_temperature)
   {
      TL = new_temperature;
      NO =  1. / (exp(optical_phonon_energy / constantss::k_B / TL ) - 1.);
      NLO =  1. 
         / (exp(polar_optical_phonon_energy / constantss::k_B / TL ) - 1.);
      beta =       
         sqrt(
            constantss::e * constantss::e * NI / 
            constantss::epsilon_0 / constantss::k_B / TL
         )
         ;
   }

   typedef boost::units::quantity<boost::units::si::frequency> rate_type;

   typedef random_generator_wrapper generator_type ;
   generator_type generator ; 
   
};


BOOST_FUSION_ADAPT_ASSOC_STRUCT(
   params_physical,
   (boost::units::quantity<constantss::energy_time>
    , hbar, scattering::constants::hbar)
   (boost::units::quantity<constantss::energy_over_temperature>
    , k_B, scattering::constants::k_B)
   (boost::units::quantity<boost::units::si::electric_charge>
    , e, scattering::constants::e)
   (boost::units::quantity<boost::units::si::temperature>
    , TL, scattering::parameters::TL)
   (boost::units::quantity<boost::units::si::mass>
    , effmass, scattering::parameters::mass)
   (boost::units::quantity<boost::units::si::alpha>
    , alpha, scattering::parameters::alpha)
   (boost::units::quantity<boost::units::si::energy>
    , energy_beta, scattering::parameters::detail::energy_beta)
   (double, Z, scattering::parameters::detail::Z)
   (boost::units::quantity<boost::units::si::density>
    , NI, scattering::parameters::detail::NI)
   (boost::units::quantity<boost::units::si::mass_density>
    , rho, scattering::parameters::rho)
   (boost::units::quantity<boost::units::si::velocity>
    , us, scattering::parameters::detail::us)
   (boost::units::quantity<boost::units::si::energy>
    , DA, scattering::parameters::detail::DA)
   (boost::units::quantity<energy_over_length>
    , DO, scattering::parameters::detail::DO)
   (double, NO, scattering::parameters::detail::NO)
   (boost::units::quantity<boost::units::si::energy>
    , optical_phonon_energy
    , scattering::parameters::detail::optical_phonon_energy)
   (boost::units::quantity<constantss::capacitance_over_length>
    , optical_permitivity, scattering::parameters::detail::optical_permitivity)
   (boost::units::quantity<constantss::capacitance_over_length>
    , static_permitivity, scattering::parameters::detail::static_permitivity)
   (double, NLO, scattering::parameters::detail::NLO)
   (boost::units::quantity<boost::units::si::energy>
    , polar_optical_phonon_energy
    , scattering::parameters::detail::polar_optical_phonon_energy)
   (boost::units::quantity<boost::units::si::time>
    , free_flight_coeff, flow::free_flight_coeff_tag)
   (params_physical::generator_type, generator, random_source_tag)
)


namespace utility
{
   template<typename T>
   struct energy_dispersion_type<params_physical, T>
   {
      typedef boost::units::quantity<boost::units::si::energy> type;
   };
   
   template<typename T>
   struct energy_dispersion_type<params_physical&, T>
   {
      typedef boost::units::quantity<boost::units::si::energy> type;
   };
}

#endif //  PARAMETER_SETUP_PHYSICAL_HPP
