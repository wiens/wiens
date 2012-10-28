#ifndef SCATTERING_OBJECTS_HPP
#define SCATTERING_OBJECTS_HPP

#include <cmath>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/value_at_key.hpp>

#include "energy_utility.hpp"

namespace scattering{

   using namespace utility ;

template<typename ParameterType>
struct ion_rate
{
   typedef typename ParameterType::rate_type return_type;
   typedef typename ParameterType::rate_type rate_type;
   const ParameterType& parameters;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::hbar>::type hbar;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::k_B>::type k_B;
      
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::TL>::type TL;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::mass>::type mass;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::alpha>::type alpha ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::NI>::type NI ; 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::Z>::type Z ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::energy_beta>::type energy_beta ;
      
   ion_rate(const ParameterType& parameters) : 
      parameters(parameters)
      , hbar(boost::fusion::at_key<constants::hbar>(parameters))
      , k_B(boost::fusion::at_key<constants::k_B>(parameters))
      , TL(boost::fusion::at_key<parameters::TL>(parameters))
      , mass(boost::fusion::at_key<parameters::mass>(parameters))
      , alpha(boost::fusion::at_key<parameters::alpha>(parameters))
      , NI(boost::fusion::at_key<parameters::detail::NI>(parameters))
      , Z(boost::fusion::at_key<parameters::detail::Z>(parameters))
      , energy_beta(boost::fusion::at_key<parameters::detail::energy_beta>
                    (parameters))
   {}

   template<typename StateType>
   inline rate_type operator()(const StateType& state) const
   {
      typename boost::fusion::result_of::value_at_key<
      ParameterType
         , parameters::detail::energy_beta>::type energy
         (energy_dispersion(parameters, state))
         ;

      return 
         2. * k_B * k_B
         / hbar / hbar
         / hbar / hbar / M_PI / M_SQRT2
         * Z * Z
         * TL * TL
         * mass * mass
         / NI
         / sqrt(mass)
         * sqrt(gamma_dispersion(parameters, energy))
         * ( 1. + 2. * alpha * energy) 
         / ( 1. + 4. * energy / energy_beta )
      ;
   }

};


template<typename ParameterType>
struct acoustic_phonon_rate
{
   typedef typename ParameterType::rate_type return_type;
   typedef typename ParameterType::rate_type rate_type;

   const ParameterType& parameters;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::hbar>::type hbar ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::k_B>::type k_B ;
 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::TL>::type TL ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::mass>::type mass ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::alpha>::type alpha ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::rho>::type rho ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::DA>::type DA ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::us>::type us ;


   acoustic_phonon_rate(const ParameterType& parameters) :
      parameters(parameters)
      , hbar(boost::fusion::at_key<constants::hbar>(parameters))
      , k_B(boost::fusion::at_key<constants::k_B>(parameters))
      , TL(boost::fusion::at_key<parameters::TL>(parameters))
      , mass(boost::fusion::at_key<parameters::mass>(parameters))
      , alpha(boost::fusion::at_key<parameters::alpha>(parameters))
      , rho(boost::fusion::at_key<parameters::rho>(parameters))
      , DA(boost::fusion::at_key<parameters::detail::DA>(parameters))
      , us(boost::fusion::at_key<parameters::detail::us>(parameters))
   {}


   template<typename StateType>
   inline rate_type operator()(const StateType& state) const
   {
      typename boost::fusion::result_of::value_at_key<
         ParameterType
         , parameters::detail::energy_beta>::type gamma
         (gamma_dispersion2(parameters, state))
         ;
      typename boost::fusion::result_of::value_at_key<
         ParameterType
         , parameters::detail::energy_beta>::type energy
         (
            2. * gamma / (1. + sqrt(4.* alpha * gamma + 1. ) )
         )
         ;

      return 
         M_SQRT2 * k_B / hbar / hbar/ hbar/ hbar / M_PI 
         * DA * DA
         * TL  
         * mass * mass
         / rho
         / us / us
         / sqrt (mass) 
         * sqrt ( gamma )
         * ( 1. + 2. * alpha * energy )
         ;
   }

};

template<typename ParameterType>
struct optical_phonon_absorption_rate
{
   typedef typename ParameterType::rate_type return_type;
   typedef typename ParameterType::rate_type rate_type ;
   const ParameterType& parameters ;
   
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::hbar>::type hbar ;
 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::mass>::type mass ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::rho>::type rho ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::alpha>::type alpha ;
 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::DO>::type DO ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::NO>::type NO ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::optical_phonon_energy>
   ::type optical_phonon_energy ;


   optical_phonon_absorption_rate(const ParameterType& parameters) : 
      parameters(parameters)
      , hbar(boost::fusion::at_key<constants::hbar>(parameters))
      , mass(boost::fusion::at_key<parameters::mass>(parameters))
      , rho(boost::fusion::at_key<parameters::rho>(parameters))
      , alpha(boost::fusion::at_key<parameters::alpha>(parameters))
      , DO(boost::fusion::at_key<parameters::detail::DO>(parameters))
      , NO(boost::fusion::at_key<parameters::detail::NO>(parameters))
      , optical_phonon_energy(
         boost::fusion::at_key<
         parameters::detail::optical_phonon_energy
         >(parameters)
      )
   {}

   template<typename StateType>
   inline rate_type operator()(const StateType& state) const
   {
      typename boost::fusion::result_of::value_at_key<
      ParameterType
         , parameters::detail::energy_beta>::type gamma
         (gamma_dispersion2(parameters, state))
         ;
      typename boost::fusion::result_of::value_at_key<
         ParameterType
         , parameters::detail::energy_beta>::type old_energy
         (
            2. * gamma / (1. + sqrt(4.* alpha * gamma + 1. ) )
         )
         ;
      typename boost::fusion::result_of::value_at_key<
         ParameterType
         , parameters::detail::optical_phonon_energy>::type energy
         (old_energy + optical_phonon_energy)
         ;

      return
         M_SQRT1_2 / M_PI / hbar / hbar
         * DO * DO 
         / rho
         / optical_phonon_energy
         * mass * sqrt (mass) 
         * NO
         * sqrt ( gamma )
         * ( 1. + alpha * energy)
      ;
   }
};


template<typename ParameterType>
struct optical_phonon_emission_rate
{
   typedef typename ParameterType::rate_type return_type;
   typedef typename ParameterType::rate_type rate_type;
   const ParameterType& parameters;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::hbar>::type hbar ;
 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::mass>::type mass ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::rho>::type rho ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::alpha>::type alpha ;
 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::DO>::type DO ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::NO>::type NO ;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::optical_phonon_energy>::type 
   optical_phonon_energy ;

   optical_phonon_emission_rate(const ParameterType& parameters) : 
      parameters(parameters)
      , hbar(boost::fusion::at_key<constants::hbar>(parameters))
      , mass(boost::fusion::at_key<parameters::mass>(parameters))
      , rho(boost::fusion::at_key<parameters::rho>(parameters))
      , alpha(boost::fusion::at_key<parameters::alpha>(parameters))
      , DO(boost::fusion::at_key<parameters::detail::DO>(parameters))
      , NO(boost::fusion::at_key<parameters::detail::NO>(parameters))
      , optical_phonon_energy(
         boost::fusion::at_key<
         parameters::detail::optical_phonon_energy>(parameters)
      )
   {}

   template<typename StateType>
   inline rate_type operator()(const StateType& state) const
   {
      typename boost::fusion::result_of::value_at_key<
         ParameterType
         , parameters::detail::energy_beta>::type gamma
         (gamma_dispersion2(parameters, state))
         ;
      typename boost::fusion::result_of::value_at_key<
         ParameterType
         , parameters::detail::energy_beta>::type old_energy
         (
            2. * gamma / (1. + sqrt(4.* alpha * gamma + 1. ) )
         )
         ;
      typename boost::fusion::result_of::value_at_key<
         ParameterType
         , parameters::detail::optical_phonon_energy>::type energy
         (old_energy + optical_phonon_energy)
         ;

      if (energy < optical_phonon_energy)
      {
         return 0;
      }
      else
      {
         energy -= optical_phonon_energy;
         return
            M_SQRT1_2 / M_PI / hbar / hbar
            * DO * DO 
            / rho
            / optical_phonon_energy
            * mass * sqrt (mass) 
            * (NO + 1.)
            * sqrt ( gamma )
            * ( 1. + alpha * energy)
            ;
      }
   }
};



template<typename ParameterType>
struct polar_optical_phonon_absorption_rate
{
   typedef typename ParameterType::rate_type return_type;
   typedef typename ParameterType::rate_type rate_type;
   const ParameterType& parameters;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::hbar>::type hbar ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::e>::type e ;
 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::mass>::type mass ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::alpha>::type alpha ;
 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::optical_permitivity>::type optical_permitivity ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::optical_permitivity>::type static_permitivity ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::NLO>::type NLO ;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::polar_optical_phonon_energy>::type 
   polar_optical_phonon_energy ;

   polar_optical_phonon_absorption_rate(const ParameterType& parameters) : 
      parameters(parameters)
      ,hbar(boost::fusion::at_key<constants::hbar>(parameters))
      ,e(boost::fusion::at_key<constants::e>(parameters))
      ,mass(boost::fusion::at_key<parameters::mass>(parameters))
      ,alpha(boost::fusion::at_key<parameters::alpha>(parameters))
      ,optical_permitivity(
         boost::fusion::at_key<
         parameters::detail::optical_permitivity>
         (parameters))
      ,static_permitivity(boost::fusion::at_key<
                          parameters::detail::static_permitivity>
                          (parameters))
      ,NLO(boost::fusion::at_key<parameters::detail::NLO>(parameters))
      ,polar_optical_phonon_energy(
         boost::fusion::at_key<parameters::detail::polar_optical_phonon_energy>
         (parameters))
   {}

   template<typename StateType>
   inline rate_type operator()(const StateType& state) const
   {
      typename boost::fusion::result_of::value_at_key<
      ParameterType
         , parameters::detail::energy_beta>::type gamma
         (gamma_dispersion2(parameters, state))
         ;
      typename boost::fusion::result_of::value_at_key<
      ParameterType
         , parameters::detail::energy_beta>::type energy
         (
            2. * gamma / (1. + sqrt(4.* alpha * gamma + 1. ) )
         )
         ;

      typename boost::fusion::result_of::value_at_key<
      ParameterType
         , parameters::detail::polar_optical_phonon_energy>::type final_energy
         (energy + polar_optical_phonon_energy)
         ;
      
      return 
         e * e 
         / hbar / hbar / 4. / M_SQRT2 / M_PI 
         * polar_optical_phonon_energy
         * (1. / optical_permitivity - 1. / static_permitivity)
         * NLO 
         * sqrt(mass)
         * ( 1 + 2. * alpha * final_energy )
         / sqrt(gamma_dispersion(parameters, energy))
         * log (abs( 
                   (sqrt(gamma_dispersion(parameters, energy))
                    + sqrt(gamma_dispersion(parameters, final_energy)) ) /
                   (sqrt(gamma_dispersion(parameters, energy))
                    - sqrt(gamma_dispersion(parameters, final_energy)) ) 
                ))
         ;
   }
};


template<typename ParameterType>
struct polar_optical_phonon_emission_rate
{
   typedef typename ParameterType::rate_type return_type;
   typedef typename ParameterType::rate_type rate_type;
   const ParameterType& parameters;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::hbar>::type hbar ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::e>::type e ;
 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::mass>::type mass ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::alpha>::type alpha ;
 
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::optical_permitivity>::type optical_permitivity ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::optical_permitivity>::type static_permitivity ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::NLO>::type NLO ;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::detail::polar_optical_phonon_energy>::type 
   polar_optical_phonon_energy ;

   polar_optical_phonon_emission_rate(const ParameterType& parameters) : 
      parameters(parameters)
      , hbar(boost::fusion::at_key<constants::hbar>(parameters))
      , e(boost::fusion::at_key<constants::e>(parameters))
      , mass(boost::fusion::at_key<parameters::mass>(parameters))
      , alpha(boost::fusion::at_key<parameters::alpha>(parameters))
      , optical_permitivity(boost::fusion::at_key<
                            parameters::detail::optical_permitivity>
                            (parameters))
      , static_permitivity(boost::fusion::at_key<
                           parameters::detail::static_permitivity>(parameters))
      , polar_optical_phonon_energy(
         boost::fusion::at_key<
         parameters::detail::polar_optical_phonon_energy>(parameters))
   {
      NLO = (boost::fusion::at_key<parameters::detail::NLO>(parameters));
   }
   
   template<typename StateType>
   inline rate_type operator()(const StateType& state) const
   {
      typename boost::fusion::result_of::value_at_key<
      ParameterType
         , parameters::detail::energy_beta>::type gamma
         (gamma_dispersion2(parameters, state))
         ;
      typename boost::fusion::result_of::value_at_key<
      ParameterType
         , parameters::detail::energy_beta>::type energy
         (
            2. * gamma / (1. + sqrt(4.* alpha * gamma + 1. ) )
         )
         ;
      typename boost::fusion::result_of::value_at_key<
      ParameterType
         , parameters::detail::polar_optical_phonon_energy>::type final_energy
         (energy - polar_optical_phonon_energy)
         ;

      if (energy < polar_optical_phonon_energy)
      {
         return 0;
      }

      return 
         e * e 
         / hbar / hbar / 4. / M_SQRT2 / M_PI 
         * polar_optical_phonon_energy
         * (1. / optical_permitivity - 1. / static_permitivity)
         * (NLO + 1.)
         * sqrt(mass)
         * ( 1 + 2. * alpha * final_energy )
         / sqrt(gamma_dispersion(parameters, energy))
         * log (abs( 
                   (sqrt(gamma_dispersion(parameters, energy))
                    + sqrt(gamma_dispersion(parameters, final_energy)) ) /
                   (sqrt(gamma_dispersion(parameters, energy))
                    - sqrt(gamma_dispersion(parameters, final_energy)) ) 
                ))
         ;
   }
};


template<typename StreamType, typename ParameterType>
StreamType& operator<<(StreamType& out, const ion_rate<ParameterType>& )
{
   out << "ion rate";
   return out;
}

template<typename StreamType, typename ParameterType>
StreamType& operator<<(StreamType& out
                       , const acoustic_phonon_rate<ParameterType>& )
{
   out << "acoustic rate";
   return out;
}

template<typename StreamType, typename ParameterType>
StreamType& operator<<(StreamType& out
                       , const optical_phonon_absorption_rate<ParameterType>& )
{
   out << "optical phonon absorption rate";
   return out;
}
template<typename StreamType, typename ParameterType>
StreamType& operator<<(StreamType& out
                       , const optical_phonon_emission_rate<ParameterType>& )
{
   out << "optical phonon emission rate";
   return out;
}

template<typename StreamType, typename ParameterType>
StreamType& operator<<(StreamType& out
                       , const 
                       polar_optical_phonon_absorption_rate<ParameterType>& )
{
   out << "polar optical phonon absorption rate";
   return out;
}

template<typename StreamType, typename ParameterType>
StreamType& operator<<(StreamType& out
                       , const 
                       polar_optical_phonon_emission_rate<ParameterType>& )
{
   out << "polar optical phonon emission rate";
   return out;
}


}
#endif // SCATTERING_OBJECTS_HPP
