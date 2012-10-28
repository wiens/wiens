#ifndef ENERGY_UTILITY_HPP
#define ENERGY_UTILITY_HPP

#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/value_at_key.hpp>

#include "scattering_tags.hpp"
#include "momentum_utility.hpp"

namespace utility
{

template<typename StateType, typename ParameterType>
typename energy_dispersion_type<ParameterType, StateType>::type
energy_dispersion(const ParameterType& parameters, const StateType& state)
{
   typename boost::fusion::result_of::value_at_key<
   ParameterType
      , scattering::constants::hbar>::type hbar
      (boost::fusion::at_key<scattering::constants::hbar>(parameters))
      ;
   typename boost::fusion::result_of::value_at_key<
   ParameterType
      , scattering::parameters::mass>::type mass 
      (boost::fusion::at_key<scattering::parameters::mass>(parameters))
      ;
   return hbar *  hbar / 2. / mass 
      * utility::momentum_square_modulus(tuple_state::query_momentum(state))
   ;
}

template<typename ParameterType, typename EnergyType>
EnergyType
gamma_dispersion(const ParameterType& parameters, const EnergyType& energy)
{
   typename boost::fusion::result_of::value_at_key<
   ParameterType
      , scattering::parameters::alpha>::type alpha 
      (boost::fusion::at_key<scattering::parameters::alpha>(parameters))
      ;
   return energy * ( 1. + alpha * energy );
}


template<typename ParameterType, typename EnergyType>
inline EnergyType
energy_from_gamma(const ParameterType& parameters, const EnergyType& gamma)
{
   return (
      2. * gamma 
      / (1. 
         + sqrt(4.
                * boost::fusion::at_key<scattering::parameters::alpha>
                (parameters) * gamma + 1. ) 
      )
   ) ;
}

template<typename ParameterType, typename StateType>
typename energy_dispersion_type<ParameterType, StateType>::type
gamma_dispersion2(const ParameterType& parameters, const StateType& state)
{
   typename boost::fusion::result_of::value_at_key<
   ParameterType
      , scattering::constants::hbar>::type hbar
      (boost::fusion::at_key<scattering::constants::hbar>(parameters))
      ;
   typename boost::fusion::result_of::value_at_key<
   ParameterType
      , scattering::parameters::mass>::type mass 
      (boost::fusion::at_key<scattering::parameters::mass>(parameters))
      ;
   return hbar *  hbar / 2. / mass 
      * utility::momentum_square_modulus(query_momentum(state))
   ;

}

}

#endif // ENERGY_UTILITY_HPP
