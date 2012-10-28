#ifndef SCATTERING_TRANSITIONS_OBJECTS_HPP
#define SCATTERING_TRANSITIONS_OBJECTS_HPP

#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/value_at_key.hpp>
#include <cmath>

struct random_source_tag ;

namespace scattering{
namespace stateselection{


template<typename ParameterType, typename StateType>
void elastic_isotropic_transition(const ParameterType& parameters
                                  , StateType& state)
{
   typedef typename
      boost::fusion::result_of::at_key<ParameterType, random_source_tag>
      ::result_type random_value_type ;
   
   random_value_type random_1
      (boost::fusion::at_key<random_source_tag>(parameters)())
      ;
   random_value_type random_2
      (boost::fusion::at_key<random_source_tag>(parameters)())
      ;
         
   random_value_type cos_phi( cos ( 2. * M_PI *  random_1 ) );
   random_value_type sin_phi( sin ( 2. * M_PI *  random_1 ) );
   
   random_value_type cos_theta( 1.0 - 2.0 * random_2 );
   random_value_type sin_theta( sqrt(1.0 - cos_theta * cos_theta) );
  

}

}

template<typename ParameterType>
struct elastic_isotropic_transition
{
   ParameterType& parameters;

   elastic_isotropic_transition(ParameterType& parameters) :
      parameters(parameters) {}

   typedef double return_type;

   typedef typename
   boost::fusion::result_of::value_at_key<ParameterType
                                          , random_source_tag>::type
   ::result_type random_value_type ;
   
   template<typename StateType>
   inline StateType operator()(const StateType& state) const
   {
      using namespace utility ;

      random_value_type random_1 //(0.3)
         (boost::fusion::at_key<random_source_tag>(parameters)())
         ;
      random_value_type random_2 //(0.6)
         (boost::fusion::at_key<random_source_tag>(parameters)())
         ;

      random_value_type cos_phi( cos ( 2. * M_PI *  random_1 ) );
      random_value_type sin_phi( sin ( 2. * M_PI *  random_1 ) );
      
      random_value_type cos_theta( 1.0 - 2.0 * random_2 );
      random_value_type sin_theta( sqrt(1.0 - cos_theta * cos_theta) );
      
      // fixme [PS]
      std::valarray<double> inner_momentum(3);


      inner_momentum[0] =
         boost::units::quantity_cast<double>(
            sqrt(momentum_square_modulus(query_momentum(state))) 
            * sin_theta * cos_phi
         ) ;
      inner_momentum[1] = boost::units::quantity_cast<double>(
         sqrt(momentum_square_modulus(query_momentum(state))) 
         * sin_theta * sin_phi
      ) ;
      inner_momentum[2] = boost::units::quantity_cast<double>(
         sqrt(momentum_square_modulus(query_momentum(state)))
         * cos_theta
      ) ;

      typedef typename query_position_type<StateType>::type position_type;
      typedef typename query_momentum_type<StateType>::type momentum_type;

      return StateType(query_position(state)
                       , momentum_type(inner_momentum 
                                       / boost::units::si::meters)
                       , std::get<2>(state)
      );
   }
};



template<typename ParameterType> // , typename EnergyType>
struct inelastic_isotropic_transition
{
   ParameterType& parameters;
   typedef boost::units::quantity<boost::units::si::energy> EnergyType;
   const EnergyType delta_energy;

   inelastic_isotropic_transition(ParameterType& parameters
                                  , EnergyType delta_energy) :
      parameters(parameters), delta_energy(delta_energy)
      ,hbar(boost::fusion::at_key<constants::hbar>(parameters))
      ,mass(boost::fusion::at_key<parameters::mass>(parameters))
   {}

   typedef double return_type;

   typedef typename
   boost::fusion::result_of::value_at_key<ParameterType
                                          , random_source_tag>::type
   ::result_type random_value_type ;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::hbar>::type hbar ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::mass>::type mass ;
   
   template<typename StateType>
   inline StateType operator()(const StateType& state) const
   {
      using namespace utility ;

      random_value_type random_1
         (boost::fusion::at_key<random_source_tag>(parameters)())
         ;
      random_value_type random_2
         (boost::fusion::at_key<random_source_tag>(parameters)())
         ;

      random_value_type cos_phi( cos ( 2. * M_PI *  random_1 ) );
      random_value_type sin_phi( sin ( 2. * M_PI *  random_1 ) );
      
      random_value_type cos_theta( 1.0 - 2.0 * random_2 );
      random_value_type sin_theta( sqrt(1.0 - cos_theta * cos_theta) );
      

      EnergyType old_gamma(energy_dispersion(parameters, state));

      EnergyType old_energy
         (
            2. * old_gamma
            / (1. 
               + sqrt(4.
                      * boost::fusion::at_key<scattering::parameters::alpha>
                      (parameters) 
                      * old_gamma
                      + 1. ) 
            )
         )
         ;
      EnergyType new_energy(old_energy + delta_energy);
      EnergyType new_gamma(gamma_dispersion(parameters, new_energy));

      // fixme [PS]
      std::valarray<double> inner_momentum(3);

      typedef boost::units::quantity<boost::units::si::wavenumber>
         momentum_length_type;
      momentum_length_type momentum_length(
         sqrt(new_gamma * 2.0 * 0.067 
              * boost::units::si::constants::codata::m_e 
              / boost::units::si::constants::codata::hbar 
              / boost::units::si::constants::codata::hbar)
      );

      inner_momentum[0] = 
         boost::units::quantity_cast<double>(
            momentum_length * sin_theta * cos_phi
         ) ;
      inner_momentum[1] = boost::units::quantity_cast<double>(
         momentum_length * sin_theta * sin_phi
      ) ;
      inner_momentum[2] = boost::units::quantity_cast<double>(
         momentum_length * cos_theta
      ) ;

      typedef typename query_position_type<StateType>::type position_type;
      typedef typename query_momentum_type<StateType>::type momentum_type;

      return StateType(query_position(state)
                       , momentum_type(inner_momentum 
                                       / boost::units::si::meters)
                       , 1.0
      );
   }
};



template<typename ParameterType>
struct polar_optic_transition
{
   ParameterType& parameters;
   typedef boost::units::quantity<boost::units::si::energy> EnergyType;
   const EnergyType delta_energy;

   polar_optic_transition(ParameterType& parameters, EnergyType delta_energy) :
      parameters(parameters), delta_energy(delta_energy)
      ,hbar(boost::fusion::at_key<constants::hbar>(parameters))
      ,mass(boost::fusion::at_key<parameters::mass>(parameters))
   {}

   typedef double return_type;

   typedef typename
   boost::fusion::result_of::value_at_key<
      ParameterType
      , random_source_tag>::type::result_type random_value_type ;

   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , constants::hbar>::type hbar ;
   typename boost::fusion::result_of::value_at_key<
      ParameterType
      , parameters::mass>::type mass ;
   
   template<typename StateType>
   inline StateType operator()(const StateType& state) const
   {
      using namespace utility ;

      random_value_type random_1
         (boost::fusion::at_key<random_source_tag>(parameters)())
         ;
      random_value_type random_2
         // (0.6)
         (boost::fusion::at_key<random_source_tag>(parameters)())
         ;


      typedef typename query_momentum_type<StateType>::type momentum_type;

      momentum_type old_momentum(query_momentum(state));
      typedef boost::units::quantity<boost::units::si::wavenumber> 
         momentum_length_type;
      momentum_length_type old_momentum_length 
         (sqrt(momentum_square_modulus(old_momentum)));
      
      boost::units::quantity<boost::units::si::wavenumber> 
         momentum_modulus(sqrt(
                             momentum_square_modulus(old_momentum)
                          )
         );

      random_value_type cos_alpha =
         boost::units::quantity_cast<std::valarray<double> > 
         (old_momentum)[0] 
         / boost::units::quantity_cast<double>(momentum_modulus) ;
      random_value_type sin_alpha = sqrt( 
         boost::units::quantity_cast<std::valarray<double>>(old_momentum)[1] 
         * boost::units::quantity_cast<std::valarray<double>>(old_momentum)[1] 
         +
         boost::units::quantity_cast<std::valarray<double>>(old_momentum)[2] 
         * boost::units::quantity_cast<std::valarray<double>>(old_momentum)[2]
      ) / boost::units::quantity_cast<double>(momentum_modulus) ;


      random_value_type cos_beta;
      random_value_type sin_beta;

      if (sin_alpha > 1.e-10)
      {
         cos_beta = 
            boost::units::quantity_cast<random_value_type>
            (boost::units::quantity_cast<std::valarray<double> >
             (old_momentum)[1] / momentum_modulus / sin_alpha);
         sin_beta = 
            boost::units::quantity_cast<random_value_type>
            (boost::units::quantity_cast<std::valarray<double>>
             (old_momentum)[2] / momentum_modulus / sin_alpha
            );
      }
      else 
      {
         cos_beta = 1.0;
         sin_beta = 0.;
      } 

      random_value_type cos_phi( cos ( 2. * M_PI *  random_1 ) );
      random_value_type sin_phi( sin ( 2. * M_PI *  random_1 ) );

      EnergyType old_gamma(energy_dispersion(parameters, state));

      EnergyType old_energy
         (
            2. * old_gamma 
            / (1. 
               + sqrt(
                  4.
                  * boost::fusion::at_key<scattering::parameters::alpha>
                  (parameters) * old_gamma 
                  + 1. ) 
            )
         )
         ;
      EnergyType new_energy(old_energy + delta_energy);
      EnergyType new_gamma(gamma_dispersion(parameters, new_energy));

       random_value_type sum_gamma(
          boost::units::quantity_cast<random_value_type>(old_gamma + new_gamma)
       );
      random_value_type sum_sqrt_gamma(
         boost::units::quantity_cast<random_value_type>(
            sqrt(old_gamma) * sqrt(new_gamma)
         )
      );

      random_value_type cos_theta (
         (
            sum_gamma - pow((sum_gamma + 2.*sum_sqrt_gamma), random_2)
            *pow((sum_gamma - 2.*sum_sqrt_gamma), (1.-random_2))
         ) 
         / (2.*sum_sqrt_gamma));

      random_value_type sin_theta( sqrt(1.0 - cos_theta * cos_theta) );

      momentum_length_type new_momentum_length(
         sqrt(new_gamma * 2.0 * 0.067 
              * boost::units::si::constants::codata::m_e
              / boost::units::si::constants::codata::hbar
              / boost::units::si::constants::codata::hbar)
      );
      
      // fixme [PS]
      std::valarray<double> inner_momentum(3);

      inner_momentum[0] =
         boost::units::quantity_cast<double>(
            new_momentum_length 
            * (cos_theta * cos_alpha - sin_theta * cos_phi * sin_alpha)
         );
      inner_momentum[1] =
         boost::units::quantity_cast<double>(
            new_momentum_length 
            * ((cos_theta * sin_alpha + sin_theta * cos_phi * cos_alpha) 
               * cos_beta - sin_theta * sin_phi * sin_beta)
         );
      inner_momentum[2] =
         boost::units::quantity_cast<double>(
            new_momentum_length 
            * ((cos_theta * sin_alpha + sin_theta * cos_phi * cos_alpha)
               * sin_beta + sin_theta * sin_phi * cos_beta)
         );
      
      typedef typename query_position_type<StateType>::type position_type;
      typedef typename query_momentum_type<StateType>::type momentum_type;

      return StateType(std::get<0>(state)
                       , momentum_type(inner_momentum 
                                       / boost::units::si::meters)
                       , std::get<2>(state)
      );
   }
};

}




#endif // SCATTERING_TRANSITIONS_OBJECTS_HPP
