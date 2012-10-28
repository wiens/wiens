#include <iostream>
#include <vector>

#include <fstream>
#include <string>

#include <map>
#include <tuple>

#include <valarray>

#include <boost/units/quantity.hpp>
#include <boost/units/systems/si.hpp>

#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/scaled_base_unit.hpp>

#include <boost/type_traits/add_reference.hpp>

#include <boost/lexical_cast.hpp>

#include "evolution_time.hpp"

#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/fold.hpp>

#include <boost/random.hpp>

#include "tuple_state.hpp"

#include "momentum_utility.hpp"
#include "customize.hpp"
#include "energy_utility.hpp"

using namespace tuple_state ;

#include "scattering_transitions_objects.hpp"
#include "scattering_objects.hpp"

#include "parameter_setup_physical.hpp"
#include "scattering_model.hpp"

namespace pd = boost::units::si;
namespace constants = boost::units::si::constants::codata;

#include "my_lua.hpp"


template<typename ParameterType, typename ItemType, typename StateType>
StateType fill_with_equil_state(ParameterType& parameters
                                , const ItemType& item, const StateType& state)
{
   typedef typename
      boost::fusion::result_of::value_at_key<ParameterType
      , random_source_tag>::type
      ::result_type random_value_type ;

   random_value_type random_component (
      -2. * log((boost::fusion::at_key<random_source_tag>(parameters)()))
   ) ;

   random_value_type theta(
      2 * M_PI * (boost::fusion::at_key<random_source_tag>(parameters)())
   );

   std::valarray<double> inner_momentum(3);
   
   typedef boost::units::quantity<pd::wavenumber>  momentum_length_type;
   momentum_length_type momentum_length(
      sqrt(
         random_component 
         / (boost::fusion::at_key<scattering::constants::hbar>(parameters)) 
         / (boost::fusion::at_key<scattering::constants::hbar>(parameters))
         * (boost::fusion::at_key<scattering::parameters::mass>(parameters)) 
         * (boost::fusion::at_key<scattering::parameters::TL>(parameters))
         * (boost::fusion::at_key<scattering::constants::k_B>(parameters))
      )

   );
      
   inner_momentum[0] =
      boost::units::quantity_cast<double>(std::get<1>(item)) ;
   inner_momentum[1] =
      boost::units::quantity_cast<double>(momentum_length * cos(theta)) ;
   inner_momentum[2] =
      boost::units::quantity_cast<double>(momentum_length * sin(theta)) ;

   return StateType(std::get<0>(item)
                    , inner_momentum / pd::meters
                    , std::get<2>(item));
}

template<typename LengthType, typename MomentumType>
struct calc_wigner_packet
{
   LengthType width;
   MomentumType mean_momentum;
   LengthType mean_position;
   LengthType momentum_width;


   calc_wigner_packet(LengthType width
                      , MomentumType mean_momentum
                      , LengthType mean_position)
      : width(width)
      , mean_momentum(mean_momentum)
      , mean_position(mean_position)
      , momentum_width(width)
   {}

   typedef typename boost::units::multiply_typeof_helper<
      LengthType
      , MomentumType>::type volume_type;
   typedef typename boost::units::power_typeof_helper<
      volume_type
      , boost::units::static_rational<-1,1> >::type result_type;
      
   result_type operator()(LengthType x, MomentumType p)
   {
      return
         1. 
         * exp ( - (x - mean_position) 
                 * (x - mean_position) 
                 / 2. / width / width )  
         * exp ( - (p - mean_momentum) 
                 * (p - mean_momentum) 
                 * 2. * momentum_width 
                 * momentum_width )
         ;
   }
};

template<typename LengthType, typename MomentumType>
struct calc_entangled_wigner_packet
{
   LengthType width;
   MomentumType mean_momentum;
   LengthType mean_position;
   LengthType momentum_width;


   calc_entangled_wigner_packet(LengthType width
                                , MomentumType mean_momentum
                                , LengthType mean_position) 
      : width(width)
      , mean_momentum(mean_momentum)
      , mean_position(mean_position)
      , momentum_width(width)
   {}

   typedef typename boost::units::multiply_typeof_helper<
      LengthType
      , MomentumType>::type volume_type;
   typedef typename boost::units::power_typeof_helper<
      volume_type
      , boost::units::static_rational<-1,1> >::type result_type;
      
   result_type operator()(LengthType x, MomentumType p)
   {
      double constant_1 (0.5);
      double constant_2 (0.5);
      return
         constant_1 * constant_1 
         * exp ( - (x - mean_position)
                 * (x - mean_position)
                 / 2. / width / width 
         )
         * exp ( - (p * p * 2. 
                    * momentum_width 
                    * momentum_width )
         )
         + 
         constant_2 * constant_2 
          * exp ( - (x + mean_position)
                  * (x + mean_position)
                  / 2. / width / width 
          )
         * exp ( - (p * p * 2. 
                    * momentum_width
                    * momentum_width )
         )
         + 
         2. * constant_1 * constant_2  
         * exp ( - (x) * (x) / 2. / width / width )  
         * exp ( - (p * p * 2. * momentum_width * momentum_width ))
         * cos(2. * mean_position * p) 
         ;
   }
};

template<typename LengthType, typename MomentumType>
struct calc_wigner_double_packet
{
   LengthType width;
   LengthType mean_position;
   MomentumType mean_momentum;
   MomentumType momentum_width;
   
   calc_wigner_double_packet(LengthType width
                             , MomentumType mean_momentum
                             , LengthType mean_position) 
      : width(width)
      , mean_position(mean_position)
      , mean_momentum(mean_momentum)
      , momentum_width(2. / width)
   {
      std::cerr << width << std::endl;
      std::cerr << momentum_width << std::endl;
   }

   typedef typename boost::units::multiply_typeof_helper<
      LengthType
      , MomentumType>::type volume_type;
   typedef typename boost::units::power_typeof_helper<
      volume_type
      , boost::units::static_rational<-1,1> >::type result_type;

   result_type operator()(LengthType x, MomentumType p)
   {
      return
         1.0 / (1. + exp(-mean_position*mean_position/width/width))
         * exp ( - (p - mean_momentum) 
                 * (p - mean_momentum) 
                 / 2. / momentum_width / momentum_width 
         )
         * (
            exp ( - (x - mean_position)
                  * (x - mean_position)
                  / 2. / width / width )
            + exp ( - (x + mean_position)
                    * (x + mean_position)
                    / 2. / width / width )
            + 2. * exp ( - x * x / 2. / width / width )
            * cos(p * 2. * mean_position )
         ) ; 
         ;
   }
};

namespace boost {
namespace units {
namespace si{

typedef derived_dimension<length_base_dimension,1,
                          mass_base_dimension,1,
                          time_base_dimension,-3,
                          current_base_dimension,-1>
::type electric_field_dimension ;

typedef unit<electric_field_dimension,system> electric_field;

}}}

template<typename T>
struct query_electric_field_type
{
   typedef boost::units::quantity<pd::electric_field, T > type;
};


template<typename SpaceType
         , typename PT
         , typename StateType
         , typename ParameterType>
inline void evolve_trajectory_2(const SpaceType& space
                                , const PT& source
                                , StateType& state
                                , ParameterType parameter)
{
   typedef double raw_quantity_type;

   // this is just a bridge to make the rest work for now ...
   // dummies have been implemented, but need further evolution
   //
   typedef typename query_electric_field_type<double>
      ::type electric_field_type;
   electric_field_type electric_field(electric_field_type::from_value(0));

   // constants currently hardcoded, modularize when time permits
   //
   typedef boost::units::quantity<pd::wavenumber> momentum_type;
   typedef typename query_position_type<StateType>::type location_type;
   momentum_type momentum_delta (constants::e * electric_field * parameter 
                                 / constants::hbar ) ;

   momentum_type new_momentum ( 
      (boost::units::quantity_cast<std::valarray<double> >
       (query_momentum(state) )[0] / pd::meters)
      + momentum_delta
   ) ;
 
   boost::units::quantity<pd::mass, double> electron_mass 
      (boost::fusion::at_key<scattering::parameters::mass>(source)) 
      ;

   boost::units::quantity<pd::energy, double> 
      gamma(utility::energy_dispersion(source, state))
      ;

   boost::units::quantity<pd::energy, double> energy
      (
         2. * gamma 
         / (1. 
            + sqrt(4.
                   * boost::fusion::at_key<scattering::parameters::alpha>
                   (source) 
                   * gamma + 1. ) 
         )
      )
      ;

    
    momentum_type mean_momentum ( 
       (
          (boost::units::quantity_cast<std::valarray<double> >
           (query_momentum(state) )[0] / pd::meters)
          + new_momentum
       ) / (2. * pd::si_dimensionless) 
    ) ;
 

    location_type location_delta(
       mean_momentum 
       * constants::hbar
       / electron_mass
       * ( parameter 
           * (2. * energy 
              * boost::fusion::at_key<scattering::parameters::alpha>(source)
              + 1.0
           )
       )
    );
 

    location_type new_position(query_position(state) + location_delta) ;

    state = std::make_tuple(new_position
                            , std::get<1>(state)
                            , std::get<2>(state)
    );
}


template<typename StreamType>
StreamType& operator<<(StreamType& out, 
                       std::tuple<
                       boost::units::quantity<pd::length>
                       , boost::units::quantity<pd::wavenumber
                       , std::valarray<double> >
                       , double
                       >& tuple)
{
   out << std::get<0>(tuple)
       << " " 
       << boost::units::quantity_cast<std::valarray<double> >
      (std::get<1>(tuple)) 
       << " " 
       << std::get<2>(tuple);
   return out;
}

template<typename StreamType>
StreamType& operator<<(StreamType& out, 
                       std::tuple<
                       boost::units::quantity<pd::length>
                       , boost::units::quantity<pd::wavenumber>
                       , double
                       >& tuple)
{
   return out;
}


template<typename ParameterType>
struct scattering_model_container
{
   typedef ParameterType parameter_type;
   typedef boost::fusion::vector<
      scattering_model<
         scattering::acoustic_phonon_rate<parameter_type>
         , scattering::elastic_isotropic_transition<parameter_type> >
      , 
      scattering_model<
         scattering::polar_optical_phonon_absorption_rate<parameter_type>
         , scattering::polar_optic_transition<parameter_type> >
      ,
      scattering_model<
         scattering::polar_optical_phonon_emission_rate<parameter_type>
         , scattering::polar_optic_transition<parameter_type> >
      > type;
};


template<typename ParameterType>
struct scattering_model_account
{
   typedef ParameterType parameter_type;
   typedef boost::fusion::map<
      boost::fusion::pair<
         scattering_model<
            scattering::acoustic_phonon_rate<parameter_type>
            , scattering::elastic_isotropic_transition<parameter_type> >
         , size_t >
      ,
      boost::fusion::pair<
         scattering_model<
         scattering::optical_phonon_emission_rate<parameter_type>
         , scattering::inelastic_isotropic_transition<parameter_type> >
         , size_t >
      ,
      boost::fusion::pair<
         scattering_model<
         scattering::optical_phonon_absorption_rate<parameter_type>
         , scattering::inelastic_isotropic_transition<parameter_type> >
         , size_t >
      ,
      boost::fusion::pair<
         scattering_model<
         scattering::polar_optical_phonon_absorption_rate<parameter_type>
         , scattering::polar_optic_transition<parameter_type> >
         , size_t >
      ,
      boost::fusion::pair<
         scattering_model<
         scattering::polar_optical_phonon_emission_rate<parameter_type>
         , scattering::polar_optic_transition<parameter_type> >
         , size_t >
   > type;
};


template<typename StreamType, typename S, typename T>
StreamType& operator<<(StreamType& out, const scattering_model<S, T>& model)
{
   out << "model: " << model.rate;
   return out;
}

template<typename StateType, typename NumericType=double>
struct evaluator
{
   const StateType& state;

   evaluator(const StateType& state) : state(state) {}

   typedef NumericType result_type;

   template<typename T>
   NumericType
   operator()(const NumericType initial_value, const T& item) // const 
   {
      return initial_value + item(state);
   }
};


template<typename StateType, typename NumericType=double>
struct evaluation_decider
{
   StateType& state;
   const NumericType limit;
   NumericType& active_rate;
   bool first;

   evaluation_decider(StateType& state
                      , const NumericType
                      limit, NumericType& active_rate) 
      : state(state), limit(limit), active_rate(active_rate), first(true) {}

   typedef NumericType result_type;

   template<typename T>
   NumericType
   operator()(const NumericType initial_value, const T& item) // const 
   {
      NumericType current_value(item(state));
      NumericType new_value(initial_value + current_value);

      std::cout << "at: "
                << item
                << std::endl;
      std::cout << "initial: " 
                << initial_value
                << " current contribution: "
                << current_value
                << std::endl;
      std::cout << "particle energy: "
                << utility::energy_from_gamma(
                   item.transition.parameters
                   , utility::gamma_dispersion2(
                      item.transition.parameters
                      , (state)) 
                )
                << std::endl;
      std::cout << "phonon energy: "
                << boost::fusion::at_key<
      scattering::parameters::detail::polar_optical_phonon_energy>
         (item.transition.parameters) 
                << std::endl;
      

      std::cout << new_value
                << " > "
                << limit
                << " ? "
                << first
                << std::endl;
      if ((new_value > limit) && first)
      {
         first = false;
         active_rate = new_value;

         StateType aaa(state);
         std::cout << "using: " << item << std::endl;
         std::cout << "\t\t\tswitching from:  " << aaa << " " 
                   << utility::energy_from_gamma(
                      item.transition.parameters
                      , utility::gamma_dispersion2(
                         item.transition.parameters, (aaa)
                      )
                   )
                   << std::endl;

         state = item.apply(state) ;

         std::cout << "\t\t\tfinal state now: "
                   << state
                   << " " 
                   << utility::energy_from_gamma(
                      item.transition.parameters
                      , utility::gamma_dispersion2(
                         item.transition.parameters, (state)
                      )
                   )
                   << std::endl;
      }
      return new_value;
   }
};

template<typename StateType
         , typename AccountingType
         , typename NumericType=double>
struct accounting_evaluation_decider
{
   StateType& state;
   const NumericType limit;
   NumericType& active_rate;
   AccountingType& account;
   bool first;

   accounting_evaluation_decider(StateType& state
                                 , AccountingType& account
                                 , const NumericType limit
                                 , NumericType& active_rate) 
      : state(state)
      , limit(limit)
      , active_rate(active_rate)
      , account(account)
      , first(true) {}

   typedef NumericType result_type;

   template<typename T>
   inline NumericType
   operator()(const NumericType initial_value, const T& item) // const 
   {
      NumericType current_value(item(state));
      NumericType new_value(initial_value + current_value);

      if ((new_value > limit) && first)
      {
         first = false;
         active_rate = new_value;

         boost::fusion::at_key<T>(account) += 1;
         state = item.apply(state) ;
      }
      return new_value;
   }
};

template<typename StateType, typename NumericType>
struct decider
{
   NumericType limit;
   const StateType& state;
   StateType final_state;
   size_t counter;

   decider(const StateType& state, NumericType limit) 
      : limit(limit), state(state), final_state(state), counter(0) 
   {} 

   template<typename T, typename S>
   bool operator()(S const& previous_value, T const& item) // const
   {
      counter++;
      bool border_reached((previous_value) > limit);

      if (border_reached)
      {
         StateType aaa(state);
         std::cout << "using: " << item << std::endl;
         std::cout << "\t\t\tswitching from:  " 
                   << aaa
                   << " " 
                   << utility::gamma_dispersion2(
                      item.transition.parameters, (aaa)
                   )
                   << std::endl;
         final_state = item.apply(state);
         std::cout << "\t\t\tfinal state now: "
                   << final_state
                   << " " 
                   << utility::gamma_dispersion2(
                      item.transition.parameters, (final_state))
                   << std::endl;
      } else
      {
         std::cout << "not using: " << item << std::endl;
      }
      
      return border_reached;
   }
};


template<typename DataStoreType
         , typename ParameterType
         , typename TimeType
         , typename AccountingType
         , typename ScatteringModelType
         , typename BookType 
         , typename ConfigType 
         , typename FlightTimeType >
void work(DataStoreType& data_store
          , ParameterType& parameters
          , TimeType timestep
          , AccountingType& accounting
          , ScatteringModelType& scattering_models
          , BookType& book
          , const ConfigType& config
          , FlightTimeType& flight_time
          , size_t max_iterations
          , size_t max_particle_count)
{
   typedef DataStoreType storage_type ;
   typedef ParameterType parameter_type ;
   using boost::fusion::at_key ;
   using boost::units::quantity ;

   typedef typename storage_type::value_type data_type ;
   typedef typename query_momentum_type<data_type>::type momentum_type ;
   typedef typename query_position_type<data_type>::type position_type ;

   typedef std::valarray<double> inner_type;
   typedef quantity<pd::wavenumber, inner_type> vector_momentum_type;
   typedef std::tuple<position_type, vector_momentum_type, double> state_type;
   state_type empty;

   quantity<pd::frequency> max_rate;
   for (auto it(data_store.begin()); it != data_store.end(); ++it)
   {
      for (size_t current_particle(0);
           current_particle < max_particle_count ;
           ++current_particle)
      {
         state_type current(fill_with_equil_state(parameters, *it, empty));
         for (size_t ticker(1); ticker < max_iterations; ++ticker)
         {
            quantity<pd::time> current_flight_duration(flight_time(current));
            quantity<pd::time> remaining_flight_duration(timestep);

            while(remaining_flight_duration > quantity<pd::time>(0) )
            {
               if (current_flight_duration < remaining_flight_duration)
               {
                  evolve_trajectory_2(data_store
                                      , parameters
                                      , current
                                      , current_flight_duration);
                  const quantity<pd::frequency> current_rate(
                     at_key<random_source_tag>(parameters)() /
                     at_key<flow::free_flight_coeff_tag>(parameters)
                  );
                  quantity<pd::frequency> calculated_rate;
                  
                  evaluation_decider<state_type, typename 
                     parameter_type::rate_type> eval(current
                                                     , current_rate
                                                     , calculated_rate);
                  accounting_evaluation_decider<state_type, typename 
                     scattering_model_account<parameter_type>::type
                     , typename parameter_type::rate_type> 
                     account_eval(current
                                  , accounting
                                  , current_rate
                                  , calculated_rate
                     );

                  quantity<pd::frequency> total_rate(
                     boost::fusion::fold(
                        scattering_models
                        , typename parameter_type::rate_type(0)
                        , account_eval
                     )
                  );
                  
                  if (total_rate > max_rate)
                     max_rate = total_rate;

                  if (total_rate >
                      1. / at_key<flow::free_flight_coeff_tag>(parameters))
                  {
                     std::cerr << "max rate too low!!! "
                               << total_rate 
                               << " " 
                               << total_rate 
                        - 1. / at_key<flow::free_flight_coeff_tag>(parameters) 
                               << std::endl;
                  }

                  remaining_flight_duration -= current_flight_duration;
               }
               else
               {
                  evolve_trajectory_2(data_store
                                      , parameters
                                      , current
                                      , remaining_flight_duration) ;
                  remaining_flight_duration = 0;
               }
            }

            momentum_type imaged_momentum(
               (boost::units::quantity_cast<std::valarray<double> >
                (std::get<1>(current))[0]) 
               / pd::meters 
            );
            
            size_t location_index( 
               std::round ( 
                  (std::get<0>(current) + 0.5 * config.L)
                  / config.L * config.x_count 
               ) 
            );
            size_t momentum_index( 
               std::round ( 
                  (imaged_momentum + 0.5 * config.p_max)
                  / config.p_max * config.y_count) 
            );
            
            book.entropy_store
               [ticker]
               [std::tuple<size_t,size_t>(location_index, momentum_index)] +=
               (std::get<2>(current));

            if (( location_index < config.x_count ) && (location_index >= 0))
            {
               book.k_integrate[ticker][location_index] += 
                  (std::get<2>(current));
               book.x_outside[ticker]++;

               if ((std::get<2>(current)) < 0)
               {
                  book.k_integrate_negative[ticker][location_index] +=
                     (std::get<2>(current));
               }
               else
               {
                  book.k_integrate_positive[ticker][location_index] +=
                     (std::get<2>(current));
               }
            }

            if (( momentum_index < config.y_count ) && (momentum_index >= 0))
            {
               book.x_integrate[ticker][momentum_index] +=
                  (std::get<2>(current));
               book.k_outside[ticker]++;

               if ((std::get<2>(current)) < 0)
               {
                  book.x_integrate_negative[ticker][momentum_index] += 
                     (std::get<2>(current)) ;
               }
               else
               {
                  book.x_integrate_positive[ticker][momentum_index] += 
                     (std::get<2>(current)) ;
               }
            }
            else
            {
            }


            if (( location_index < config.x_count ) && (location_index >= 0))
               if (( momentum_index < config.y_count ) 
                   && (momentum_index >= 0))
               {
                  double old_value(fabs(
                                      book.timed_output_container[ticker]
                                      [location_index]
                                      [momentum_index]
                                   )
                  )
                  ;
                  book.timed_output_container[ticker]
                     [location_index]
                     [momentum_index]
                     += (std::get<2>(current));


                  if ((old_value != 0) 
                      && (old_value > fabs(book.timed_output_container[ticker]
                                           [location_index]
                                           [momentum_index]
                          )
                      )
                  )
                  {
                     book.timed_cancellation_map[ticker]
                        [location_index]
                        [momentum_index]
                        += fabs(std::get<2>(current));
                  }

                  book.total_outside[ticker]++;

                  book.entropy[ticker] += 
                     (std::get<2>(current)) * (std::get<2>(current)) ;
                  book.participating[ticker]++;
               }
         }
      }
   }
}



template<typename momentum_type, typename position_type>
struct config
{
   position_type L ;
   momentum_type p_max ; 
   size_t x_count ; // space 
   size_t y_count ; // momentum
   boost::units::quantity<pd::time> timestep ;
   size_t max_iterations ;
   size_t max_particle_count ;

   config() :
      L(800 * 1e-9*pd::meter)
      , p_max(6e9 / pd::meter) 
      , x_count(400)
      , y_count(3000)
      , timestep(100.0e-15 * pd::seconds)
      , max_iterations(10)
      , max_particle_count (30)
   {}
} ;

template<typename entropy_store_type
         , typename k_integrate_type
         , typename x_integrate_type
         , typename k_outside_type
         , typename x_outside_type
         , typename k_integrate_negative_type
         , typename k_integrate_positive_type
         , typename x_integrate_negative_type
         , typename x_integrate_positive_type
         , typename timed_output_container_type
         , typename timed_cancellation_map_type
         , typename total_outside_type
         , typename entropy_type
         , typename participating_type
         >
struct book
{
   book(entropy_store_type& entropy_store,
        k_integrate_type& k_integrate,
        x_integrate_type& x_integrate,
        k_outside_type& k_outside,
        x_outside_type& x_outside,
        k_integrate_negative_type& k_integrate_negative,
        k_integrate_positive_type& k_integrate_positive,
        x_integrate_negative_type& x_integrate_negative,
        x_integrate_positive_type& x_integrate_positive,
        timed_output_container_type& timed_output_container,
        timed_cancellation_map_type& timed_cancellation_map,
        total_outside_type& total_outside,
        entropy_type& entropy,
        participating_type& participating
   ) :
      entropy_store(entropy_store)
      ,k_integrate(k_integrate)
      ,x_integrate(x_integrate)
      ,k_outside(k_outside)
      ,x_outside(x_outside)
      ,k_integrate_negative(k_integrate_negative)
      ,k_integrate_positive(k_integrate_positive)
      ,x_integrate_negative(x_integrate_negative)
      ,x_integrate_positive(x_integrate_positive)
      ,timed_output_container(timed_output_container)
      ,timed_cancellation_map(timed_cancellation_map)
      ,total_outside(total_outside)
      ,entropy(entropy)
      ,participating(participating)
   {}

   entropy_store_type& entropy_store ;
   k_integrate_type& k_integrate ;
   x_integrate_type& x_integrate ;
   k_outside_type& k_outside ;
   x_outside_type& x_outside ;
   k_integrate_negative_type& k_integrate_negative ;
   k_integrate_positive_type& k_integrate_positive ;
   x_integrate_negative_type& x_integrate_negative ;
   x_integrate_positive_type& x_integrate_positive ;
   timed_output_container_type& timed_output_container ;
   timed_cancellation_map_type& timed_cancellation_map ;
   total_outside_type& total_outside ;
   entropy_type& entropy ;
   participating_type& participating ;
} ;

template<typename SourceType, typename ConfigType>
void set_config(SourceType& source, ConfigType& config)
{
   config.L = ( lua::extract<double>(source, "L")
                * boost::units::si::meters  ) ;
   config.p_max = ( lua::extract<double>(source, "p_max") 
                    / pd::meter ) ;
   config.x_count = ( lua::extract<long>(source, "x_count") ) ;
   config.y_count = ( lua::extract<long>(source, "y_count") ) ;
   config.timestep = ( lua::extract<double>(source, "timestep") 
                       * pd::seconds) ;
   config.max_iterations = ( lua::extract<long>(source, "max_iterations") ) ;
   config.max_particle_count =
      ( lua::extract<long>(source, "max_particle_count") ) ;
}

template<typename SourceType, typename ParametersType>
void set_parameters(SourceType& source, ParametersType& parameters)
{
   
   parameters.NI = lua::extract<double>(source, "NI") 
      / boost::units::si::cubic_meters ;
   parameters.TL = lua::extract<double>(source, "TL")
      * boost::units::si::kelvin ;
   parameters.Z = lua::extract<double>(source, "Z") ;

   parameters.effmass =
      constantss::m_e * lua::extract<double>(source, "effmass") ;
   parameters.alpha = lua::extract<double>(source, "alpha") 
      / constantss::e / boost::units::si::volt ; 
   parameters.distance = lua::extract<double>(source, "distance") 
      * boost::units::si::meters ;

   // recompute these
   //
   parameters.beta =       
      sqrt(
         constantss::e * constantss::e * parameters.NI / 
         constantss::epsilon_0 / constantss::k_B / parameters.TL
      )
      ;
   parameters.energy_beta = 
      constantss::hbar * constantss::hbar 
      * parameters.beta * parameters.beta / 2. / parameters.effmass
      ;
      
   parameters.rho = lua::extract<double>(source, "rho")
      * boost::units::si::kilogram_per_cubic_meter;
   parameters.us = lua::extract<double>(source, "us")
      * boost::units::si::meters / boost::units::si::seconds ;
   parameters.DA = lua::extract<double>(source, "DA")
      * constantss::e * boost::units::si::volt ;

   parameters.beta =       
      sqrt(
         constantss::e * constantss::e * parameters.NI / 
         constantss::epsilon_0 / constantss::k_B / parameters.TL
      )
      ;
      
   parameters.energy_beta = lua::extract<double>(source, "energy_beta")
      * boost::units::si::joule ;

   parameters.DO = lua::extract<double>(source, "DO")
      * boost::units::si::joule / boost::units::si::meters ;
   
   parameters.optical_phonon_energy = 
      lua::extract<double>(source, "optical_phonon_energy")
      * constantss::e * boost::units::si::volt ;

   parameters.NO =  1. / (exp(parameters.optical_phonon_energy 
                              / constantss::k_B / parameters.TL ) - 1.) ;
   
   parameters.polar_optical_phonon_energy =
      lua::extract<double>(source, "polar_optical_phonon_energy")
      * constantss::e * boost::units::si::volt ;

   parameters.optical_permitivity = constantss::epsilon_0 
      * lua::extract<double>(source, "optical_permitivity") ;
   parameters.static_permitivity = constantss::epsilon_0 
      * lua::extract<double>(source, "static_permitivity") ;
   parameters.NLO =  1. / (exp(parameters.polar_optical_phonon_energy 
                               / constantss::k_B / parameters.TL ) - 1.) ;

   parameters.free_flight_coeff = 
      (lua::extract<double>(source, "free_flight_coeff") 
       / boost::units::si::hertz) ;
}


int main(int argc, char** argv)
{
   if (argc < 3)
   {
      std::cerr << "error! insufficient parameters" << std::endl ;
      std::cerr << "please use: "
                << argv[0] 
                << " config_file parameter_file"
                << std::endl ;
      exit(-2) ;
   }

   using boost::fusion::at_key ;
   using boost::units::quantity ;

   typedef double numeric_type;
   typedef std::vector< numeric_type > inner_matrix_type;
   typedef std::vector< inner_matrix_type > matrix_type;
   
   typedef quantity<pd::length>   position_type;
   typedef quantity<pd::wavenumber> momentum_type;

   typedef calc_wigner_double_packet<position_type, momentum_type>
      calculator_type;

   position_type mean_position(0 * pd::meters);
   momentum_type mean_momentum(0 / pd::meters);

   typedef params_physical parameter_type;
   parameter_type parameters;

   lua::lua_value_wrapper parameters_source(argv[2]) ;
   set_parameters(parameters_source, parameters) ;

   position_type width(
      sqrt(
         at_key<scattering::constants::hbar>(parameters)
         * at_key<scattering::constants::hbar>(parameters) 
         / 4. / at_key<scattering::parameters::mass>(parameters)
         / at_key<scattering::constants::k_B>(parameters)
         / at_key<scattering::parameters::TL>(parameters)
      )
   );

   std::cout << width << std::endl;
   calc_wigner_packet<position_type, momentum_type> 
      calcme(width, mean_momentum, mean_position);
   
   config<momentum_type, position_type> cfg ;
   
   lua::lua_value_wrapper config_source(argv[1]) ;
   set_config(config_source, cfg) ;

   std::vector<double> entropy(cfg.max_iterations, 0);
   std::vector<size_t> participating(cfg.max_iterations, 0);
   typedef std::map<std::tuple<size_t, size_t>, double > entropy_store_type;
   std::vector<entropy_store_type> entropy_store(cfg.max_iterations);

   typedef std::tuple<position_type
      , momentum_type
      , calculator_type::result_type> data_type;
   typedef std::vector<data_type> storage_type;
   storage_type data_store;

   typedef std::map<long, std::map<long, double > > inner_record_type ; 
   std::map<long, inner_record_type> timed_output_container ;
   std::map<long, inner_record_type> timed_cancellation_map ;

   std::cerr 
      << "maximum extent: " 
      << "q: " <<  1. / cfg.x_count * cfg.L - 0.5 * cfg.L 
      << "\t" 
      << "p : " << 1. /cfg.y_count * cfg.p_max - cfg.p_max * 0.5
      << "\t"
      << 1. / cfg.x_count * cfg.L
      << "\t"
      << 1. / cfg.y_count * cfg.p_max
      << std::endl;

   {
      std::ofstream outputfile(std::string("initial"));

      for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
      {
         for (size_t count_y(0) ; count_y < cfg.y_count; 
              ++count_y)
         {

            if (abs
                (calcme(
                    count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L
                    , count_y * 1. /cfg.y_count * cfg.p_max - cfg.p_max * 0.5)
                ) > 1e-30)
               data_store.push_back(
                  data_type(count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L ,
                            count_y * 1. /cfg.y_count * cfg.p_max 
                            - cfg.p_max * 0.5
                            , (
                               calcme(
                                  count_x * 1. /cfg.x_count * cfg.L 
                                  - 0.5 * cfg.L
                                  , count_y * 1. /cfg.y_count * cfg.p_max
                                  - cfg.p_max * 0.5
                               )
                            )
                  )
               );
            timed_output_container[0][count_x][count_y] +=
               cfg.max_particle_count 
               * boost::units::quantity_cast<double>(
                  calcme(
                     count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L
                     , count_y * 1. /cfg.y_count * cfg.p_max - cfg.p_max * 0.5)
               );
            outputfile << boost::units::quantity_cast<double>(
               calcme(count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L
                      , count_y * 1. /cfg.y_count * cfg.p_max 
                      - cfg.p_max * 0.5
               )
            ) << " ";
            entropy[0] += cfg.max_particle_count 
               * boost::units::quantity_cast<double>(
                  calcme(
                     count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L
                     , count_y * 1. /cfg.y_count * cfg.p_max - cfg.p_max * 0.5)
               )
               * boost::units::quantity_cast<double>(
                  calcme(count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L
                         , count_y * 1. /cfg.y_count * cfg.p_max
                         - cfg.p_max * 0.5)
               ) ;

            entropy_store[0][std::tuple<size_t, size_t>(count_x,count_y)] = 
               cfg.max_particle_count
               * boost::units::quantity_cast<double>(
                  calcme(count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L
                         , count_y * 1. /cfg.y_count * cfg.p_max 
                         - cfg.p_max * 0.5)
               )
               ;
         }
         outputfile << std::endl;
      }
      outputfile.close();
   }
   std::cout << "active particles: " << data_store.size() << std::endl;
   participating[0] = data_store.size() * cfg.max_particle_count; 
   
   flow::evolution_time<parameter_type> flight_time(parameters);

   scattering::acoustic_phonon_rate<parameter_type> acoustic_rate(parameters);

   scattering::elastic_isotropic_transition<parameter_type> acoustic_transition
   (parameters);
   scattering::inelastic_isotropic_transition<parameter_type> 
   optical_absorption_transition(
      parameters
      , at_key<scattering::parameters::detail::optical_phonon_energy>
      (parameters)
   );
   scattering::inelastic_isotropic_transition<parameter_type>
   optical_emission_transition(
      parameters, 
      -1. * at_key<scattering::parameters::detail::optical_phonon_energy>
      (parameters)
   );
   scattering::polar_optic_transition<parameter_type>
   polar_optical_absorption_transition(
      parameters
      , at_key<scattering::parameters::detail::polar_optical_phonon_energy>
      (parameters)
   );
   scattering::polar_optic_transition<parameter_type> 
   polar_optical_emission_transition(
      parameters
      , -1. 
      * at_key<scattering::parameters::detail::polar_optical_phonon_energy>
      (parameters)
   );

   scattering::polar_optic_transition<parameter_type>
   polar_optical_absorption_transition2(
      parameters
      , at_key<scattering::parameters::detail::polar_optical_phonon_energy>
      (parameters) * .97
   );
   scattering::polar_optic_transition<parameter_type> 
   polar_optical_emission_transition2(
      parameters
      , -1.
      * at_key<scattering::parameters::detail::polar_optical_phonon_energy>
      (parameters) * .97
   );

   scattering::polar_optical_phonon_absorption_rate<parameter_type>
   polar_optical_absorption_rate(parameters);
   scattering::polar_optical_phonon_emission_rate<parameter_type>
   polar_optical_emission_rate(parameters);

   scattering::optical_phonon_absorption_rate<parameter_type>
   optical_absorption_rate(parameters);
   scattering::optical_phonon_emission_rate<parameter_type>
   optical_emission_rate(parameters);

   scattering_model_container<parameter_type>::type scattering_models(
      scattering_model_generate(acoustic_rate, acoustic_transition)
      ,
      scattering_model_generate(polar_optical_absorption_rate
                                , polar_optical_absorption_transition)
      , 
      scattering_model_generate(polar_optical_emission_rate
                                , polar_optical_emission_transition)
   );

   scattering_model_account<parameter_type>::type accounting(
      boost::fusion::make_pair<
      scattering_model<
      scattering::acoustic_phonon_rate<parameter_type>
      , scattering::elastic_isotropic_transition<parameter_type> > >(0)
      ,
      boost::fusion::make_pair<
      scattering_model<
      scattering::optical_phonon_emission_rate<parameter_type>
      , scattering::inelastic_isotropic_transition<parameter_type> > >(0)
      ,
      boost::fusion::make_pair<
      scattering_model<
      scattering::optical_phonon_absorption_rate<parameter_type>
      , scattering::inelastic_isotropic_transition<parameter_type> > >(0)
      ,
      boost::fusion::make_pair<
      scattering_model<
      scattering::polar_optical_phonon_absorption_rate<parameter_type>
      , scattering::polar_optic_transition<parameter_type> > >(0)
      ,
      boost::fusion::make_pair<
      scattering_model<
      scattering::polar_optical_phonon_emission_rate<parameter_type>
      , scattering::polar_optic_transition<parameter_type> > >(0)
   );
   std::vector<std::vector<double> > 
   k_integrate(cfg.max_iterations, std::vector<double>(cfg.x_count, 0));
   std::vector<std::vector<double> >
   k_integrate_positive(cfg.max_iterations
                        , std::vector<double>(cfg.x_count, 0)
   );
   std::vector<std::vector<double> > 
   k_integrate_negative(cfg.max_iterations
                        , std::vector<double>(cfg.x_count, 0)
   );

   std::vector<std::vector<double> >
   x_integrate(cfg.max_iterations, std::vector<double>(cfg.y_count, 0));
   std::vector<std::vector<double> >
   x_integrate_positive(cfg.max_iterations
                        , std::vector<double>(cfg.y_count, 0));
   std::vector<std::vector<double> >
   x_integrate_negative(cfg.max_iterations
                        , std::vector<double>(cfg.y_count, 0));

   std::vector<size_t> x_contributions;
   std::vector<size_t> k_contributions;
   std::vector<size_t> total_contributions;

   std::vector<size_t> k_outside(cfg.max_iterations);
   std::vector<size_t> x_outside(cfg.max_iterations);
   std::vector<size_t> total_outside(cfg.max_iterations);

   std::cout << "grid size: " << cfg.x_count * cfg.y_count << std::endl;

   std::vector<size_t> ky(cfg.y_count,0);
   std::vector<size_t> kz(cfg.y_count,0);

   {
      size_t ticker = 0;
      std::ofstream x_condensed(std::string("condensedx_")
                                + boost::lexical_cast<std::string>(ticker));
      for (size_t count_y(0) ; count_y < cfg.y_count; ++count_y)
      {
         x_condensed << (count_y * 1. /cfg.y_count * cfg.p_max 
                         - cfg.p_max * 0.5) ;

         double sum(0);
         for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
         {
            sum += timed_output_container[ticker][count_x][count_y];
         }
         x_condensed << "\t" << sum / cfg.max_particle_count << std::endl;
      }
      x_condensed.close();

   }

   typedef book<
      decltype( entropy_store ),
      decltype( k_integrate ),
      decltype( x_integrate ),
      decltype( k_outside ),
      decltype( x_outside ),
      decltype( k_integrate_negative ),
      decltype( k_integrate_positive ),
      decltype( x_integrate_negative ),
      decltype( x_integrate_positive ),
      decltype( timed_output_container ),
      decltype( timed_cancellation_map ),
      decltype( total_outside ),
      decltype( entropy ),
      decltype( participating )
      > book_type ;

   book_type books(
      entropy_store,
      k_integrate,
      x_integrate,
      k_outside,
      x_outside,
      k_integrate_negative,
      k_integrate_positive,
      x_integrate_negative,
      x_integrate_positive,
      timed_output_container,
      timed_cancellation_map,
      total_outside,
      entropy,
      participating
   ) ;

   work(data_store
        , parameters
        , cfg.timestep
        , accounting
        , scattering_models
        , books
        , cfg
        , flight_time
        , cfg.max_iterations
        , cfg.max_particle_count
   ) ;

   std::ofstream entropy_file(std::string("entropy")); 
   for (size_t ticker(0); ticker < cfg.max_iterations; ++ticker)
   {
      std::ofstream outputfile(std::string("propagation_")
                               + boost::lexical_cast<std::string>(ticker));
      for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
      {
         for (size_t count_y(0) ; count_y < cfg.y_count; ++count_y)
         {
            outputfile << timed_output_container[ticker][count_x][count_y]
               / cfg.max_particle_count << " ";
         }
         outputfile << std::endl;
      }
      outputfile.close();

      std::ofstream cancellation_outputfile(std::string("cancellation_map_") 
                                            + boost::lexical_cast<std::string>
                                            (ticker));
      for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
      {
         for (size_t count_y(0) ; count_y < cfg.y_count; ++count_y)
         {
            cancellation_outputfile 
               << timed_cancellation_map[ticker][count_x][count_y]
               / cfg.max_particle_count
               << " ";
         }
         outputfile << std::endl;
      }
      cancellation_outputfile.close();


      std::ofstream k_condensed(std::string("condensedk_") 
                                + boost::lexical_cast<std::string>(ticker));
      for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
      {
         k_condensed << (count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L);

         double sum(0);
         for (size_t count_y(0) ; count_y < cfg.y_count; ++count_y)
         {
            sum += timed_output_container[ticker][count_x][count_y];
         }
         k_condensed << "\t" << sum / cfg.max_particle_count  << std::endl;
      }
      k_condensed.close();

      std::ofstream x_condensed(std::string("condensedx_")
                                + boost::lexical_cast<std::string>(ticker));
      for (size_t count_y(0) ; count_y < cfg.y_count; ++count_y)
      {
         x_condensed << (count_y * 1. /cfg.y_count * cfg.p_max
                         - cfg.p_max * 0.5) ;

         double sum(0);
         for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
         {
            sum += timed_output_container[ticker][count_x][count_y];
         }
         x_condensed << "\t" << sum / cfg.max_particle_count << std::endl;
      }
      x_condensed.close();

      std::ofstream x_integrated(std::string("integratedx_")
                                 + boost::lexical_cast<std::string>(ticker));
      for (size_t count_y(0) ; count_y < cfg.y_count; ++count_y)
      {
         x_integrated << (count_y * 1. /cfg.y_count * cfg.p_max
                          - cfg.p_max * 0.5) 
                      << "\t"
                      << x_integrate[ticker][count_y]
            / cfg.max_particle_count
                      << std::endl;
      }
      x_integrated.close();

      std::ofstream k_integrated(std::string("integratedk_")
                                 + boost::lexical_cast<std::string>(ticker));
      for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
      {
         k_integrated << (count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L)
                      << "\t"
                      << k_integrate[ticker][count_x]
            / cfg.max_particle_count
                      << std::endl;
      }
      k_integrated.close();

      {
         std::ofstream x_positive(std::string("positivex_")
                                  + boost::lexical_cast<std::string>(ticker));
         for (size_t count_y(0) ; count_y < cfg.y_count; ++count_y)
         {
            x_positive << (count_y * 1. /cfg.y_count * cfg.p_max
                           - cfg.p_max * 0.5)
                       << "\t"
                       << x_integrate_positive[ticker][count_y]
               / cfg.max_particle_count
                       << std::endl;
         }
         x_positive.close();

         std::ofstream k_positive(std::string("positivek_")
                                  + boost::lexical_cast<std::string>(ticker));
         for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
         {
            k_positive << (count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L)
                       << "\t"
                       << k_integrate_positive[ticker][count_x]
               / cfg.max_particle_count
                       << std::endl;
         }
         k_positive.close();

      }

      {
         std::ofstream x_negative(std::string("negativex_")
                                  + boost::lexical_cast<std::string>(ticker));
         for (size_t count_y(0) ; count_y < cfg.y_count; ++count_y)
         {
            x_negative << (count_y * 1. /cfg.y_count * cfg.p_max
                           - cfg.p_max * 0.5)
                       << "\t"
                       << x_integrate_negative[ticker][count_y] 
               / cfg.max_particle_count
                       << std::endl;
         }
         x_negative.close();

         std::ofstream k_negative(std::string("negativek_")
                                  + boost::lexical_cast<std::string>(ticker));
         for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
         {
            k_negative << (count_x * 1. /cfg.x_count * cfg.L - 0.5 * cfg.L)
                       << "\t"
                       << k_integrate_negative[ticker][count_x]
               / cfg.max_particle_count
                       << std::endl;
         }
         k_negative.close();
      }


      {
         double global_entropy(0);
         for (entropy_store_type::const_iterator 
                 it(entropy_store[ticker].begin()); 
              it != (entropy_store[ticker].end());
              ++it)
         {
            global_entropy += (it->second) * (it->second);
         }

         double entrpy(0);
         for (size_t count_x(0) ; count_x < cfg.x_count; ++count_x)
         {
            for (size_t count_y(0) ; count_y < cfg.y_count; ++count_y)
            {
               entrpy += timed_output_container[ticker][count_x][count_y]
                  * timed_output_container[ticker][count_x][count_y];
            }
         }
         entropy_file << ticker << "\t" << 1. * ticker * cfg.timestep 
                      << "\t" << entrpy / cfg.max_particle_count 
                      << "\t" << entropy[ticker] / cfg.max_particle_count 
                      << "\t" << global_entropy / cfg.max_particle_count
                      << "\t" << participating[ticker]
                      << std::endl;

      }
   }
   entropy_file.close();

   std::cout << "k outside" << std::endl;
   for (std::vector<size_t>::iterator it(k_outside.begin());
        it!= k_outside.end();
        ++it)
   {
      std::cout << *it << " " ;
   }
   std::cout << std::endl;

   std::cout << "x outside" << std::endl;
   for (std::vector<size_t>::iterator it(x_outside.begin());
        it!= x_outside.end();
        ++it)
   {
      std::cout << *it << " " ;
   }
   std::cout << std::endl;

   std::cout << "total outside" << std::endl;
   for (std::vector<size_t>::iterator it(total_outside.begin());
        it!= total_outside.end();
        ++it)
   {
      std::cout << *it << " " ;
   }
   std::cout << std::endl;

   std::cout << "total transition accounting: " << std::endl;

   std::cout << "acoustic: " << std::endl;
   std::cout << at_key< 
      scattering_model<
         scattering::acoustic_phonon_rate<parameter_type>
         , scattering::elastic_isotropic_transition<parameter_type> > > 
   (accounting)
      << std::endl;

   std::cout << "polar optical absorption:" << std::endl;
   std::cout << at_key<
      scattering_model<
         scattering::polar_optical_phonon_absorption_rate<parameter_type>
         , scattering::polar_optic_transition<parameter_type> > > 
   (accounting) 
      << std::endl;
   std::cout << "polar optical emission:" << std::endl;
   std::cout << at_key<
      scattering_model<
         scattering::polar_optical_phonon_emission_rate<parameter_type>
         , scattering::polar_optic_transition<parameter_type> > > 
   (accounting) 
      << std::endl;


   return 0;
}
