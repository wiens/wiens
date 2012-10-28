#ifndef EVOLUTION_TIME
#define EVOLUTION_TIME

#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/include/value_at_key.hpp>


struct random_source_tag;

namespace flow
{

struct free_flight_coeff_tag;

template<typename T>
double free_flight_coeff(T);

template <typename ParameterType>
struct evolution_time
{
   ParameterType& parameters;
   evolution_time(ParameterType& parameters) : parameters(parameters) {}

   typedef typename boost::fusion::result_of::value_at_key<
      ParameterType, 
      free_flight_coeff_tag>::type return_type;

   template<typename ParticleType>
   return_type operator()(const ParticleType& particle) 
   {
      return -log(boost::fusion::at_key<random_source_tag>(parameters)()) 
         * (boost::fusion::at_key<free_flight_coeff_tag>(parameters)) ;
   }
};

}

#endif // EVOLUTION_TIME
