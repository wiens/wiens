#ifndef MOMENTUM_UTILITY_HPP
#define MOMENTUM_UTILITY_HPP

namespace utility
{

// template<typename StateType>
// struct query_momentum_type
// {
//    typedef StateType type;
// };
// 
// template<typename StateType> 
// typename query_momentum_type<StateType>::type 
// query_momentum(const StateType& state)
// {
//    return state;
// }

template<typename T, typename S>
struct energy_dispersion_type
{
   typedef S type;
};

template<typename MomentumType>
struct query_momentum_square_modulus_type
{
   typedef MomentumType type;
};

template<typename MomentumType>
struct momentum_square_modulus_calculate
{
   typename query_momentum_square_modulus_type<MomentumType>::type 
   apply(const MomentumType& momentum)
   {
      return momentum * momentum ;
   }
};

template<typename MomentumType> 
typename query_momentum_square_modulus_type<MomentumType>::type 
momentum_square_modulus(const MomentumType& momentum)
{
   return momentum * momentum ;
}

}

#endif //  MOMENTUM_UTILITY_HPP
