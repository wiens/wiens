#ifndef TUPLE_STATE_HPP
#define TUPLE_STATE_HPP

namespace tuple_state {

template<typename T>
struct query_momentum_type;

template<typename S, typename T, typename U>
struct query_momentum_type<std::tuple<S,T,U> >
{
   typedef T type;
};

template<typename S, typename T, typename U>
typename query_momentum_type<std::tuple<S,T,U> >::type
query_momentum(const std::tuple<S,T,U>& state_tuple)
{
   return std::get<1>(state_tuple);
}

template<typename T>
struct query_position_type { typedef T type; };

template<typename S, typename T, typename U>
struct query_position_type<std::tuple<S,T,U> >
{
   typedef S type;
};

template<typename S, typename T, typename U>
typename query_position_type<std::tuple<S,T,U> >::type
query_position(const std::tuple<S,T,U>& state_tuple)
{
   return std::get<0>(state_tuple);
}

template<typename T>
T query_position(T t)
{
   return t ;
}

}

#endif // TUPLE_STATE_HPP
