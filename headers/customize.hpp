#ifndef CUSTOMIZE_HPP
#define CUSTOMIZE_HPP

#include <valarray>

namespace utility
{
   template<typename S, typename T>
   struct query_momentum_square_modulus_type<boost::units::quantity<S, T> >
   {
      typedef typename boost::units::multiply_typeof_helper<
         S,S
         >::type dimension_type;
      typedef boost::units::quantity<dimension_type, T> type;
   };

   template<typename S, typename T>
   struct query_momentum_square_modulus_type<
      boost::units::quantity<S, std::valarray<T> > >
   {
      typedef typename boost::units::multiply_typeof_helper<
         S,S
         >::type dimension_type;
      typedef typename std::valarray<T>::value_type condensed_type;
      
      typedef boost::units::quantity<dimension_type, condensed_type> type;
   };


   template<typename S, typename T> 
   struct momentum_square_modulus_calculate<
      boost::units::quantity<S, std::valarray<T> >
      >
   {
      typename query_momentum_square_modulus_type<
         boost::units::quantity<S, std::valarray<T> > >::type 
      apply(const boost::units::quantity<S, std::valarray<T> > & momentum)
      {
         typedef typename query_momentum_square_modulus_type<
         boost::units::quantity<S, std::valarray<T> > >::type
            return_type;
         return return_type::from_value(
            ((boost::units::quantity_cast<std::valarray<T> >
              (momentum * momentum)).sum()
            )
         );
      }
   };

   template<typename S, typename T> 
   typename query_momentum_square_modulus_type<
      boost::units::quantity<S, std::valarray<T> > >::type 
   momentum_square_modulus(const boost::units::quantity<S, std::valarray<T> > &
                           momentum)
   {
      typedef typename query_momentum_square_modulus_type<
         boost::units::quantity<S, std::valarray<T> > >
         ::type return_type;
      return return_type::from_value(
         ((boost::units::quantity_cast<std::valarray<T> >
           (momentum * momentum)).sum()
         )
      );
   }
}


template<typename StreamType, typename T>
StreamType& operator<< (StreamType& out, const std::valarray<T>& val)
{
   out << "[ " ;
   for (size_t i(0); i<val.size(); ++i)
      out << val[i] << ' ';
   out << ']' ;
   return out;
}

namespace std{

template<class Y>
std::ostream& operator<<(std::ostream& os,const std::valarray<Y>& val)
{
   os << "[ " ;
   for (size_t i(0); i<val.size(); ++i)
      os << val[i] << ' ';
   os << ']' ;
   return os;
}

}


#endif // CUSTOMIZE_HPP
