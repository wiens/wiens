#ifndef SCATTERING_MODEL
#define SCATTERING_MODEL


template<typename RateType, typename TransitionType>
struct scattering_model
{
   RateType rate;
   TransitionType transition;

   scattering_model(const RateType& rate, const TransitionType& transition) :
      rate(rate), transition(transition)
   {}

   template<typename StateType>
   inline typename RateType::return_type operator()(const StateType& state) 
      const
   {
      return rate(state);
   }

   template<typename StateType>
   inline StateType apply(const StateType& state) const
   {
      return transition(state);
   }
};


template<typename RateType, typename TransitionType>
inline scattering_model<RateType, TransitionType> 
scattering_model_generate(const RateType& rate
                          , const TransitionType& transition)
{
   return scattering_model<RateType, TransitionType>(rate, transition);
}

#endif // SCATTERING_MODEL
