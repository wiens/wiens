#ifndef SCATTERING_TAGS_HPP
#define SCATTERING_TAGS_HPP

namespace scattering {
   namespace constants {
      struct hbar ;
      struct k_B ;
      struct e ;
   }
   
   namespace parameters {
      
      struct TL ;
      struct mass ;
      struct alpha ;
      struct rho ;
      
      namespace detail {
         struct NI ;
         struct Z ;
         struct energy_beta ;

         struct DA ;
         struct us ;

         struct DO ;
         struct NO ;
         struct optical_phonon_energy ;

         struct polar_optical_phonon_energy ;
         struct optical_permitivity ;
         struct static_permitivity ;
         struct NLO ;

      }
   }
}


namespace flow
{

   struct free_flight_coeff_tag;

}

struct random_source_tag;


#endif // #define SCATTERING_TAGS_HPP
