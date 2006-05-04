#ifndef H_MULTIREGIONS


namespace MultiRegions{


  // multiregion stuff here
  enum TransState{
    NotSet,      ///< No transformed state set 
    Local,       ///< Local  Modal space array contains "true" expansion values
    Continuous,  ///< Continuous Modal space array contains "true" expansion values
    LocalCont,   ///< Both local and continuous space array contains "true" Expansion values 
  };


}// end of namespace

#define H_MULTIREGIONS
#endif
