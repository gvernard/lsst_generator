#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>

#include "gerlumph.hpp"
#include "auxiliary_functions.hpp"

int main(int argc,char* argv[]){
  std::string json_input_filename = "input_cpp.json";


  // Set number of LSST filters (6 of course!)
  int Nfilters = 6; 

  // Set generic parameters
  genericParameters gen = readGenericParameters(json_input_filename);
  double Rein = gen.Rein; // in [10^14 cm]
  double vtot = gen.vtot; // in [10^14 cm/day]
  std::vector<std::string> map_ids = gen.ids;
  std::vector<double> lrest = gen.lrest; // in [nm]

  // Set the time vectors from LSST strategy
  std::vector< std::vector<double> > times = readTimeVectors(json_input_filename);

  // Set the 6 LSST filter profiles
  std::vector<factoryProfilePars> profile_pars = createProfileParsFromInput(json_input_filename);





  // Starting convolutions and output
  for(int i=0;i<map_ids.size();i++){
    MagnificationMap map(map_ids[i],Rein);

    // I need to define pixSizePhys (map dependent) for each profile before creating it
    std::vector<BaseProfile*> profiles(Nfilters);
    for(int j=0;j<Nfilters;j++){
      profile_pars[j].pixSizePhys = map.pixSizePhys;
      profiles[j] = FactoryProfile::getInstance()->createProfile(profile_pars[j]);
    }

    double profMaxOffset = 1000; // dummy for the moment
    EffectiveMap emap(profMaxOffset,&map);
    Kernel kernel(map.Nx,map.Ny);

    for(int j=0;j<Nfilters;j++){
      kernel.setKernel(profiles[j]);
      map.convolve(&kernel,&emap);
    }

  }


  return 0;
}
