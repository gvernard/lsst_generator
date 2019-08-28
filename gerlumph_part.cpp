#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <fstream>

#include "gerlumph.hpp"
#include "auxiliary_functions.hpp"

int main(int argc,char* argv[]){
  if( argc != 2 ){
    std::cout << "1 command line argument is required!" << std::endl;
    std::cout << "1st: the path to the json file containing all the parameters for the run." << std::endl;
    return 0;
  }
  std::string json_input_filename = argv[1];


  // Set LSST parameters (times and depths)
  lsstParameters lsst(json_input_filename);
  std::cout << "lsst ok" << std::endl;

  // Set generic parameters
  genericParameters gen(json_input_filename);
  std::cout << "gen ok" << std::endl;

  // Set the 6 LSST filter profiles
  std::vector<factoryProfilePars> profile_pars = createProfileParsFromInput(json_input_filename);
  std::cout << "profile ok" << std::endl;


  // Set the velocities
  velocityParameters vp(json_input_filename);
  velocityComponents vel(gen.Nlc);
  vel.createVelocitiesK04(321,vp.ra,vp.dec,vp.sigma_l,vp.sigma_s,vp.sigma_disp,vp.zl,vp.zs,vp.Dl,vp.Ds,vp.Dls);
  if( gen.velocities ){
    vel.writeVelocities(gen.path_2_output+"velocities.dat");
  }
  std::vector<double> vtot(gen.Nlc);
  std::vector<double> phi_vtot(gen.Nlc);
  for(int i=0;i<gen.Nlc;i++){
    vtot[i]     = vel.tot[i].v;
    phi_vtot[i] = vel.tot[i].phi;
    
  }
  std::cout << "velocities ok" << std::endl;




  // Map loop
  for(int ii=0;ii<gen.ids.size();ii++){
    double Rein = 13.5*sqrt(gen.mass[ii]*vp.Dls*vp.Ds/vp.Dl); // in 10^14 cm
    std::cout << "Rein is: " << Rein << " x10^14 cm" << std::endl;
    MagnificationMap map(gen.ids[ii],Rein);

    // Profile loop: I need to define pixSizePhys (map dependent) for each profile before creating it
    std::vector<BaseProfile*> profiles(lsst.Nfilters);
    for(int j=0;j<lsst.Nfilters;j++){
      profile_pars[j].pixSizePhys = map.pixSizePhys;
      profiles[j] = FactoryProfile::getInstance()->createProfile(profile_pars[j]);
    }

    double profMaxOffset = profiles[lsst.Nfilters-1]->Nx/2;
    EffectiveMap emap(profMaxOffset,&map);
    Kernel kernel(map.Nx,map.Ny);


    // Set light curves
    std::cout << "setting light curves" << std::endl;
    LightCurveCollection mother(gen.Nlc,&emap);
    mother.createVelocityLocations(213,lsst.tmax,vtot,phi_vtot);
    //mother.createRandomLocations(213,1000);



    // Convolution and extraction loop
    std::vector<LightCurveCollection> all_filters_full_raw;
    std::vector<LightCurveCollection> all_filters_sampled_raw;
    for(int j=0;j<lsst.Nfilters;j++){
      LightCurveCollection dum_full = mother;
      all_filters_full_raw.push_back(dum_full);
      LightCurveCollection dum_sampled = mother;
      all_filters_sampled_raw.push_back(dum_sampled);
    }

    std::cout << "starting convolutions" << std::endl;
    for(int j=0;j<lsst.Nfilters;j++){
      kernel.setKernel(profiles[j]);
      map.convolve(&kernel,&emap);

      all_filters_full_raw[j].extractFull();
      all_filters_sampled_raw[j].extractStrategy(vtot,lsst.times[j]);
    }
    std::cout << "convolutions done" << std::endl;


    // Output
    std::cout << "starting output" << std::endl;
    if( gen.full_data ){
      std::cout << "writing uncompressed data" << std::endl;
      //      writeUncompressedData(gen.path_2_output,lsst,mother,all_filters_full_raw,all_filters_sampled_raw);
      writeUncompressedDataNew(gen.path_2_output,lsst,mother,all_filters_full_raw,all_filters_sampled_raw);
    }
    if( gen.degraded_data ){
      std::cout << "writing degraded data" << std::endl;
      writeCompressedData(gen.path_2_output,lsst,mother,all_filters_full_raw,all_filters_sampled_raw);
    }
    FILE* fh = fopen((gen.path_2_output+"time_and_length.dat").c_str(),"w");
    fprintf(fh,"%15s: %13.6e [%15s]\n","pixel size",map.pixSizePhys,"x10^14cm");
    fprintf(fh,"%15s: %13.6e [%15s]\n","Rein",Rein,"x10^14cm");
    fprintf(fh,"%15s: %13.6e [%15s]\n","time0",lsst.tmin,"MJD");
    fprintf(fh,"%15s: %13.6e [%15s]\n","duration",lsst.tmax,"days");
    fclose(fh);


  }


  std::cout << "Done." << std::endl;
  return 0;
}
