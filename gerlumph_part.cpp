#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <fstream>

#include "json/json.h"

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
  std::vector<double> rhalfs = calculateRhalf(json_input_filename);
  std::cout << "profile ok" << std::endl;


  // Set the velocities
  velocityParameters vp(json_input_filename);
  velocityComponents vel(gen.Nlc);
  vel.createVelocitiesK04(321,vp.ra,vp.dec,vp.sigma_l,vp.sigma_s,vp.sigma_disp,vp.epsilon,vp.zl,vp.zs,vp.Dl,vp.Ds,vp.Dls);
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
    //std::vector<BaseProfile*> profiles(lsst.Nfilters);
    //for(int j=0;j<lsst.Nfilters;j++){
    //  profile_pars[j].pixSizePhys = map.pixSizePhys;
    //  profiles[j] = FactoryProfile::getInstance()->createProfile(profile_pars[j]);
    //}
    std::vector<BaseProfile*> profiles = createProfilesFromInput(json_input_filename,map.pixSizePhys);


    // set convolution kernel
    double profMaxOffset = profiles[lsst.Nfilters-1]->Nx/2;
    EffectiveMap emap(profMaxOffset,&map);
    Kernel kernel(map.Nx,map.Ny);


    // Rotate map (i.e. rotate the velocity vectors in an opposite way)
    for(int i=0;i<gen.Nlc;i++){
      phi_vtot[i] -= gen.gamma_angle[i];
    }


    // Set light curves
    std::cout << "setting light curves" << std::endl;
    LightCurveCollection mother(gen.Nlc,&emap);
    mother.createVelocityLocations(213,lsst.tmax,vtot,phi_vtot);
    //mother.createRandomLocations(213,1000);
    //    mother.A[0].x = 0;
    //    mother.A[0].y = 0;
    //    mother.B[0].x = 10000 - profMaxOffset;
    //    mother.B[0].y = 10000 - profMaxOffset;


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
    std::cout << "writing uncompressed data" << std::endl;
    writeUncompressedDataNew(gen.path_2_output,lsst,mother,all_filters_full_raw,all_filters_sampled_raw);

    // Write velocities
    if( gen.velocities ){
      vel.writeVelocities(gen.path_2_output+"velocities.dat");
    }

    // Writing start and end points of each light curve on the magnification map
    if( gen.start_end ){
      std::cout << "writing start and end points" << std::endl;
      FILE* fh_points = fopen((gen.path_2_output+"xy_start_end.dat").c_str(),"w");
      for(int i=0;i<mother.Ncurves;i++){
	fprintf(fh_points,"%8.3f %8.3f %8.3f %8.3f\n",mother.A[i].x,mother.A[i].y,mother.B[i].x,mother.B[i].y);
      }
      fclose(fh_points);
    }

    // Write file with generic parameters, necessary to convert the x-values of the theoretical light curves to time
    Json::Value out;
    out["pixSizePhys"] = map.pixSizePhys; // in 10^14 cm
    out["Rein"]        = Rein; // in 10^14 cm
    out["time0"]       = lsst.tmin; // in MJD
    out["duration"]    = lsst.tmax; // in days
    out["offset"]      = profMaxOffset; // in pixels
    out["profile_type"] = profile_pars[0].type;
    Json::Value sizes = Json::Value(Json::arrayValue);
    for(int j=0;j<lsst.Nfilters;j++){
      sizes.append(rhalfs[j]);
    }
    out["r_half"] = sizes;
    
    std::ofstream jsonfile(gen.path_2_output+"parameters.json");
    jsonfile << out;
    jsonfile.close();
  }


  std::cout << "Done." << std::endl;
  return 0;
}
