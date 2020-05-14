#include "auxiliary_functions.hpp"

#include <fstream>
#include <cmath>
#include "json/json.h"

lsstParameters::lsstParameters(const std::string filename){
  Json::Value root;
  std::ifstream fin(filename);
  fin >> root;
  
  this->tmax = root["lsst"]["years"].asDouble()*365.25; // maximum length of observations in days
  
  // read the dates file for each filter
  std::string path_2_dates = root["path_2_dates"].asString();
  Json::Value filters = root["filters"];
  std::string dum;
  double c1,c2;
  for(int j=0;j<filters.size();j++){
    this->filters.push_back(filters[j].asString());
    std::vector<double> time;
    std::vector<double> depth;
    std::ifstream date_file( path_2_dates + filters[j].asString() + "_dates.dat" );
    std::getline(date_file,dum);
    while( date_file >> c1 >> c2 ){
      time.push_back(c1);
      depth.push_back(c2);
    }
    date_file.close();
    this->times.push_back(time);
    this->depths.push_back(depth);
  }
  
  // get the tmin and rescale times
  double tmin = 2820735; // a JD far in the future
  for(int j=0;j<filters.size();j++){
    if( times[j].size() > 0 ){
      if( times[j][0] < tmin ){
	tmin = times[j][0];
      }
    }
  }    
  for(int j=0;j<filters.size();j++){
    for(int i=0;i<times[j].size();i++){
      times[j][i] -= tmin;
    }
  }
  this->tmin = tmin;
  
  for(int i=0;i<this->filters.size();i++){
    this->errbase[i] = root["lsst"]["errbase"+this->filters[i]].asDouble();
  }
}

genericParameters::genericParameters(const std::string filename){
  Json::Value root;
  std::ifstream fin(filename);
  fin >> root;
  
  Json::Value maps = root["maps"];
  for(int j=0;j<maps.size();j++){
    this->ids.push_back( maps[j]["id"].asString() );
    this->mass.push_back( maps[j]["mass"].asDouble() );
    this->gamma_angle.push_back( maps[j]["gamma_angle"].asDouble() + 90.0 ); // This is the orientation of the shear, i.e. 2*phi, converted directly from Position Angle (east of north) to a normal cartesian system
  }
  
  Json::Value lrest = root["lrest"];
  for(int j=0;j<lrest.size();j++){
    this->lrest.push_back( lrest[j].asDouble() );
  }
  
  this->Nlc           = root["output"]["numlc"].asInt();
  this->seed          = root["output"]["seed"].asInt();
  this->velocities    = root["output"]["velocities"].asBool();
  this->start_end     = root["output"]["start_end"].asBool();
  this->path_2_output = root["path_2_output"].asString();
}

velocityParameters::velocityParameters(const std::string filename){
  Json::Value root;
  std::ifstream fin(filename);
  fin >> root;
  Json::Value vel = root["system"];
  
  this->ra         = vel["ra"].asDouble();
  this->dec        = vel["dec"].asDouble();
  this->sigma_l    = vel["sigma_l"].asDouble();
  this->sigma_s    = vel["sigma_s"].asDouble();
  this->sigma_disp = vel["sigma_disp"].asDouble();
  this->epsilon    = vel["epsi"].asDouble();
  this->zl         = vel["zl"].asDouble();
  this->zs         = vel["zs"].asDouble();
  this->Dl         = vel["Dl"].asDouble();
  this->Ds         = vel["Ds"].asDouble();
  this->Dls        = vel["Dls"].asDouble();
}


double m52snr(double dm){
  // (dm is the difference between mag and mag 5 sigma of, the depth)
  // find the SNR for a star of magnitude m obsreved
  // under conditions of 5-sigma limiting depth m5.  This assumes
  // Gaussianity and might not be strictly true in bluer filters.
  // see table 2 and eq 5 in astroph/0805.2366
  double snr = 5.0*pow(-0.4*dm,10);
  return 2.5*log10(1.0+1.0/snr);
}

std::vector<double> calculateRhalf(const std::string filename){
  Json::Value root;
  std::ifstream fin(filename);
  fin >> root;

  std::vector<double> rhalfs;
  std::string profile_type = root["profile"]["type"].asString();
  if( profile_type == "parametric" ){
    double r0 = root["profile"]["r0"].asDouble();
    for(int j=0;j<root["lrest"].size();j++){
      double r = r0*pow(root["lrest"][j].asDouble()/root["profile"]["l0"].asDouble(),root["profile"]["nu"].asDouble());
      rhalfs.push_back(r);
    }
  } else if( profile_type == "ssdisc" ){
    double b = pow(root["profile"]["mbh"].asDouble(),2.0);
    double c = root["profile"]["fedd"].asDouble()/root["profile"]["eta"].asDouble();
    for(int j=0;j<root["lrest"].size();j++){
      double a = pow(root["lrest"][j].asDouble(),4.0);
      double r = 0.0097*pow(a*b*c,1.0/3.0); // in [10^14 cm]
      rhalfs.push_back(r);
    }     
  } else {
    for(int j=0;j<root["lrest"].size();j++){
      rhalfs.push_back(0);
    }    
  }

  return rhalfs;
}


std::vector<BaseProfile*> createProfilesFromInput(const std::string filename,double pixSizePhys){
  Json::Value root;
  std::ifstream fin(filename);
  fin >> root;
  
  std::string profile_type  = root["profile"]["type"].asString();
  std::string profile_shape = root["profile"]["shape"].asString();
  double incl   = root["profile"]["incl"].asDouble();
  double orient = root["profile"]["orient"].asDouble();
  Json::Value lrest  = root["lrest"];
  std::vector<BaseProfile*> profiles(lrest.size());

  if( profile_type != "custom" ){
    if( profile_type == "parametric" ){
      for(int j=0;j<lrest.size();j++){
	double l = lrest[j].asDouble();
	double rhalf = root["profile"]["r0"].asDouble()*pow(l/root["profile"]["l0"].asDouble(),root["profile"]["nu"].asDouble());
	if( profile_shape == "uniform" ){
	  profiles[j] = new UniformDisc(pixSizePhys,rhalf/0.707,incl,orient);
	  //profiles[j] = new UniformDisc(pixSizePhys,rhalf,incl,orient);
	} else if( profile_shape == "gaussian" ){
	  profiles[j] = new Gaussian(pixSizePhys,rhalf/1.18,incl,orient);
	  //profiles[j] = new Gaussian(pixSizePhys,rhalf,incl,orient);
	}
      }
    } else if( profile_type == "ssdisc" ){
      double b = pow(root["profile"]["mbh"].asDouble(),2.0);
      double c = root["profile"]["fedd"].asDouble()/root["profile"]["eta"].asDouble();
      for(int j=0;j<lrest.size();j++){
	double l  = lrest[j].asDouble();
	double a = pow(l,4.0);
	double rhalf = 0.0097*pow(a*b*c,1.0/3.0); // in [10^14 cm]
	if( profile_shape == "uniform" ){
	  profiles[j] = new UniformDisc(pixSizePhys,rhalf/0.707,incl,orient);
	} else if( profile_shape == "gaussian" ){
	  profiles[j] = new Gaussian(pixSizePhys,rhalf/1.18,incl,orient);
	}
      }
    }
  } else {
    factoryProfilePars custom_pars;
    custom_pars.type   = profile_type;
    custom_pars.shape  = profile_shape;
    custom_pars.incl   = incl;
    custom_pars.orient = orient;
    custom_pars.pixSizePhys = pixSizePhys;
    for(int j=0;j<lrest.size();j++){
      custom_pars.lrest = lrest[j].asDouble();
      custom_pars.filename = root["path_2_custom"].asString() + root["filters"][j].asString() + ".fits";
      custom_pars.profPixSizePhys = root["profile"]["profPixSizePhys"][j].asDouble();
      profiles[j] = FactoryProfile::getInstance()->createProfile(custom_pars);	
    }
  }

  return profiles;
}



std::vector<factoryProfilePars> createProfileParsFromInput(const std::string filename){
  Json::Value root;
  std::ifstream fin(filename);
  fin >> root;

  std::vector<factoryProfilePars> profile_pars_vector;
  factoryProfilePars profile_pars;
  profile_pars.type   = root["profile"]["type"].asString();
  profile_pars.shape  = root["profile"]["shape"].asString();
  profile_pars.incl   = root["profile"]["incl"].asDouble();
  profile_pars.orient = root["profile"]["orient"].asDouble();

  if( profile_pars.type == "parametric" ){
    profile_pars.pars_parametric = {
      root["profile"]["r0"].asDouble(),
      root["profile"]["l0"].asDouble(),
      root["profile"]["nu"].asDouble()
    };
  } else if( profile_pars.type == "ssdisc" ){
    profile_pars.pars_parametric = {
      root["profile"]["mbh"].asDouble(),
      root["profile"]["fedd"].asDouble(),
      root["profile"]["eta"].asDouble()
    };
  }

  Json::Value lrest = root["lrest"];
  for(int j=0;j<lrest.size();j++){
    profile_pars.lrest = root["lrest"][j].asDouble();
    if( profile_pars.type == "custom" ){
      profile_pars.filename = root["path_2_custom"].asString() + root["filters"][j].asString() + ".fits";
      profile_pars.profPixSizePhys = root["profile"]["profPixSizePhys"][j].asDouble();
    }
    profile_pars_vector.push_back(profile_pars);
  }
    
  return profile_pars_vector;
}


void writeUncompressedData(std::string path,lsstParameters lsst,LightCurveCollection& mother,const std::vector<LightCurveCollection>& full,const std::vector<LightCurveCollection>& sampled){
  // Write sampled curves
  for(int j=0;j<lsst.Nfilters;j++){
    LightCurveCollection sample;
    sample.Ncurves = mother.Ncurves;
    sample.lightCurves = new LightCurve*[sample.Ncurves];
    sample.A = (point*) malloc(sample.Ncurves*sizeof(point)); // just because the destructor complains when there is no pointer to A or B to delete
    sample.B = (point*) malloc(sample.Ncurves*sizeof(point));
    for(int i=0;i<mother.Ncurves;i++){
      int Nsamples = sampled[j].lightCurves[i]->Nsamples;
      sample.lightCurves[i] = new LightCurve(Nsamples);
      for(int k=0;k<Nsamples;k++){
	sample.lightCurves[i]->t[k]  = lsst.tmin + sampled[j].lightCurves[i]->t[k];
	sample.lightCurves[i]->m[k]  = lsst.errbase[j] - 2.5*log10(sampled[j].lightCurves[i]->m[k]);
	sample.lightCurves[i]->dm[k] = m52snr(sample.lightCurves[i]->m[k]-lsst.depths[j][k]);
      }
    }
    sample.writeCurves(path,"table"+lsst.filters[j]+"_");
  }
  
  // Write theoretical curves
  for(int i=0;i<mother.Ncurves;i++){
    std::string full_file_name = path + "tablet_" + std::to_string(i) + ".dat";
    FILE* fh = fopen(full_file_name.c_str(),"w");
    for(int k=0;k<full[0].lightCurves[i]->Nsamples;k++){
      fprintf(fh," %11.6e",(double)k);
      for(int j=0;j<lsst.Nfilters;j++){
	fprintf(fh," %11.6e",lsst.errbase[j]-2.5*log10(full[j].lightCurves[i]->m[k]));
      }
      fprintf(fh,"\n");
    }
    fclose(fh);
  }
}

void writeCompressedData(std::string path,lsstParameters lsst,LightCurveCollection& mother,const std::vector<LightCurveCollection>& full,const std::vector<LightCurveCollection>& sampled){
  std::cout << "creating collection with all filters per light curve" << std::endl;
  std::vector<int> Nfull(mother.Ncurves);
  LightCurveCollection all_filters_full = mother;
  for(int i=0;i<mother.Ncurves;i++){
    Nfull[i] = full[0].lightCurves[i]->Nsamples;
    int Ntot = Nfull[i]*lsst.Nfilters;
    
    all_filters_full.lightCurves[i] = new LightCurve(Ntot);
    int start = 0;
    for(int j=0;j<lsst.Nfilters;j++){
      //      for(int k=start;k<Ntot;k++){
      for(int k=start;k<start+Nfull[i];k++){
	all_filters_full.lightCurves[i]->m[k] = lsst.errbase[j] - 2.5*log10( full[j].lightCurves[i]->m[k-start] );
      }
      start += Nfull[i];
    }
  }
  all_filters_full.writeCurvesDegraded<unsigned char>(path,"");
  


  std::string line1,line2,line3,line4,line5;
  std::vector<std::string> theo_limits;
  for(int i=0;i<mother.Ncurves;i++){
    std::ifstream file(path+"comp_p_"+std::to_string(i)+".dat");
    std::getline(file,line1);
    std::getline(file,line2);
    std::getline(file,line3);
    theo_limits.push_back(line3);
    file.close();
  }



  std::vector<int> Nf(lsst.Nfilters);
  for(int j=0;j<lsst.Nfilters;j++){
    Nf[j] = lsst.times[j].size();
  }
  int Ntot = std::accumulate(Nf.begin(),Nf.end(),0);
  LightCurveCollection all_filters_sampled = mother;
  for(int i=0;i<mother.Ncurves;i++){
    //    std::cout << " >>>>>>>>>>>>>> " << i << std::endl;
    all_filters_sampled.lightCurves[i] = new LightCurve(Ntot);
    int start = 0;
    for(int j=0;j<lsst.Nfilters;j++){
      //      std::cout << " >>>>>>>>>>>>>> filter " << j << std::endl;
      for(int k=0;k<Nf[j];k++){
	all_filters_sampled.lightCurves[i]->t[start + k]  = lsst.tmin + sampled[j].lightCurves[i]->t[k];
	all_filters_sampled.lightCurves[i]->m[start + k]  = lsst.errbase[j] - 2.5*log10( sampled[j].lightCurves[i]->m[k] );
	all_filters_sampled.lightCurves[i]->dm[start + k] = m52snr( all_filters_sampled.lightCurves[i]->m[start + k] - lsst.depths[j][k] );
	//	std::cout << all_filters_sampled.lightCurves[i].m[k] << std::endl;
      }
      start += Nf[j];
    }
  }
  all_filters_sampled.writeCurvesDegraded<unsigned char,unsigned short int,unsigned char>(path,"");
  //all_filters_sampled.writeCurvesDegraded<unsigned short int,unsigned short int,unsigned short int>(path,"");
  


  // Read and overwrite params.dat files
  std::string basic;
  for(int j=0;j<lsst.Nfilters;j++){
    basic += " " + std::to_string(Nf[j]);
  }
  for(int i=0;i<mother.Ncurves;i++){
    std::ifstream file(path+"comp_p_"+std::to_string(i)+".dat");
    std::getline(file,line1);
    std::getline(file,line2);
    std::getline(file,line3);
    std::getline(file,line4);
    std::getline(file,line5);
    line1 = basic + " " + std::to_string(Nfull[i]*lsst.Nfilters);
    line3 = theo_limits[i];
    file.close();
    
    std::ofstream out(path+"comp_p_"+std::to_string(i)+".dat");
    out << line1 << std::endl << line2 << std::endl << line3 << std::endl << line4 << std::endl << line5 << std::endl;
    out.close();
  }


}

void writeUncompressedDataNew(std::string path,lsstParameters lsst,LightCurveCollection& mother,const std::vector<LightCurveCollection>& full,const std::vector<LightCurveCollection>& sampled){

  // Write sampled curves
  for(int j=0;j<lsst.Nfilters;j++){
    std::string filter_mag = path + "filter_"+lsst.filters[j]+"_mag.dat";
    FILE* fh_filter_mag = fopen(filter_mag.c_str(),"w");
    std::string filter_dmag = path + "filter_"+lsst.filters[j]+"_dmag.dat";
    FILE* fh_filter_dmag = fopen(filter_dmag.c_str(),"w");
    
    for(int i=0;i<mother.Ncurves;i++){
      for(int k=0;k<sampled[j].lightCurves[0]->Nsamples;k++){
	fprintf(fh_filter_mag,"%13.6e",lsst.errbase[j] - 2.5*log10(sampled[j].lightCurves[i]->m[k]));
	fprintf(fh_filter_dmag,"%13.6e",m52snr(sampled[j].lightCurves[i]->m[k]-lsst.depths[j][k]));
      }
      fprintf(fh_filter_mag,"\n");
      fprintf(fh_filter_dmag,"\n");
    }
  }

  // Write length of theoretical light curves
  std::string theo_length = path + "theo_length.dat";
  FILE* fh_len = fopen(theo_length.c_str(),"w");
  for(int i=0;i<mother.Ncurves;i++){
    fprintf(fh_len,"%d\n",full[0].lightCurves[i]->Nsamples);
  }

  // Write theoretical curves
  std::string theo_mag = path + "theo_mag.dat";
  FILE* fh_mag = fopen(theo_mag.c_str(),"w");
  for(int i=0;i<mother.Ncurves;i++){
    for(int k=0;k<full[0].lightCurves[i]->Nsamples;k++){
      for(int j=0;j<lsst.Nfilters;j++){
	fprintf(fh_mag," %13.6e",lsst.errbase[j]-2.5*log10(full[j].lightCurves[i]->m[k]));
      }
      fprintf(fh_mag,"\n");
    }
  }
  fclose(fh_mag);

}
