#include <string>
#include <vector>
#include <fstream>

#include "json/json.h"
#include "gerlumph.hpp"



std::vector< std::vector<double> > readTimeVectors(const std::string filename){
  Json::Value root;
  std::ifstream fin(filename);
  fin >> root;

  std::vector< std::vector<double> > out_times;
  Json::Value times = root["times"];
  for(int j=0;j<times.size();j++){
    std::vector<double> out_time;
    Json::Value time = times[j];
    for(int j=0;j<time.size();j++){
      out_time.push_back( time[j].asDouble() );
    }
    out_times.push_back(out_time);
  }

  return out_times;
}

struct genericParameters {
  std::vector<std::string> ids;
  std::vector<double> lrest;
  double Rein;
  double vtot;
  int Nlc;
  int seed;
  bool full_data;
  bool degraded_data;
};

genericParameters readGenericParameters(const std::string filename){
  Json::Value root;
  std::ifstream fin(filename);
  fin >> root;

  genericParameters pars;
  
  Json::Value maps = root["maps"];
  for(int j=0;j<maps.size();j++){
    pars.ids.push_back( maps[j]["id"].asString() );
  }

  Json::Value lrest = root["lrest"];
  for(int j=0;j<lrest.size();j++){
    pars.lrest.push_back( lrest[j].asDouble() );
  }

  pars.Rein          = root["Rein"].asDouble();
  pars.vtot          = root["vtot"].asDouble();
  pars.Nlc           = root["output"]["Nlc"].asInt();
  pars.seed          = root["output"]["seed"].asInt();
  pars.full_data     = root["output"]["full_data"].asBool();
  pars.degraded_data = root["output"]["degraded_data"].asBool();
  
  return pars;
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
      root["profile"]["s0"].asDouble(),
      root["profile"]["l0"].asDouble(),
      root["profile"]["n"].asDouble()
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
      profile_pars.filename = root["profile"]["filenames"][j].asString();
      profile_pars.profPixSizePhys = root["profile"]["profPixSizePhys"][j].asDouble();
    }
    profile_pars_vector.push_back(profile_pars);
  }
    
  return profile_pars_vector;
}
