#ifndef AUXILIARY_HPP
#define AUXILIARY_HPP

#include <vector>
#include <string>
#include "gerlumph.hpp"

class lsstParameters {
public:
  int Nfilters = 6;
  std::vector< std::vector<double> > times;
  std::vector< std::vector<double> > depths;
  std::vector<std::string> filters;
  double tmin;
  double tmax;
  double errbase[6];

  lsstParameters(const std::string filename);
};

class genericParameters {
public:
  std::vector<std::string> ids;
  std::vector<double> mass;
  std::vector<double> lrest;    // in [nm]
  double Rein;
  int Nlc;
  int seed;
  bool full_data;
  bool degraded_data;
  bool velocities;
  std::string path_2_output;

  genericParameters(const std::string filename);
};

class velocityParameters {
public:
  double ra;         // in hours
  double dec;        // in deg
  double sigma_l;    // km/s
  double sigma_s;    // km/s
  double sigma_disp; // km/s
  double zl;
  double zs;
  double Dl;
  double Ds;
  double Dls;

  velocityParameters(const std::string filename);
};

double m52snr(double dm);
std::vector<factoryProfilePars> createProfileParsFromInput(const std::string filename);
void writeUncompressedData(std::string path,lsstParameters lsst,LightCurveCollection& mother,const std::vector<LightCurveCollection>& full,const std::vector<LightCurveCollection>& sampled);
void writeUncompressedDataNew(std::string path,lsstParameters lsst,LightCurveCollection& mother,const std::vector<LightCurveCollection>& full,const std::vector<LightCurveCollection>& sampled);
void writeCompressedData(std::string path,lsstParameters lsst,LightCurveCollection& mother,const std::vector<LightCurveCollection>& full,const std::vector<LightCurveCollection>& sampled);

#endif /* AUXILIARY_HPP */
