//
//  sigModel.h
//  cygnoNeutrino
//
//  Created by Stefano Piacentini on 19/10/22.
//

#ifndef sigModel_h
#define sigModel_h

#include <BAT/BCModel.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Math/ProbFunc.h"

class sigMOD : public BCModel
{
public:
    
  sigMOD(const std::string& name, std::string datafile,
	 std::string bkgfile, std::string sigfile);
  
  ~sigMOD(){};
  
  double LogLikelihood(const std::vector<double>& pars);
  
  double LogAPrioriProbability(const std::vector<double>& pars);
  
  std::vector<double> getData(std::string datafile);
  std::vector<double> getBkg(std::string bkgfile);
  std::vector<double> getSig(std::string sigfile);
  std::map<int,double> getLogFactMap(std::string logfactfile);

  void SetModeEval(double a){fModeLogEval = a;}
  
private:
 long  double nbmax, nsmax;
  int nbins;
  
  std::vector<double> data;
  std::vector<double> bkg;
  std::vector<double> sig;  
  std::map<int,double> mapLogFactorial;
  
  
  double bkgNorm;
  double sigNorm;
  
  
  int nsigevts;
  int nbkgevts;
  
  double bkg_prior_lambda;    
  
  std::vector<double> getVectorFromFile(std::string filename);
  std::map<int,double> getMapFromFile(std::string filename);

  double fModeLogEval;
  
};

#endif


/* sigModel_h */
