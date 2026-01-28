//
//  bkgModel.cpp
//  cygnoNeutrino
//
//  Created by Stefano Piacentini on 19/10/22.
//

#ifndef bkgModel_h
#define bkgModel_h

#include <BAT/BCModel.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Math/ProbFunc.h"

class bkgMOD : public BCModel
{
public:
    
  bkgMOD(const std::string& name, std::string datafile,
	 std::string bkgfile);
  
  ~bkgMOD(){};
  
  double LogLikelihood(const std::vector<double>& pars);
  
  double LogAPrioriProbability(const std::vector<double>& pars);
  
  std::vector<double> getData(std::string datafile);
  std::vector<double> getBkg(std::string bkgfile);
  std::map<int,double> getLogFactMap(std::string logfactfile); 

  void SetModeEval(double a){fModeLogEval = a;} 

private:
  double nbmax;
  int nbins;
  
  int nsigevts;
  int nbkgevts;
  
  std::vector<double> data;
  std::vector<double> bkg;

  std::map<int,double> mapLogFactorial; 
  
  double bkgNorm;
  
  double bkg_prior_lambda;
  
  
  std::vector<double> getVectorFromFile(std::string filename);

  std::map<int,double> getMapFromFile(std::string filename);       

  double fModeLogEval; 
    
};

#endif /* bkgModel_h */
