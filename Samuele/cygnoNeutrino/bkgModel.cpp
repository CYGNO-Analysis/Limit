//
//  PMT_standard.cpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//

#include "bkgModel.h"
#include <TMath.h>
#include <BAT/BCMath.h>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <sstream>

bkgMOD::bkgMOD(const std::string& name,
               std::string datafile, std::string bkgfile) : BCModel(name),fModeLogEval(0)
{
    std::cout<<"Initializating BAT bkg model '"<<name<<"'"<<std::endl;
    
    data = getData(datafile);
    bkg  = getBkg(bkgfile);
    mapLogFactorial = getLogFactMap("mapLog.txt"); 
    
    bkgNorm = 0;
    for(unsigned int i =0; i<bkg.size(); i++) {
        bkgNorm += bkg[i];
    }
    
    nsigevts = data[0];
    nbkgevts = data[1];
    bkg_prior_lambda = nbkgevts;
    //bkg_prior_lambda = bkg_amount; // to be modified?
    nbmax = 5*bkg_prior_lambda;
    
    if(bkg.size()!=(data.size()-2)) {
        throw std::runtime_error("Length of data file and bkg file not matching.\n");
    }
    nbins = (int)data.size() - 2;
    
    std::cout<<"Number of bins: "<<nbins<<std::endl<<std::flush;
    
    //AddParameter("nb",nbkgevts-7*sqrt(nbkgevts),  nbkgevts+4*sqrt(nbkgevts)+5*nsigevts, "nb", "[counts]");
    AddParameter("nb",0, 3*(nbkgevts+nsigevts), "nb", "[counts]");
    
}

double bkgMOD::LogLikelihood(const std::vector<double>& pars) {
    

  double LL = 0.;
    
    for(int i = 0; i<nbins; i++) {
        double lambda_i = pars[0] * (bkg[i] / bkgNorm);
        LL += BCMath::LogPoisson((long long int)(data[i+2]), lambda_i);
	//LL += BCMath::ApproxLogFact(int(data[i+2]));
	//LL += mapLogFactorial[data[i+2]];
    }

    if(fModeLogEval!=0){
      LL-=fModeLogEval;
    }
    
    
    return LL;
}

double bkgMOD::LogAPrioriProbability(const std::vector<double>& pars) {
    double LL = 0.;
    
    // prior on nb
    LL += BCMath::LogPoisson(int(pars[0]), bkg_prior_lambda);
    //LL += BCMath::LogGaus(int(pars[0]), bkg_prior_lambda, 10 * sqrt(bkg_prior_lambda));
    
    return LL;
}


std::vector<double> bkgMOD::getData(std::string datafile) {
    std::vector<double> ret;
    
    ret = getVectorFromFile(datafile);
    
    return ret;
}


std::vector<double> bkgMOD::getBkg(std::string bkgfile) {
    std::vector<double> ret;
    
    ret = getVectorFromFile(bkgfile);
    
    return ret;
}

std::map<int,double> bkgMOD::getLogFactMap(std::string logfactfile){
  std::map<int,double> ret;

  ret = getMapFromFile(logfactfile);

  return ret;
}

std::vector<double> bkgMOD::getVectorFromFile(std::string filename) {
    std::vector<double> ret;
    
    
    
    std::ifstream myfile;
    myfile.open(filename);
    
    std::string line;
    
    if(myfile.is_open()) {
        while(std::getline(myfile, line)) {
     
            double tmp = stod(line);
            
            ret.push_back(tmp);
        }
    } else {
        throw std::runtime_error("Could not open the file\n");
    }
    myfile.close();
    
    //throw std::runtime_error("DEBUG"); //DEBUG
    
    return ret;
}



std::map<int,double> bkgMOD::getMapFromFile(std::string filename) {

  std::map<int,double> ret;
  
  std::ifstream myfile;
  myfile.open(filename);
  
  std::string line; 

  int temp1;
  double temp2;
  
  if(myfile.is_open()) {
        while(std::getline(myfile, line)) {

	  myfile>> temp1 >> temp2;

	  ret[temp1]=temp2;
	  
        }
    } else {
        throw std::runtime_error("Could not open the file\n");
    }
    myfile.close();

    return ret;
}
