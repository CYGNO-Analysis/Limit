//
//  sigModel.cpp
//  cygnoNeutrino
//
//  Created by Stefano Piacentini on 19/10/22.
//

#include <TMath.h>
#include <BAT/BCMath.h>
#include <cmath>

#include "sigModel.h"

sigMOD::sigMOD(const std::string& name,
               std::string datafile, std::string bkgfile,
               std::string  sigfile) : BCModel(name),fModeLogEval(0)
{
  std::cout<<"Initializating BAT bkg model '"<<name<<"'"<<std::endl;
    
    data = getData(datafile);
    bkg  = getBkg(bkgfile);
    sig  = getSig(sigfile);
    mapLogFactorial = getLogFactMap("mapLog.txt");
    
    bkgNorm = 0;
    for(unsigned int i =0; i<bkg.size(); i++) {
        bkgNorm += bkg[i];
    }
    
    sigNorm = 0;
    for(unsigned int i =0; i<sig.size(); i++) {
        sigNorm += sig[i];
    }
    
    nsigevts = data[0];
    nbkgevts = data[1];
    bkg_prior_lambda = nbkgevts;
    //bkg_prior_lambda = bkg_amount; // to be modified?

    nbmax = 5*bkg_prior_lambda; //number of bkg events constraind to 5 times the number of theoretical bkg events  

    //nbmax=2*nbkgevts;
    
    int normData = 0;
    
    for(unsigned int i = 2; i<data.size(); i++) {
        normData += data[i];
    }

    nsmax = normData; //number of signal events constraind to sum of events
    
    //nsmax=nsigevts*7;
    
    
    
    if(bkg.size()!=(data.size()-2)) {
        throw std::runtime_error("Length of data file and bkg file not matching.\n");
    } else if(sig.size()!=(data.size()-2)) {
        throw std::runtime_error("Length of data file and sig file not matching.\n");
    }
    nbins = (int)data.size() - 2;
    
    // DEBUG
    std::cout<<"Number of bins: "<<nbins<<std::endl<<std::flush;

    //a > b ? a : b
    
    //AddParameter("nb", nbkgevts-7*sqrt(nbkgevts),  nbkgevts+4*sqrt(nbkgevts)+5*nsigevts, "nb", "[counts]");  //number of bkg events constraind to 5 times the number of theoretical bkg events 
    //AddParameter("ns", 0, nsigevts+14*nsigevts, "ns", "[counts]");  //number of signal events constraind to sum of events

    AddParameter("nb", 0,  3*(nbkgevts+nsigevts), "nb", "[counts]");  //number of bkg events constraind to 5 times the number of theoretical bkg events 
    AddParameter("ns", 0, nbkgevts+nsigevts, "ns", "[counts]");  //number of signal events constraind to sum of events         
    
}

double sigMOD::LogLikelihood(const std::vector<double>& pars) {
    
    double LL = 0.;
    
    for(int i = 0; i<nbins; i++) {
        double lambda_i = pars[0] * (bkg[i] / bkgNorm) + pars[1] * (sig[i]/sigNorm);
        LL += BCMath::LogPoisson((long long int)(data[i+2]), lambda_i);
	//std::cout << LL << std::endl;

	//LL += BCMath::ApproxLogFact(int(data[i+2]));
	//LL += mapLogFactorial[data[i+2]];
    }

    
    if(fModeLogEval!=0){
      LL-=fModeLogEval;
    }
    
    
    return LL;
}

double sigMOD::LogAPrioriProbability(const std::vector<double>& pars) {
    
  double LL = 0.;
    
    // prior on nb
    LL += BCMath::LogPoisson(int(pars[0]), bkg_prior_lambda);
    //LL += BCMath::LogGaus(int(pars[0]), bkg_prior_lambda, 10 * sqrt(bkg_prior_lambda));
    
    // uniform prior on ns
    if(pars[1]<0 || pars[1]>nsmax) {
        LL += log(0.0);
    } else {
        LL += log(1.0/nsmax);
    }
    
    return LL;
}



std::vector<double> sigMOD::getData(std::string datafile) {
    std::vector<double> ret;
    
    ret = getVectorFromFile(datafile);
    
    return ret;
}


std::vector<double> sigMOD::getBkg(std::string bkgfile) {
    std::vector<double> ret;
    
    ret = getVectorFromFile(bkgfile);
    
    return ret;
}

std::vector<double> sigMOD::getSig(std::string sigfile) {
    std::vector<double> ret;
    
    ret = getVectorFromFile(sigfile);
    
    return ret;
}

std::map<int,double> sigMOD::getLogFactMap(std::string logfactfile){
  std::map<int,double> ret;

  ret = getMapFromFile(logfactfile);

  return ret;
}

std::vector<double> sigMOD::getVectorFromFile(std::string filename) {
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

std::map<int,double> sigMOD::getMapFromFile(std::string filename) {

  std::map<int,double> ret;
  
  std::ifstream myfile;
  myfile.open(filename);
  
  std::string line; 

  int temp1;
 long   double temp2;
  
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
