#ifndef old_model_h
#define old_model_h

#include <BAT/BCModel.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
//#include <omp.h> //for parallelization
#include <numeric>


#include <BAT/BCMath.h>
#include "Math/ProbFunc.h"

class oldmodel : public BCModel
{
public:

    // Constructor
    oldmodel(const std::string& name);
    // Destructor
    ~oldmodel(){};
    
    double LogLikelihood(const std::vector<double>& pars);
    double LogAPrioriProbability(const std::vector<double>& pars);
    
    // Utility
    std::vector<double> getTable(std::string fname);
    
private:
    
    std::vector<double> data;
    std::vector<double> bkg;//, bkg_err;
    std::vector<double> sig;
    
    double mubmax   = 0.0;
    double bkg_mean = 0.0;
    double bkg_sigma = 0.0;
    double norm_bkg = 0.0;
    double norm_sig = 0.0;

};


#endif /* old_model_h */
