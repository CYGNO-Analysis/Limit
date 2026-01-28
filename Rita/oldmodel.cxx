#include "oldmodel.h"

oldmodel::oldmodel(const std::string& name)
    : BCModel(name)
{
    std::cout<<"Starting fit for '"<<name<<"' oldmodel"<<std::endl;
    
    // Loading input files
    data = getTable("./input/data.txt");
    bkg  = getTable("./input/background.txt");
    sig  = getTable("./input/background.txt");
    
    // Computing mubmax
    mubmax = std::accumulate(data.begin(), data.end(), 0.0);
    //mubmax = 5000000.;
    std::cout<<"DEBUG: mubmax = "<<mubmax<<std::endl;
    
    // [FIXME]: put this as an argument of the constructor
    bkg_mean = 130591;
    bkg_sigma = 2270;
    
    norm_bkg = std::accumulate(bkg.begin(), bkg.end(), 0.0);
    norm_sig = std::accumulate(sig.begin(), sig.end(), 0.0);
    
    // rescale bkg
    std::transform(bkg.begin(),bkg.end(), bkg.begin(), [&] (double a) {
        return a/norm_bkg;
    });
    
    // define bkg uncertainty in each bin
    //double bkg_percentage_uncertainty = 0.05;
    //std::transform(bkg.begin(),bkg.end(), std::back_inserter(bkg_err), [&] (double a) {
    //    return a * bkg_percentage_uncertainty;
    //});
    
    // rescale sig
    std::transform(sig.begin(),sig.end(), sig.begin(), [&] (double a) {
        return a/norm_sig;
    });

    
    // Add background events parameter
    AddParameter("mub", 0., 3.*mubmax, "mu_b", "[events]");
    
    // Add signal events parameter
    AddParameter("mus", 0., 3.*mubmax, "mu_s", "[events]");
    // Uncomment following line if you need a bkgonly fit [to do in a smarter way]
    //GetParameter("mus").Fix(0.0);
    
    // To parallelize the chains [TO BE CHECKED]
    //omp_set_dynamic(0);
    //omp_set_num_threads(nth);
}


double oldmodel::LogLikelihood(const std::vector<double>& pars) {
    double LL = 0.0;
    double mub = pars[0];
    double mus = pars[1];
    
    unsigned int N = data.size();
    
    for(unsigned int i =0; i<N; i++) {
        double lambda_i = mub * bkg[i] + mus * sig[i];
        //LL += BCMath::LogGaus((int)data[i], lambda_i);
        LL += BCMath::LogPoisson((int)data[i], lambda_i);
    }
    
    return LL;
}

double oldmodel::LogAPrioriProbability(const std::vector<double>& pars) {
    double LP = 0.0;
    double mub = pars[0];
    
    // prior on mub
    //LP += BCMath::LogPoisson(mub, bkg_mean);
    LP += BCMath::LogGaus(mub, bkg_mean, bkg_sigma);
    
    // prior on mus
    LP += 1./1000;
    return LP;
}


// Utility

std::vector<double> oldmodel::getTable(std::string fname){
    std::vector<double> ret;
    std::ifstream inFile(fname.c_str());
    
    std::string line;
    int vcounter = 0;
    while(std::getline(inFile, line)) {
        if (vcounter!=0) { //skip first row
            std::stringstream to_split(line.c_str());
            std::string element;
            
            int hcounter = 0;
            while(std::getline(to_split, element, '\t')) {
                if (hcounter!=0 and element != "") { //skip first column}
                    element.erase(remove_if(element.begin(), element.end(), isspace), element.end());
                    //std::cout<<"DEBUG '"<<element<<"'"<<std::endl<<std::flush;
                    ret.push_back(std::stod(element));
                }
                hcounter++;
            }
        }
        vcounter++;
    }
    return ret;
}
