#include <BAT/BCLog.h>

#include <iostream>
#include <string>
#include <stdexcept>

#include "oldmodel.h"


int main(int argc, char *argv[]) {
    
    int Nch, NIter;
    Nch = 8;
    NIter = 10000000;
    
    // Choosing the model
    std::string title    = "test";
    oldmodel m(title);
    
    //throw std::runtime_error("DEBUG");
    
    std::string res_dir = "./results/test/";
    
    // Creating the result folder if it does not exist
    int com = std::system(("mkdir "+res_dir).c_str());
    if(com != 0) {
        std::cout<<com<<std::endl;
    }
    
    BCLog::OpenLog(res_dir + title + "_log.txt", BCLog::detail, BCLog::detail);

    m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    m.SetPrecision(BCEngineMCMC::kMedium);
    
    // Set number of iterations of prerun phase (ideally infinite == until convergence)
    m.SetNIterationsPreRunMax(100000000);

    // Set the number of the MCMC integration iterations and chains
    m.SetNIterationsRun(NIter);
    m.SetNChains(Nch);
    
    
    // Set initial point (not needed, but it can make the convergence faster)
    // std::vector<double> x0;
    // x0.push_pack(0.0);
    // m.SetInitialPositions(x0);
    
    m.WriteMarkovChain(res_dir+m.GetSafeName() + "_mcmc.root", "RECREATE");//, true, true);
    
    // Run MCMC, marginalizing posterior
    m.MarginalizeAll();
    // Definisci le variabili per i limiti
double lower_bound_mu, upper_bound_mu;

// Ottieni il quantile inferiore (5% per un intervallo di confidenza al 90%)
lower_bound_mu = m.GetMarginalized(0).GetQuantile(0.05);

// Ottieni il quantile superiore (95% per un intervallo di confidenza al 90%)
upper_bound_mu = m.GetMarginalized(0).GetQuantile(0.95);

// Stampa i limiti
std::cout << "Intervallo di confidenza al 90% per mu: ["
          << lower_bound_mu << ", " << upper_bound_mu << "]" << std::endl;



    // Print summary plots
    m.PrintParameterPlot(res_dir+m.GetSafeName() + "_parameters.pdf");
    m.PrintCorrelationPlot(res_dir+m.GetSafeName() + "_correlation.pdf");
    m.PrintCorrelationMatrix(res_dir+m.GetSafeName() + "_correlationMatrix.pdf");
    m.PrintKnowledgeUpdatePlots(res_dir+m.GetSafeName() + "_update.pdf");

    // Print results of the analysis into a text file
    m.PrintSummary();
    
    // Close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
    
}

