# Samuele's code
Neutrino feasibility study described in his [thesis](https://www.arxiv.org/pdf/2408.03760) (Chapter 7). Estimation with direction and energy dependence of a CYGNO-like 30 cubic meter detector of the sensitivity to neutrinos 

## Code flow
1. Retrieve from CYGNO simulation the files referring to internal and or external background (tipically ROOT files with nTuple) 
2. Compile and run the background template codes (high statistics required)
3. Compile and run the signal spectrum template codes (high statistics required)
4. Starting from templates, generate toy MC fake dataset (possibly code generating toys not present in this repo)
5. Run the Bayesian analysis code

## Main code content
### spectrum.txt
Contains the Solar neutrino energy spectrum (they need to be multiplied by 10<sup>11</sup>)

### TemplateBackgroundED.cpp
This creates the background template starting from the Geant4 ROOT simulation files 

### TemplateSignal.cpp
This creates signal template starting from spectrum.txt. Beware: it is important to choose the energy/angular resolution as this determines the size of the bins and the data smearing

### cygnoNeutrino/bkgModel.h 
This will contain the background model

### cygnoNeutrino/sigModel.h 
This will contain the signal model

### cygnoNeutrino/runAll.sh 
This shows how to correctly launch the code for the Bayesian analysis. 

## Final Analysis flow
- Compile inside cygnoNeutrino folder with ```make```
- Arguments for the executable: toy.txt containing the toyMC data, background and signal (to be converted in xt files with code possibly present in the folder)
- Example on how to run

``` ./runfit ../Out_V8_CosThetaFlat_HSBkg_NID/toyMC_txt/${i}/${j}.txt ../Out_V8_CosThetaFlat_HSBkg_NID/template_txt/background.txt ../Out_V8_CosThetaFlat_HSBkg_NID/template_txt/signal.txt ${j} ${i} NID```
