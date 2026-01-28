# Samuele's code
Neutrino feasibility study described in his [thesis](https://www.arxiv.org/pdf/2408.03760) (Chapter 7). Estimation with direction and energy dependence of a CYGNO-like 30 cubic meter detector of the sensitivity to neutrinos 

## Code flow
1. Retrieve from CYGNO simulation the files referring to internal and or external background (tipically ROOT files with nTuple) 
2. Compile and run the background template codes (high statistics required)
3. Compile and run the signal spectrum template codes (high statistics required)
4. Starting from templates, generate toy MC fake dataset (possibly code generating toys not present in this repo)
5. Run the Bayesian analysis code

## Content
### spectrum.txt
Contains the Solar neutrino energy spectrum (they need to be multiplied by 10<sup>11</sup>






-spectrum.txt

contiene il flusso di neutrini solari in funzione dell'energia, credo ci sia un fattore 10^11 da moltiplicare

-TemplateBackgroundED.cpp

crea il template del background a partire dagli spettri simulati da CYGNO30

-TemplateSignal.cpp

crea i template di segnale a partire dallo spettro dei neutrini solari considerando l'interazione N.B. per queste due cose è importante la risoluzione che determina la dimensione dei bin (1sigma) e fa lo smearing


-Cygno neutrino
|
|-bkgModel.h ha il modello di background
|
|-sigModel.h ha il modello del segnale
|
|-runAll.sh è interessante xk ti spiega come lanciare il codice per fare l'analisi bayesiana (argomenti)
|



5) Girare l'analisi bayesiana che si compila con make, gli argomenti sono il toy.txt template background (che va convertito in txt, ci dovrebbero essere i codici), template segnale in txt (dettagli da runAll.sh)
./runfit ../Out_V8_CosThetaFlat_HSBkg_NID/toyMC_txt/${i}/${j}.txt ../Out_V8_CosThetaFlat_HSBkg_NID/template_txt/background.txt ../Out_V8_CosThetaFlat_HSBkg_NID/template_txt/signal.txt ${j} ${i} NID
