#!/bin/bash

for i in $(ls ../Out_V8_CosThetaFlat_HSBkg/toyMC_txt | grep "$1\b")
do
    for j in {1..10000}   
    do
        ./runfit ../Out_V8_CosThetaFlat_HSBkg/toyMC_txt/${i}/${j}.txt ../Out_V8_CosThetaFlat_HSBkg/template_txt/background.txt ../Out_V8_CosThetaFlat_HSBkg/template_txt/signal.txt ${j} ${i} ED
    done
done
