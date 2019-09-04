//
// Created by xos81802 on 21/10/2018.
//

#include <iostream>
#include <math.h>
#include "Histogram.h"


Histogram::Histogram (std::string filename) : filename(filename) {
    counts.resize(13);
    squared.resize(13);
    std::fill(counts.begin(), counts.end(), 0.0);
    std::fill(squared.begin(), squared.end(), 0.0);
}



void Histogram::addToRunningDistribution(const std::vector<double> &contactsDistribution){

    // normalize
    double total=0;
    for(unsigned int i=0; i< 13; i++){
        total += contactsDistribution[i];
    }
    double invTotal = 1.0/total;


    for(unsigned int i=0; i< 13; i++){
        double value = contactsDistribution[i]*invTotal;
        counts[i] += value;
        squared[i] += value*value;
    }

    counter+=1;
}



void Histogram::printResults(){

    double invCounter = 1.0/(double)counter;
    std::vector<double> averages(13);
    std::vector<double> variances(13);
    double sum = 0;

    for(unsigned int i=0; i< 13; i++){
        averages[i] = counts[i]*invCounter;
        sum+=averages[i];
        variances[i] = squared[i]*invCounter - averages[i]*averages[i];
        std::cout << "SQUARED AVE " << squared[i]*invCounter << " AVE-SQUARED " << averages[i]*averages[i] << std::endl;
    }

    for(unsigned int i=0; i< 13; i++){
        std::cout << i << " " << averages[i] << " - " << std::sqrt(variances[i]) << std::endl;
    }
    std::cout << " SUM : " << sum << std::endl;
}



