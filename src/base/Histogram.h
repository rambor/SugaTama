//
// Created by xos81802 on 21/10/2018.
//

#ifndef SUGATAMA_HISTOGRAM_H
#define SUGATAMA_HISTOGRAM_H

#include <string>
#include <vector>

class Histogram {

    std::vector<double> counts;
    std::vector<double> squared;
    const std::string filename;
    unsigned int counter=0;

public:
    Histogram(std::string name);

    void addToRunningDistribution(const std::vector<double> &contactsDistribution);
    void printResults();
};


#endif //SUGATAMA_HISTOGRAM_H
