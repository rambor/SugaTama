//
// Created by xos81802 on 18/09/2018.
//

#ifndef SUGATAMA_LOGGER_H
#define SUGATAMA_LOGGER_H

#include <string>
#include <vector>
#include <iostream>

class Logger {

    const unsigned int limit;
    unsigned int counter;

    std::vector<float> temperatures;
    std::vector<float> klDivergences;
    std::vector<unsigned int> workingLimits;

    std::string outputfilename;
    std::string datafilename;
    std::vector<std::string> outputlines;

    float * pTemp;// = &tempDuringRun.front();
    float * pDivergence;// = &divergenceDuringRun.front();
    unsigned int * pWorkingLimit;// = &workingLimitDuringRun.front();

public:
    Logger(std::string outputfilename);

    Logger(std::string outputfilename, unsigned int limit);

    void logIt(std::string line);

    void writeToFile();

    void updateParams(float temp, float kl, unsigned wl);
    void setDataFilename(std::string name){ this->datafilename = name;}
};


#endif //SUGATAMA_LOGGER_H
