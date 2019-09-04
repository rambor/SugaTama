//
// Created by xos81802 on 18/09/2018.
//

#include "Logger.h"

Logger::Logger(std::string filename) : outputfilename(filename), limit(0) {


}

Logger::Logger(std::string filename, unsigned int limit) : outputfilename(filename), limit(limit) {

    temperatures.resize(limit);
    klDivergences.resize(limit);
    workingLimits.resize(limit);

    pTemp = &temperatures.front();
    pDivergence = &klDivergences.front();
    pWorkingLimit = &workingLimits.front();

    counter = 0;
}

void Logger::logIt(std::string line) {

    outputlines.emplace_back(line);
}


void Logger::updateParams(float temp, float kl, unsigned wl){
    *(pTemp+counter) = temp;
    *(pDivergence+counter) = kl;
    *(pWorkingLimit+counter) = wl;
}


void Logger::writeToFile() {

    const char *outputFileName;
    std::string tempname = outputfilename + ".log";
    outputFileName = tempname.c_str() ;

    FILE * pFile;
    pFile = fopen(outputFileName, "w");

    unsigned int lines = outputlines.size();
    for (unsigned int i=0; i < lines; i++){
        std::fprintf(pFile, outputlines[i].c_str());
    }

}