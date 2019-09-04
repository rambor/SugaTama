//
// Created by xos81802 on 20/09/2018.
//

#include "CEConfiguration.h"

CEConfiguration::CEConfiguration(){};

CEConfiguration::CEConfiguration(float binWidth, unsigned int totalInPotential) {
    this->totalInPotential = totalInPotential;
    this->binWidth = binWidth;
    subUnitAssignments.resize(this->totalInPotential); // equal to the number of subunits - 1, hold integer referring to bin assigment
    assignedDistances.resize(this->totalInPotential);  // equal to the number of subunits - 1
}

void CEConfiguration::populateAssignment(unsigned int bin, unsigned int assignedBin) {
    subUnitAssignments[bin] = assignedBin;
    assignedDistances[bin] = 0.5*binWidth + bin*binWidth;
}