//
// Created by xos81802 on 20/09/2018.
//

#ifndef SUGATAMA_CECONFIGURATION_H
#define SUGATAMA_CECONFIGURATION_H



#include <vector>

class CEConfiguration {
    std::vector<unsigned int> subUnitAssignments; // equal to the number of subunits - 1, hold integer referring to bin assigment
    std::vector<float> assignedDistances;  // equal to the number of subunits - 1
    unsigned int totalInPotential;
    float binWidth;
    float score;

public:
    CEConfiguration();
    CEConfiguration(float binWidth, unsigned int totalInPotential);

    void populateAssignment(unsigned int bin, unsigned int assignedBin);
    unsigned int getBinAssignment(unsigned int vectorIndex){ return subUnitAssignments[vectorIndex];}

    float getScore() const { return score;}
    void setScore(float score){ this->score = score;}
};


#endif //SUGATAMA_CECONFIGURATION_H
