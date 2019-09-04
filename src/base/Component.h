//
// Created by xos81802 on 14/07/2018.
//

#ifndef SUGATAMA_COMPONENT_H
#define SUGATAMA_COMPONENT_H

#include <string>
#include <vector>
#include <math.h>
#include <regex>
#include <iostream>
#include <set>
#include "../EulerTour/EulerTour.h"

#ifdef __cplusplus
extern "C" {
#endif
#include <libqhull/qhull_a.h>
#ifdef __cplusplus
}
#endif

class Model;
class Component {

private:

    std::string id;
    std::set<unsigned int> anchors;            // can be empty
    std::set<unsigned int> centeredAnchors;            // can be empty, anchor is bead that is centered on the residued

    std::vector<unsigned int> resids;          // can be empty
    std::vector<std::string> chains;  // can be empty
    std::vector<float> potentialFunction;  // can be empty

    std::set<unsigned int> beads_in_use;       // never empty
    std::set<unsigned int> best;       // never empty

    std::vector<unsigned int> cvxPoints;       // never empty
    std::vector<double> contactsDistributionOfModel;
    EulerTour tour;

    std::set<unsigned int> hull;

    unsigned int totalComponents;
    unsigned int targetNumberOfBeads=0;
    float targetVolume, invTargetVolume;
    Model * pModel;  // for the euler tour
    float currentCVXVolume;
    bool empty = true;
    double totalContactSum;
    unsigned int anchorCount=0;
    float percentageStep;

    double contactsKLDivergence;

public:

    Component(std::string id, float volume, Model *pModel);
    Component(std::string id, float volume, Model *pModel, bool contiguous);


    bool checkID(std::string str);
    float numberOfLatticePointsPotential(float value);

    unsigned int removeLatticePoint(unsigned int index);
    unsigned int addLatticePoint(unsigned int index); // add to beads_in_use and

    bool indexInUse(unsigned int index);

    void addResid(unsigned int index, std::string chain);

    std::string const getID() const { return id;}
    unsigned int getTotalResids() { return resids.size();}
    unsigned int getResidByIndex(unsigned int index) { return resids[index];}
    std::string getChainsByIndex(unsigned int index) { return chains[index];}
    void addAnchor(unsigned int index);

    // desired number of beads but exclude the anchors
    void setTargetNumberOfLatticePoints(float value);
    void addCenteredAnchors(unsigned int index);
    unsigned int getTargetNumberOfLatticePoints();

    bool anchorsEmpty(){ return empty;}

    std::set<unsigned int> * getAnchors(){ return &anchors; }
    std::set<unsigned int> * getCenteredAnchors(){ return &centeredAnchors; }
    std::set<unsigned int> * getBeadsInUse() { return &beads_in_use; }
    std::set<unsigned int> * getHull(){ return &hull;}
    std::vector<double> * getContactsDistribution(){ return &contactsDistributionOfModel;}

    unsigned int getTotalNumberOfBeads() { return beads_in_use.size(); }

    unsigned int getTotalNumberOfComponents(){ return tour.getNumberOfComponents();}

    bool inUse(unsigned int index){ return !(beads_in_use.find(index) == beads_in_use.end()); }
    bool isAnchor(unsigned int index);
    bool isCenteredAnchor(unsigned int index);
    float potential();
    void printAnchors();
    void printCenteredAnchors();
    void printContactsPerCenteredAnchor();

    void printSet();
    void printBest();
    void writeToFile(std::string nameOf);
    void writeAnchorsToFile(std::string nameOf);
    void writeCenteredAnchorsToFile();
    void populatePotential(float percentage);

    void copyBestToInUse();
    float calculateCVXVolume();
    void setBest();

    void setTotalContactSum(double value){ totalContactSum = value;}
    double getTotalContactSum(){ return totalContactSum;}
    void printConstraints();
    unsigned int getAnchorCount(){return anchorCount;}
    void setCurrentCVXVolume();
    void setCurrentCVXVolume(float value){ this->currentCVXVolume=value;}
    float getCurrentCVXVolume(){ return currentCVXVolume;}

    //float getCurrentVolumeEnergy(){ return std::abs(currentCVXVolume-targetVolume)*invTargetVolume;}
    float calculateVolumeEnergy(float testVolume){ return std::abs(testVolume-targetVolume)*invTargetVolume;}


    double setContactsKLDivergence(double val){ return contactsKLDivergence = val;}
    double getContactsKLDivergence(){ return contactsKLDivergence;}


};

#endif //SUGATAMA_COMPONENT_H
