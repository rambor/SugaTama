//
// Created by xos81802 on 10/07/2018.
//

#ifndef SUGATAMA_ANNEAL_H
#define SUGATAMA_ANNEAL_H

#include "thirdparty/Eigen/Core"
#include "thirdparty/Eigen/Geometry"
#include "thirdparty/power_sasa.h"

#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <Data.h>
#include <random>
#include <Component.h>
#include <SubUnit.h>
#include "Model.h"
#include "Aligner/KDE.h"
#include <IntraSubUnitContact.h>
#include <CEPotential.h>

#ifdef __cplusplus
extern "C" {
#endif
#include <libqhull/qhull_a.h>
#ifdef __cplusplus
}
#endif

class Data;

//struct Trial{ double value; std::vector<unsigned int> indices;};

class Anneal {

    struct Trial{ double value; std::vector<unsigned int> indices;
        Trial(double val, std::vector<unsigned int>  & vec) : value(val) {
            //indices = std::vector<unsigned int>(vec);
            indices = std::move(vec);
        }

        bool operator<(const Trial & a) const{
            return value < a.value;
        }
    };

    struct find_pt_by_key : std::unary_function< KDE::ProbabilityBead, bool>{
        find_pt_by_key(unsigned int keyToFind) : key(keyToFind){}
        bool operator () (KDE::ProbabilityBead p) { return p.index == key; }
    private:
        unsigned int key;
    };

    float percentAddRemove, beta, eta, lambda, alpha, mu, stov, asaAcceptanceRate, complementASAAcceptanceRate, intASAAcceptanceRate, intComplementASAAcceptanceRate;
    float highT, highTempStartForCooling, interconnectivityCutOff;
    unsigned int lowerV, upperV, highTempRounds, ccmultiple;
    bool isRefine=false;
    unsigned int totalNumberOfPhasesForSeededModeling=1;
    std::string filenameprefix;
    std::vector<double> contactsDistribution;

private:

    enum TLogLevel {logDETAIL, logWARNING, logINFO, logDEBUG};

    std::vector<Component> components;
    unsigned int totalComponents=0, distributionlimit=13;
    float contactCutOff, violation_limit, probe_radius, delta_r = 1.4f;
    unsigned int maxbin, totalBins;

    std::vector<Eigen::Vector3f> totalCoordinatesInObject; // for symmetry related object
    std::vector<IntraSubUnitContact> intraSubUnitContacts;

    vector3 translationVector;
    std::vector<CEPotential> cePotentials;

    float surfaceToVolume(const unsigned int total, std::vector<float> & weights, std::vector<Eigen::Vector3f> & coordinates);

    float changeInSurfaceToVolumeAdd( unsigned int beadToAdd, std::set<unsigned int> *beads_in_use, Model * pModel );

    void addToContactsDistribution(unsigned int addMe, std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, Model *pModel);
    void removeFromContactsDistribution(unsigned int removeMe, std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, Model *pModel);

    void populateContactsDistribution(std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, Model *pModel);
    void calculateModelPrDistributionDirect(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, const unsigned int workingLimit, Model *pModel, Data *pData);

    double connectivityPotential(unsigned int numberOfComponents);
    double contactPotential(std::set<unsigned int> *beads_in_use, Model *pModel);

    void updateASATemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp);
    void updateASAConstantTemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp);
    void updateASAModTemp(unsigned int index, float evalMax, float acceptRate, double &temp, double &inv_temp);

    void addForSurfaceToVolume(unsigned int toAdd, std::vector<unsigned int> & unsortedIndices, std::vector<float> & weights, std::vector<Eigen::Vector3f> & coordinates, Model * pModel);
    float surfaceAreaSymObject(const unsigned int workingLimit, const unsigned int totalSubUnits, float sa);

    void removeForSurfaceToVolume(unsigned int toRemove, unsigned int workingLimit, std::vector<unsigned int> & unsortedIndices, std::vector<float> & weights, std::vector<Eigen::Vector3f> & coordinates);
    void swapForSurfaceToVolume(unsigned int indexToFind, unsigned int indexToAdd, std::vector<unsigned int> & unsortedIndices, std::vector<Eigen::Vector3f> & coordinates, Model * pModel);
    void printParameters(std::vector<float> * accept, std::vector<double> * temp, std::vector<float> * divergence, std::vector<unsigned int> * wl);

    float calculateKLDivergenceAgainstPDBPR(std::vector<unsigned int> &modelPR, std::vector<double> &targetPR);

    void recalculateDeadLimit(unsigned int workingLimit, std::vector<unsigned int> &bead_indices, std::set<unsigned int> &hull, Model * pModel, unsigned int totalBeadsInSphere);

    void getHullPoints(std::set<unsigned int> &hullpts, std::set<unsigned int> &beads_in_use, Model * pModel);
    void populateLayeredDeadlimitUsingSet(std::set<unsigned int> & beads_in_use, std::set<unsigned int> & hull, Model * pModel);

    vector3 rotate(vector3 const & vec, double cosx, double sinx, double cosz, double sinz);

/**
 *
 * calculate Kulback Liebler divergence of contacts distribution
 * @param distribution
 * @return
 */
    inline double calculateKLDivergenceContactsDistribution(std::vector<double> &distribution) {

        double kl=0.0d, totalCounts = 0.0;
        for(unsigned int i=0; i< 13; i++){
            totalCounts += distribution[i];
        }


        for (unsigned int i=0; i < distributionlimit; i++){
            // i know every value in working_probability up to zeroBin is nonzero
            double prob = contactsDistribution[i];  // bounded by experimental Shannon Number
            double tempvalue = distribution[i];
            if (prob > 0){
                if (tempvalue > 0){
                    kl += prob * std::log(prob / tempvalue * totalCounts);
                } else { // if tempvalue is Zero (empty)
                    kl += 1.1;
                }
            }
        }

        return kl;
    }

    unsigned int removeFromPrSym(unsigned int const &removeMeSubUnitIndex, std::vector<unsigned int> &beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> &prBins,
                        Model *pModel, Data *pData);

    unsigned int removeFromPrSymHelical(unsigned int const &removeMeSubUnitIndex, std::vector<unsigned int> &beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> &prBins,
                                 Model *pModel, Data *pData);

    unsigned int addToPrSym(const unsigned int addMeSubUnitIndex, std::vector<unsigned int> &beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> &prBins,
                   Model *pModel, Data *pData);

    unsigned int addToPrSymHelical(const unsigned int addMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData);

    bool checkForRepeats(std::vector<unsigned int> beads) ;
    bool checkSetAndVector(unsigned int workingLimit, std::vector<unsigned int> * indices, std::set<unsigned int> * beads_in_use);

    double calculateKLEnergySymmetry(std::vector<unsigned int> *subUnit_indices, std::vector<unsigned int> *binCount, const unsigned int indicesWorkingLimit, unsigned int &violation, Model *pModel, Data *pData);
    float calculateCVXVolumeSymmetry(std::vector<unsigned int> *subUnit_indices, const unsigned int indicesWorkingLimit, Model *pModel);

    void updateContactsDistribution(std::set<unsigned int> *beads_in_use_tree, Model *pModel);
    void populateContactsDistributionSym(std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, Model *pModel);



    void addToPrDirect(const unsigned int addMe, std::vector<unsigned int> & beadsInUse, const unsigned int upperLimit, std::vector<unsigned int> & prBins, Model * pModel, Data * pData );

    void removeLatticePositionToModelDirect(
            std::vector<unsigned int> & bead_indices,
            std::vector<unsigned int> & pBinCount,
            unsigned int * pWorkingLimit,
            const unsigned int * pLatticePointToRemove, Model * pModel, Data * pData);

    void removeFromPrDirect(unsigned int removeMe, std::vector<unsigned int> & beadsInUse, const unsigned int upperLimit, std::vector<unsigned int> & prBins, Model * pModel, Data * pData);

    void restoreRemovingLatticePointFromBackUpDirect(
            unsigned int * pWorkingLimit,
            std::vector<unsigned int> * pBinCountBackUp,
            std::vector<unsigned int> * pBinCount);

    void restoreAddingFromBackUpDirect(
                                                      unsigned int * pWorkingLimit,
                                                      std::vector<unsigned int> * pBinCountBackUp,
                                                      std::vector<unsigned int> * pBinCount);

    void fillPrBinsAndAssignTotalBinDirect(Model * pModel, Data * pData);

    void addLatticPositionToModelDirect(std::vector<unsigned int> * pIndices,
                                        unsigned int * pWorkingLimit,
                                        std::vector<unsigned int>::iterator * pItIndex);


    void calculateModelPrDistributionSymX(std::vector<unsigned int> *subUnit_indices, std::vector<SubUnit> *subUnits, std::vector<unsigned int> *binCount, const unsigned int indicesWorkingLimit, double &potential, Model *pModel, Data *pData);

    double intraSubUnitContact(std::vector<SubUnit> * subUnits, std::vector<unsigned int> * cvxIndices, Model *pModel);

    double removeFromPrSymX(unsigned int const &removeMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, std::vector<SubUnit> & subUnits, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData);
    void swapInSubUnits(unsigned int indexOfCurrent, unsigned int indexOfNew, std::vector<SubUnit> & subUnits, Model * pModel);
    void removeFromSubUnits(unsigned int indexOfPositionToRemove, std::vector<SubUnit> & subUnits);
    double addToPrSymX(std::vector<unsigned int> & beadsInUse, std::vector<SubUnit> *subUnits, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData);
    double removeLatticePositionToModelSymX(
            std::vector<unsigned int> & bead_indices,
            std::vector<SubUnit> & subUnits,
            std::vector<unsigned int> & modelPrBins,
            unsigned int * pWorkingLimit,
            const unsigned int * pLatticePointToRemove, Model * pModel, Data *pData);

    void undoRemoveFromSubUnits(std::vector<SubUnit> &subUnits);

    void addToSubUnits(unsigned int indexOfPositionToAdd, std::vector<SubUnit> & subUnits, Model * pModel);
    void undoAddToSubUnits(std::vector<SubUnit> &subUnits);
    float symXRotateTranslate(unsigned int round, unsigned int index, float step_limit, std::vector<SubUnit> & subunits, std::vector<CEPotential> & cePotentials, Data * pData,  Model * pModel, bool printIt, std::string name);
    double moments(std::vector<SubUnit> *subUnits, std::vector<unsigned int> *indices, unsigned int totalIndices, float targetRg, Model *pModel);
    float calculateCOMPotential(const vector3 * const com, std::vector<SubUnit> & subunits, std::vector<CEPotential> & cePotentials);

    void logToCeFile(char * text);


    void logger(std::string description, std::string value) {

        unsigned int len = 40 - description.size();
        std::string firsthalf = std::string(len, ' ');
        printf("%s%s : %s\n", firsthalf.c_str(), description.c_str(), value.c_str());

    }


    /**
     * Formats Number to String for standard output
     * decimals greater than 7 defaults to scientific notation
     *
     * @param number
     * @param decimals
     * @return
     */
    std::string formatNumber(float number, int decimals = 2) {
        char buffer [50];
        switch(decimals){
            case 1 :
                sprintf (buffer, "%.1f", number);
                break;
            case 2 :
                sprintf (buffer, "%.2f", number);
                break;
            case 3 :
                sprintf (buffer, "%.3f", number);
                break;
            case 4 :
                sprintf (buffer, "%.4f", number);
                break;
            case 5 :
                sprintf (buffer, "%.5f", number);
                break;
            case 6 :
                sprintf (buffer, "%.6f", number);
                break;
            case 7 :
                sprintf (buffer, "%.7f", number);
                break;
            default:
                sprintf (buffer, "%.2E", number);
        }

        return std::string(buffer);
    }

public:

    Anneal(float highT,
           float percent,
           unsigned int highTempRounds,
           std::string prefix,
           float alpha,
           float beta,
           float eta,
           float lambda,
           float mu,
           unsigned int multiple,
           float accRate,
           float interconn
    );


    float calculateCVXHULLVolume(char * flags, std::vector<unsigned int> *bead_indices, unsigned int upTo, Model *pModel);
    float calculateCVXHULLVolumeSet(char *flags, std::set<unsigned int> *bead_indices_tree, unsigned int upTo, Model *pModel);

    void estimateInterSubUnitDistancesX(Model *pModel, Data *pData);
    float estimateInterSubUnitEnergyX(std::vector<SubUnit> & subunits, Model *pModel, Data *pData);

    bool canRemoveIfNotCenteredAnchor(unsigned int index);

    bool createInitialModelCVXHull(Model *pModel, Data *pData, std::string name);
    bool createInitialModelCVXHullDirect(Model *pModel, Data *pData, std::string name);
    bool createInitialModelSymmetry(Model *pModel, Data *pData);
    bool createInitialModelSymmetryEx(Model *pModel, Data *pData);

    bool createInitialModelHelicalSymmetry(Model *pModel, Data *pData);
    bool createInitialModelSymmetryX(Model *pModel, Data *pData);
    bool createInitialModelCVXHullSeeded(Model *pModel, Data *pData, std::string name);
    bool createSeedFromPDB(Model *pModel, Data *pData, std::string name, std::string PDBFilename, unsigned int totalPhases);
    bool createSeedFromPDBDirect(Model *pModel, Data *pData, std::string name, std::string PDBFilename, unsigned int numberOfUniqueConnectedPhases);


    float connectivityPotentialPhases(unsigned int mainConnectivity);
    float biCameralVolumePotential(unsigned int lowerBound, unsigned int upperBound, float value);

    std::string refineSymModel(Model *pModel, Data *pData, std::string nameTo);
    std::string refineSymModelHelical(Model *pModel, Data *pData, std::string nameTo);
    std::string refineHomogenousBodyASAHybridEx(Model *pModel, Data *pData, std::string outputname);
    std::string refineHomogenousBodyASAHybridDirect(Model *pModel, Data *pData, std::string outputname);
    std::string refineHomogenousBodyASACVXSeeded(Model *pModel, Data *pData, std::string outputname);
    std::string refineSymModelRefine(Model *pModel, Data *pData, std::string nameTo);
    std::string refineHomogenousBodyInputEx(Model *pModel, Data *pData, std::string outputname);
    std::string refineHomogenousBodyASAHybridSymmetryX(Model *pModel, Data *pData, std::string nameTo);

    void ceMapOptimization(Model *pModel, Data *pData, unsigned int topN, unsigned int minN, unsigned int maxN, std::vector<KDE::ProbabilityBead> & lattice);
    void ceMapOptimizationSym(Model *pModel, Data *pData, unsigned int topN, const unsigned int minN, const unsigned int maxN, std::vector<KDE::ProbabilityBead> & lattice);

    bool initializeModelToRefineSym(Model *pModel, Data *pData, std::string name, std::string PDBFilename);
    // bool setAnchorPoints(std::string anchorFileName, std::string pdbFile, Model *pModel);

    bool initializeModelToRefine(Model *pModel, Data *pData, std::string name, std::string PDBFilename);
    bool setAnchorPoints(std::string anchorFileName, std::string pdbFile, Model *pModel);
    void setContactsDistribution(std::vector<unsigned int> & contactsDistributionSeed);

    unsigned int getMaxBinUsedInAnnealing(){ return maxbin;}
    unsigned int gettotalBinsDerivedFromData(){ return totalBins;}
    float getHighTempStartForCooling(){return highTempStartForCooling;}
    float getPercentAddRemove(){ return percentAddRemove;}
    float getAlpha(){ return alpha;}
    float getLambda(){ return lambda;}
    float getMu(){ return mu;}
    float getEta(){ return eta;}
    float getBeta(){ return beta;}
    float getStoV(){return stov;}

    double contactPotentialFunction(unsigned int nc);
    double addToContactsPotential(unsigned int addMe, double currentContactsSum, std::set<unsigned int> *beads_in_use, Model *pModel);
    double removeFromContactsPotential(unsigned int removeMe, double currentContactsSum, std::set<unsigned int> *beads_in_use, Model *pModel);

    const std::vector<double> &getContactsDistribution() const;
    void writeContactsDistributionToFile(float bin_width, float bead_radius, float qmax);
    void updateContactsDistributionToFile(std::string filename);

    inline void beadToPoint(pointT *ptestPoint, Bead *pBead) {
        ptestPoint[0] = pBead->getX();
        ptestPoint[1] = pBead->getY();
        ptestPoint[2] = pBead->getZ();
    }

    void calculateModelPrDistributionSym(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, unsigned int beadIndiciesWorkingLimit, unsigned int &violation, Model *pModel, Data *pData);
    void calculateModelPrDistributionSymCE(std::vector<unsigned int> *bead_indices, std::set<unsigned int> & subUnit_indices, std::map<unsigned int, std::vector<vector3> > & map , std::vector<double> & contactsDistributionOfModel, std::vector<unsigned int> *binCount, unsigned int beadIndiciesWorkingLimit, unsigned int &violation, Model *pModel, Data *pData);
    void calculateModelPrDistributionSymHelical(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, unsigned int beadIndiciesWorkingLimit, unsigned int &violation, Model *pModel, Data *pData);
    void calculateModelPrDistributionSymEX(std::vector<vector3> & subUnit_indices, std::vector<unsigned int> *binCount, std::vector<double> & contactsDistributionOfModel, const unsigned int indicesWorkingLimit, Model *pModel, Data *pData, unsigned int & violations);

    void calculateModelPrDistribution(std::vector<unsigned int> *bead_indices, std::vector<unsigned int> *binCount, unsigned int workingLimit, unsigned int totalBeadsInSphere, Model *pModel, Data *pData) ;
    void calculateModelParametersSymmetry(std::set<unsigned int> *subUnit_indices, Model *pModel);
    unsigned int getViolations(std::vector<unsigned int> *bead_indices, unsigned int beadIndiciesWorkingLimit, Model *pModel);
    unsigned int getViolationsFromSet(std::set<unsigned int> *subUnit_indices, std::vector<double> & contactsDistributionOfModel, unsigned int indicesWorkingLimit, Model *pModel);
    unsigned int getViolationsHelical(std::vector<unsigned int> *bead_indices, unsigned int beadIndiciesWorkingLimit, Model *pModel, Data *pData);
    unsigned int minimizeViolations(std::vector<unsigned int> *bead_indices, unsigned int beadIndiciesWorkingLimit, Model *pModel);
    void printSymModel(std::vector<unsigned int> *subUnit_indices, const unsigned int indicesWorkingLimit, Model *pModel, Data *pData);

    unsigned int getUseableNeighborFromSet(std::set<unsigned int> *beads_in_use, Model *pModel, unsigned int & selectedIndex);
    unsigned int getConnectedNeighborFromSet(std::set<unsigned int> *beads_in_use, Model *pModel, unsigned int & selectedIndex);

    unsigned int getUseableNeighborFromRestrictedSet(std::set<unsigned int> *beads_in_use, Model *pModel, unsigned int & selectedIndex);

    unsigned int numberOfContactsFromSet(std::set<unsigned int> *beads_in_use, Model *pModel, unsigned int selectedIndex);
    unsigned int numberOfContactsFromSetSym(std::set<unsigned int> *beads_in_use, Model *pModel, unsigned int selectedIndex);

    float calculateRgSym(std::vector<unsigned int> &beadsInUse, unsigned int const &workingLimit, Model *pModel);

    void fillPrBinsAndAssignTotalBin(Model * pModel, Data * pData);

    unsigned int getUseableNeighborForCOM(Model *pModel, unsigned int & com);
    unsigned int getTranslateToPointNearCOM(Model *pModel, unsigned int & com);
    /**
     *
     */
    inline void addLatticPositionToModel(std::vector<unsigned int> * pIndices,
                                                 unsigned int * pWorkingLimit,
                                                 std::vector<unsigned int>::iterator * pItIndex){


        //std::copy(pIndices->begin(), pIndices->end(), pBackUpState->begin()); // make backup
        // make the swap at the border (workingLimit)
        std::iter_swap(pIndices->begin() + *pWorkingLimit, *pItIndex); // this swaps to working position and changes backup
        // increment workingLimit to include new position
        *pWorkingLimit += 1;
        std::sort(pIndices->begin(), pIndices->begin() + *pWorkingLimit);
    }


    /**
     * addMe must be present in beadsInUse
     * pBin is a pointer to the vector of binned distances in Model
     * prBin is the vector holding counts per bin
     * upperLimit is usually the working Limit of the set
     */
    inline void addToPr(const unsigned int addMe, const std::vector<unsigned int> & beadsInUse, const unsigned int upperLimit, const unsigned short int * const pBin, const unsigned int & totalBeadsInUniverse, std::vector<unsigned int> & prBins ){

        // beadsInUse must be sorted
        unsigned int row;
        unsigned int row2;


        unsigned int i=0;
        // add down column
        const unsigned int * const ptr = beadsInUse.data();

        while(ptr[i] < addMe && i < upperLimit){
            row = ptr[i];
            row2 = (row*totalBeadsInUniverse - (row*(row+1)/2)) + (addMe - row);// to prevent underFlow subtract 1 last
            prBins[ *(pBin + row2 - 1 ) ]++;
            i++;
        }
        // out of first loop, i = addMe
        i++;
        // Add across row
        row2 = addMe*(unsigned int)totalBeadsInUniverse - addMe*(addMe+1)/2 - addMe;
        while (i < upperLimit){
            prBins[ *(pBin + (row2 + ptr[i]) - 1 ) ]++;
            i++;
        }
    }


    inline void restoreAddingFromBackUp(std::vector<unsigned int> * pIndices,
                                                const std::vector<unsigned int> * pBackUpState,
                                                unsigned int * pWorkingLimit,
                                                const std::vector<unsigned int> * pBinCountBackUp,
                                                std::vector<unsigned int>::iterator * pBinCountBegin){

        std::copy(pBackUpState->begin(), pBackUpState->end(), pIndices->begin());
        *pWorkingLimit -= 1;
        std::copy(pBinCountBackUp->begin(), pBinCountBackUp->end(), *pBinCountBegin); //copy to bin count
    }


    /**
     *
     */
    inline void removeLatticePositionToModelByIndex(
            std::vector<unsigned int> & bead_indices,
            std::vector<unsigned int> & pBinCount,
            unsigned short int * const pBin,
            unsigned int * pWorkingLimit,
            unsigned int totalBeadsInSphere,
            const unsigned int latticeIndexToRemove){

        const auto itIndex = bead_indices.begin() + latticeIndexToRemove;
        // remove original from P(r)
        // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
        removeFromPr(*itIndex, bead_indices, *pWorkingLimit, pBin, totalBeadsInSphere, pBinCount);
        // reduce the workingLimit
        // if wl = 10
        // 0 1 2 3 4 5 6 7 8 9 10
        // remove 4, wl-=1
        // 0 1 2 3 9 5 6 7 8 4 10
        // sort to 9
        // 0 1 2 3 5 6 7 8 9 4 10
        //
        *pWorkingLimit -= 1;
        // swap selected point to the workingLimit
        std::iter_swap(itIndex, bead_indices.begin() + *pWorkingLimit);
        // still need to sort, swap changes the order
        // to save some time, sort should only be from swapped point to workingLimit
        std::sort(bead_indices.begin(), bead_indices.begin() + *pWorkingLimit);
    }

/**
 *
 */
    inline void removeLatticePositionToModel(
            std::vector<unsigned int> & bead_indices,
            std::vector<unsigned int> & pBinCount,
            unsigned short int * const pBin,
            unsigned int * pWorkingLimit,
            unsigned int totalBeadsInSphere,
            const unsigned int * pLatticePointToRemove){

        auto pBeginIt = bead_indices.begin();
        auto itIndex = std::find(pBeginIt, pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
        // remove original from P(r)
        // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
        removeFromPr(*pLatticePointToRemove, bead_indices, *pWorkingLimit, pBin, totalBeadsInSphere, pBinCount);
        // reduce the workingLimit
        // if wl = 10
        // 0 1 2 3 4 5 6 7 8 9 10
        // remove 4, wl-=1
        // 0 1 2 3 9 5 6 7 8 4 10
        // sort to 9
        // 0 1 2 3 5 6 7 8 9 4 10
        //
        *pWorkingLimit -= 1;
        // swap selected point to the workingLimit
        std::iter_swap(itIndex, pBeginIt + *pWorkingLimit);
        // still need to sort, swap changes the order
        // to save some time, sort should only be from swapped point to workingLimit
        std::sort(pBeginIt, pBeginIt + *pWorkingLimit);
    }

    /**
      * requires a sorted beadsInuse vector
      * removeMe is the index of the bead from origination
    */
    inline void removeFromPr(const unsigned int removeMe, std::vector<unsigned int> & beadsInUse, const unsigned int upperLimit, unsigned short int * const pBin, const unsigned int & totalBeadsInUniverse, std::vector<unsigned int> & prBins ){
        // beads in Use is sorted
        unsigned int row;
        unsigned int row2;

        unsigned int i=0;

        //const int * ptr = &(beadsInUse.front());
        const unsigned int * ptr = beadsInUse.data();

        while(ptr[i] < removeMe && i < upperLimit){
            row = ptr[i];
            row2 = (row*totalBeadsInUniverse - (row*(row+1)/2)) + (removeMe - row);// - 1;
            prBins[ *(pBin + row2 - 1 ) ]--;
            i++;
        }

        // when loop ends beadsInUse[i] => removeMe
        i++;
        // remove across row
        row2 = removeMe*totalBeadsInUniverse - removeMe*(removeMe+1)/2 - removeMe; // assume this is always positive
        while (i < upperLimit){
            //prBins[ *(pBin + row2 + beadsInUse[i] - 1 ) ]--;
            prBins[ *(pBin + (row2 + ptr[i]) - 1 ) ]--;
            i++;
        }
    }


    float totalInterSubUnitDistances(const vector3 * const comvector, std::vector<SubUnit> & subunits, Model *pModel);

    inline void restoreRemovingLatticePointFromBackUp(std::vector<unsigned int>::iterator * pBeginIt,
                                                              unsigned int * pWorkingLimit,
                                                              std::vector<unsigned int> * pBinCountBackUp,
                                                              std::vector<unsigned int>::iterator * pBinCountBegin){
        // value we want to recover is at wl = 9
        // wl += 1 => 10
        // 0 1 2 3 9 5 6 7 8 4 10
        // sort to wl
        // 0 1 2 3 4 5 6 7 8 9 10
        //
        *pWorkingLimit += 1;
        // if I have 5000 lattice points copy from backup is O(n) versus sorting to number of lattice points in working limit
        // sorting is n*log(n) with n = 200?  should be much smaller
        std::sort(*pBeginIt, *pBeginIt+*pWorkingLimit); // swapped index is at workingLimit
        std::copy(pBinCountBackUp->begin(), pBinCountBackUp->end(), *pBinCountBegin); //copy to bin count
    }


    unsigned int removeLatticePositionToModelSym(
                                        std::vector<unsigned int> & bead_indices,
                                        std::vector<unsigned int> & pBinCount,
                                        unsigned int * pWorkingLimit,
                                        const unsigned int * pLatticePointToRemove, Model * pModel, Data *pData);

    unsigned int removeLatticePositionToModelSymHelical(
            std::vector<unsigned int> & bead_indices,
            std::vector<unsigned int> & pBinCount,
            unsigned int * pWorkingLimit,
            const unsigned int * pLatticePointToRemove, Model * pModel, Data *pData);



    inline float updateCVXIndices(const unsigned int & subUnitWorkingLimit, std::vector<unsigned int> &inUse,
                                          std::vector<unsigned int> &cvxIndices, Model *pModel) {

        unsigned int workingLimit = subUnitWorkingLimit;
        char flags[] = "qhull FA";
        coordT points[3*subUnitWorkingLimit];
        std::vector<unsigned int> active_indices(subUnitWorkingLimit);

        for (unsigned int i = 0; i < subUnitWorkingLimit; i++) {
            unsigned int pInt = inUse[i];
            beadToPoint(&points[i*3], pModel->getBead(pInt));
            active_indices[i] = pInt;
        }

        // needs to be optimized
        qh_new_qhull(3, workingLimit, points, 0, flags, nullptr, nullptr);

        auto volume_test = (unsigned int)(qh totvol);

        vertexT * vertices = qh vertex_list;
        unsigned int totalV = (unsigned int)(qh num_vertices);

        // only move CVX hull points
        cvxIndices.resize(totalV);
        for (unsigned int v = 0; v < totalV; v++) { //
            cvxIndices[v] = active_indices[qh_pointid(vertices->point)];
            vertices = vertices->next;
        }
        qh_freeqhull(true);

        return (float)volume_test;
    }

};

/**
 * calculates connectivity potential as a simple quadratic
 * @param mainConnectivity
 * @return
 */
inline float Anneal::connectivityPotentialPhases(unsigned int mainConnectivity){
    float currentConnectivityPotential = (mainConnectivity-1.0f)*(mainConnectivity-1.0f);
    float temp;
    for(unsigned int i=0; i < components.size(); i++) {
        //std::cout << i << " component tour => "<< components[i].getTotalNumberOfComponents() << std::endl;
        temp = components[i].getTotalNumberOfComponents()-1.0f;
        currentConnectivityPotential += temp*temp;
    }
    return currentConnectivityPotential;
}


inline float Anneal::biCameralVolumePotential(unsigned int lowerBound, unsigned int upperBound, float value){

    float diff = 0.01;
    if (value < lowerBound) {
        diff += (lowerBound - value)/(float)lowerBound;
    } else if (value > upperBound){
        diff += (value - upperBound)/(float)upperBound;
    }
    return diff;//*diff;
}





#endif //SUGATAMA_ANNEAL_H
