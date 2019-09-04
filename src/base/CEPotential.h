//
// Created by xos81802 on 20/09/2018.
//

#ifndef SUGATAMA_CEPOTENTIAL_H
#define SUGATAMA_CEPOTENTIAL_H

#include <SubUnit.h>
#include <iostream>

class CEPotential {
//    SubUnit * pSubUnit1;
//    SubUnit * pSubUnit2;

    unsigned int subUnit1Index, subUnit2Index;

    std::vector<float> probabilitiesPerBin;
    std::vector<float> distances;
    std::vector<float> cdf;

    float binWidth;
    float alpha=0.63;
    unsigned int totalBins;

    unsigned int currentBin;
    float currentDistance, upperLimit, lowerLimit;

    void updateCDF();

public:

    CEPotential();
    CEPotential (SubUnit & psub1, SubUnit & psub2, float binWidth, unsigned int totalBins);
    CEPotential (int sub1, int sub2, float distance, float lower, float upper);
/**
     * Copy Constructor
     * @param p2
     */
    CEPotential(const CEPotential &p2) {
        this->currentBin = p2.currentBin;
        this->currentDistance = p2.currentDistance;
        this->upperLimit = p2.upperLimit;
        this->lowerLimit = p2.lowerLimit;

//        this->pSubUnit1 = &(*p2.pSubUnit1);
//        this->pSubUnit2 = &(*p2.pSubUnit2);

//        this->pSubUnit1 = *(&p2.pSubUnit1);
//        this->pSubUnit2 = *(&p2.pSubUnit2);

        this->subUnit1Index = p2.subUnit1Index;
        this->subUnit2Index = p2.subUnit2Index;

        this->totalBins = p2.totalBins;
        this->probabilitiesPerBin.resize(this->totalBins);
        this->distances.resize(this->totalBins);
        this->cdf.resize(this->totalBins);

        this->binWidth = p2.binWidth;
        this->alpha=p2.alpha;

        std::copy(p2.probabilitiesPerBin.begin(),
                  p2.probabilitiesPerBin.end(),
                  this->probabilitiesPerBin.begin());

        std::copy(p2.distances.begin(),
                  p2.distances.end(),
                  this->distances.begin());

        std::copy(p2.cdf.begin(),
                  p2.cdf.end(),
                  this->cdf.begin());
    }

    void updateProbabilities(std::vector<unsigned int> & counts, float totalCounts);
    void seedProbabilities(float distance);

    unsigned int setConfiguration();

    //SubUnit * getSubUnit1(){return pSubUnit1;}
    unsigned int getSubUnit1(){ return subUnit1Index;}
    unsigned int getSubUnit2(){ return subUnit2Index;}

    float getCurrentDistance(){ return currentDistance;}

    //SubUnit * getSubUnit2(){return pSubUnit2;}

    void printDistribution(int round);

    float getMostProbableDistance(){
        float max=probabilitiesPerBin[0];
        unsigned int maxBin=0;
        for (unsigned int i=1; i<this->totalBins; i++){
            if (max < probabilitiesPerBin[i]){
                max = probabilitiesPerBin[i];
                maxBin = i;
            }
        }
        return distances[maxBin];
    };




    void setMostProbableDistance(){
//        float max=probabilitiesPerBin[0];
        float weightedAverage=0.0f;
//        unsigned int maxBin=0;
//        for (unsigned int i=1; i < this->totalBins; i++){
//
//            weightedAverage += distances[i]*probabilitiesPerBin[i];
//
//            if (max < probabilitiesPerBin[i]){
//                max = probabilitiesPerBin[i];
//                maxBin = i;
//            }
//        }


        for (unsigned int i=0; i < this->totalBins; i++){
            weightedAverage += distances[i]*probabilitiesPerBin[i];
        }


//        currentDistance = distances[maxBin];
//        upperLimit = currentDistance + 0.5*binWidth;
//        lowerLimit = currentDistance - 0.5*binWidth;

        currentDistance = weightedAverage;
        upperLimit = currentDistance + 0.5f*binWidth;
        lowerLimit = currentDistance - 0.5f*binWidth;

    };

    void setTest(float value){

        currentDistance = std::floor(value/binWidth)*binWidth + 0.5f*binWidth;

        //double oldprob;
        float sum = 0.0f;
        for(unsigned int i=0; i<totalBins; i++){
            // get old probability by index
            probabilitiesPerBin[i] = 0.000001;
            sum += probabilitiesPerBin[i];
        }
        sum -= probabilitiesPerBin[std::floor(value/binWidth)];
        probabilitiesPerBin[ std::floor(value/binWidth)] = 0.999;
        sum += probabilitiesPerBin[std::floor(value/binWidth)];

        // renormalization
        float invsum = 1.0f/sum;
        for(unsigned int i=0; i<totalBins; i++){
            probabilitiesPerBin[i] *= invsum;
        }

        updateCDF();

        upperLimit = currentDistance + 0.5f*binWidth;
        lowerLimit = currentDistance - 0.5f*binWidth;
    };

    bool keepIt(float dis);

    inline float calculateEnergy(float distance){

        bool right = distance < upperLimit;
        bool left = distance > lowerLimit;

        if (right && left){
            return 0.0f;
        } else {
            if (right){
                float diff = lowerLimit - distance;
                return diff*diff;
            } else {
                float diff = distance - upperLimit;
                return diff*diff;
            }
        }
    }

    inline float calculatePotentialEnergy(float distance){
        float diff = distance - currentDistance;
//        return diff*diff;
        if (distance > lowerLimit && distance < upperLimit){
            //return diff*diff;
            return 0;
        } else {
//            float diffU = distance - upperLimit;
//            float diffL = lowerLimit - distance;
//            float diff = ( diffU > diffL) ? diffU : diffL;
            return diff*diff;
        }
    }
};


#endif //SUGATAMA_CEPOTENTIAL_H
