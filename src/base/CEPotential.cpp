//
// Created by xos81802 on 20/09/2018.
//

#include "CEPotential.h"

CEPotential::CEPotential(){};


CEPotential::CEPotential (int sub1, int sub2, float distance, float lower, float upper){
    //    this->pSubUnit1 = &psub1;
    //    this->pSubUnit2 = &psub2;

    this->subUnit1Index = sub1;
    this->subUnit2Index = sub2;

    this->binWidth = binWidth;
    currentDistance = distance;
    upperLimit = currentDistance + upper;
    lowerLimit = currentDistance - lower;

    this->totalBins=1;
    this->probabilitiesPerBin.resize(totalBins);
    this->distances.resize(totalBins);
    this->cdf.resize(totalBins);

    this->probabilitiesPerBin[0] = 1.0;
    this->distances[0] = distance;

}

CEPotential::CEPotential (SubUnit & psub1, SubUnit & psub2, float binWidth, unsigned int totalBins){

//    this->pSubUnit1 = &psub1;
//    this->pSubUnit2 = &psub2;

    this->subUnit1Index = psub1.getSubUnitIndex();
    this->subUnit2Index = psub2.getSubUnitIndex();

    this->binWidth = binWidth;
    this->totalBins = totalBins;

    /*
     * initialize probabilities
     */
    float prob = 1.0f/(float)totalBins;
    this->probabilitiesPerBin.resize(totalBins);
    this->distances.resize(totalBins);
    this->cdf.resize(totalBins);

//    probabilitiesPerBin[0] = 0.002;
//    probabilitiesPerBin[1] = 0.99;
//    probabilitiesPerBin[2] = 0.002;
//    probabilitiesPerBin[3] = 0.002;
//    probabilitiesPerBin[4] = 0.002;

    for(unsigned int i=0; i<totalBins; i++){
        probabilitiesPerBin[i] = prob;
        distances[i] = 0.5f*binWidth + i*binWidth;
    };
    this->updateCDF();
}

/*
 * use CDF to sample
 * returns currentBin that was used to set Configuration
 */
unsigned int CEPotential::setConfiguration() {

    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> randomIndex(0,1.0f);

    float value = randomIndex(gen);
/*
 * CDF Example value = 0.09
 * [0] =>   0 < x <= 0.1  -> 0.1
 * [1] => 0.1 < x <= 0.3  -> 0.3
 * [2] => 0.3 < x <= 0.5  -> 0.5
 * [3] => 0.5 < x <= 0.7  -> 0.7
 * [4] => 0.7 < x <= 1.0  -> 1.0
 *
 */
    for(unsigned int i=0; i<totalBins; ++i){
        if ( value <= cdf[i] ) {
            currentBin = i;
            break;
        }
    }

    currentDistance = distances[currentBin];
    upperLimit = currentDistance + 0.5f*binWidth;
    lowerLimit = currentDistance - 0.5f*binWidth;
    return currentBin;
};



bool CEPotential::keepIt(float distance) {

    if (distance >= lowerLimit && distance <= upperLimit){
        return true;
    }
    return false;
}


void CEPotential::seedProbabilities(float distance){
    // which bin does distance belong to?
    unsigned int assignedBin = 0;
    for(unsigned int i=0; i<totalBins; i++){
        float lower = i*binWidth;
        float upper = (i+1)*binWidth;
        if (distance > lower && distance <= upper){
            assignedBin = i;
            break;
        }
    };

    std::vector<unsigned int> counts(totalBins);
    std::fill(counts.begin(), counts.end(), 1);
    counts[assignedBin] = 10;
    float totalCounts = 0;
    for(unsigned int i=0; i<totalBins; i++){
        totalCounts += counts[i];
    }
    this->updateProbabilities(counts, totalCounts);
}


void CEPotential::updateProbabilities(std::vector<unsigned int> &counts, float totalCounts) {

    if (totalBins > 1){
        float inv = 1.0f/totalCounts;
        //double oldprob;
        float sum = 0.0f;
        for(unsigned int i=0; i<totalBins; i++){
            // get old probability by index
            probabilitiesPerBin[i] = alpha*counts[i]*inv + (1.0f-alpha)*probabilitiesPerBin[i];
            sum += probabilitiesPerBin[i];
        }

        // renormalization
        float invsum = 1.0f/sum;
        for(unsigned int i=0; i<totalBins; i++){
            probabilitiesPerBin[i] *= invsum;
        }

        this->updateCDF();
    }
}

void CEPotential::printDistribution(int round){
//    std::string tempHeader = "CEPotential\n";
//    char buffer[80];
//    int cstring = snprintf(buffer, 80, "%4i SUBUNITS => %d <-> %d \n", round, pSubUnit1->getSubUnitIndex(), pSubUnit2->getSubUnitIndex());
//    tempHeader.append(buffer);

    float expectedValue=0;
    std::cout << " SUBUNITS " << subUnit1Index << " <-> " << subUnit2Index <<std::endl;
    for(unsigned int i=0; i<totalBins; i++){
        std::printf("  => %3d %5.3f %.4f \n", i, distances[i], probabilitiesPerBin[i]);
        expectedValue += distances[i]*probabilitiesPerBin[i];
    }

    FILE * pFile;
    pFile = fopen("ce_probabilities.txt", "a");

    fprintf(pFile, "%4i SUBUNITS => %d <-> %d \n", round, subUnit1Index, subUnit2Index);
    for(unsigned int i=0; i<totalBins; i++){
        std::string index = std::to_string(i);
        fprintf(pFile, "   %s %6.2f %.4f \n", index.c_str(), distances[i], probabilitiesPerBin[i]);
    }
    fprintf(pFile, " EXPECTATION : %.2f\n", expectedValue);
    fclose(pFile);
}

void CEPotential::updateCDF() {

/*
 * CDF Example value = 0.09
 * [0] =>   0 < x <= 0.1  -> 0.1
 * [1] => 0.1 < x <= 0.3  -> 0.3
 * [2] => 0.3 < x <= 0.5  -> 0.5
 * [3] => 0.5 < x <= 0.7  -> 0.7
 * [4] => 0.7 < x <= 1.0  -> 1.0
 *
 */

    float sum=0.0f;
    for(unsigned int i=0; i<totalBins; i++){
        sum += probabilitiesPerBin[i];
        cdf[i] = sum;
    }

    if (sum > 1.0001f){
        float diff = 1.0f - sum;
        std::cout << " TOO LARGE CDF " << sum << " LOG10 " << std::log10(diff) <<std::endl;
        exit(0);
    }
}