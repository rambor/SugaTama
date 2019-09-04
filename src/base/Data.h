//
// Created by xos81802 on 30/04/2018.
//

#ifndef SUGATAMA_DATA_H
#define SUGATAMA_DATA_H

#include <string>
#include <vector>
#include <cfloat>
#include "Datum.h"
#include <iostream>
#include <cstdio>
#include <utility>
#include <algorithm>

#include <stdexcept>
#include <fstream>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>


class Data {

private:
    std::string filename;
    unsigned int totalDataPoints, totalPhases, zeroBin;
    std::vector<Datum> intensities; //intensities
    std::vector<Datum> workingSet;
    std::vector<float> invVarianceWorkingSet;
    std::vector<unsigned int> pointsPerBin;
    std::vector<float> lowerBounds;
    std::vector <float> qvalues;
    float volume, rg, rave;
    float qmin, qmax, ns_dmax, dmax, shannon_bins, weight=1.0;
    float bin_width, invBinSize;

    bool isIntensity=false, isPr=true;

    std::vector<float> probability_per_bin;
    std::vector<float> bin_coefficients;
    std::vector<float> working_probability_per_bin;
    std::vector<float> sinc_qr;

    struct RealSpace {
        float r, pr, sigma;
        RealSpace(float rvalue, float pvalue, float svalue) : r(rvalue), pr(pvalue), sigma(svalue){}
    };

    std::vector<RealSpace> pofr;

    bool checkPofRFile(std::string file);
    bool checkIofQFile(std::string file);
    void createCrossValidationDataset();
    void extractIntensityData();
    void extractPrData();
    void parseBins();
    void normalize(float norm);
    void normalizePrBins();
    void normalizePofR(unsigned int count);
    void setDataBinSize(unsigned int bins);
    unsigned int getPointsPerBin(float ioversigma, unsigned int totalInBin);

    long int findInArray(const std::string & value, std::vector<std::string> * strings);

    std::string getFileExt(const std::string& s);

    void logger(std::string description, std::string value) {

        unsigned int len = 40 - description.size();
        std::string firsthalf = std::string(len, ' ');
        printf("%s%s : %s\n", firsthalf.c_str(), description.c_str(), value.c_str());

    }

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
            default:
                sprintf (buffer, "%.6f", number);
        }

        return std::string(buffer);
    }

public:

    Data(std::string filename);
    Data(std::string filename, float fdmax, float fqmax);

    Data(const Data &dataToCopy); // copy constructor

    ~Data(){
        score = nullptr;
        delete score;
    }


    Data & operator=(const Data & dataToCopy) {

        Data tmp(dataToCopy);
        tmp.swap(*this);
////        if (isIntensity){
//            score = new ReciSpaceScore(probability_per_bin.size(), workingSet, sinc_qr, invVarianceWorkingSet);
//        } else if (isPr){
//            score = new RealSpaceScore(zeroBin, invBinSize, working_probability_per_bin);
//        }

        return *this;
    }


    /**
     * Rule of 3.5, define copy, destructor and assignment operator
     * @param other
     */
    void swap(Data & other) noexcept {

        std::swap(filename, other.filename);
        std::swap(totalDataPoints, other.totalDataPoints);
        std::swap(totalPhases, other.totalPhases);
        std::swap(zeroBin, other.zeroBin);
        std::swap(volume, other.volume);
        std::swap(rg, other.rg);
        std::swap(rave, other.rave);
        std::swap(qmin, other.qmin);
        std::swap(qmax, other.qmax);
        std::swap(ns_dmax, other.ns_dmax);
        std::swap(dmax, other.dmax);
        std::swap(shannon_bins, other.shannon_bins);
        std::swap(weight, other.weight);
        std::swap(bin_width, other.bin_width);
        std::swap(invBinSize, other.invBinSize);
        std::swap(isIntensity, other.isIntensity);
        std::swap(isPr, other.isPr);

        std::swap(intensities, other.intensities);
        std::swap(workingSet,other.workingSet);
        std::swap(invVarianceWorkingSet,other.invVarianceWorkingSet);
        std::swap(pointsPerBin,other.pointsPerBin);
        std::swap(lowerBounds,other.lowerBounds);
        std::swap(qvalues,other.qvalues);
        std::swap(probability_per_bin,other.probability_per_bin);
        std::swap(bin_coefficients,other.bin_coefficients);
        std::swap(working_probability_per_bin,other.working_probability_per_bin);
        std::swap(sinc_qr,other.sinc_qr);
        std::swap(pofr,other.pofr);

        // pointer to the score that is on the heap
        std::swap(score, other.score);
    }

//    Data(Data && other) {
//        swap(*this, other);
//    }



    bool getIsIntensity() const;
    bool getIsPr() const;
    float getRg(){ return rg;}
    double getBinWidth(){return (double)bin_width;}
    const float getDmax() const {return dmax;}

    const std::string getDataFilename() const { return filename; }

    unsigned int getZeroBin(){return zeroBin;}

    float getProbabilityPerBin(int bin){return probability_per_bin[bin];}

    void setScoringFunction(unsigned int maxBin);

    double getScore(std::vector<unsigned int> &modelPR);

    unsigned int getShannonBins(){return (unsigned int)shannon_bins;}

    void printKLDivergence(std::vector<unsigned int> &modelPR);

    void printICalc(std::vector<unsigned int> &modelPR);

    float getVolume() {return volume;}

    /**
     * converts distances to bins, can be larger than dmax
     * used for populating vector of pre-assigned bins only
     */
    inline unsigned short int convertToBin(float distance){

        float ratio = distance/bin_width;
        auto floored = (float)floor(ratio);
        unsigned short int binlocale;
        /*
         * lower < distance <= upper
         * bin = 0 is the distance between two beads
         * due to floating point errors, sometimes distance will be slightly larger than 1
         * need to subtract 1 to make it zero
         */
        float diff = std::abs(ratio-floored);

        /*
         * Since we are comparing floats for binning, need to look at relative epsilon
         * We assume ratio is always >= floored
         * bead distances that are equal to the upper limit of a bin-width need to be rounded down and decremented by 1
         */
        if (diff <= ratio*100*FLT_EPSILON ){ // for ratios that are at the border
            binlocale = (floored > 0) ? ((unsigned short int)floored - 1) : 0;

//            if (binlocale > 0){
//                std::cout << bin_width <<  " ratio " << ratio << "  " << floored << " => " << binlocale << " distance " << distance << " " << dmax << std::endl;
//                printf("    ratio : %.5E floored => %i  TIME : %.4f\n", ratio, floored);
//                exit(0);
//            }
        } else {
            binlocale = (unsigned short int)floored;
        }

        return binlocale;
    }

    /**
* converts distances to bins, can be larger than dmax
*/
    inline unsigned int convertToBinUnsignedInt(float distance){

        float ratio = distance/bin_width;
        float floored = floor(ratio);
        unsigned int binlocale=0;
        /*
         * lower < distance <= upper
         * bin = 0 is the distance between two beads
         * due to floating point errors, sometimes distance will be slightly larger than 1
         * need to subtract 1 to make it zero
         */
        float diff = std::abs(ratio-floored);
        if (diff <= ratio*100*FLT_EPSILON ){
            binlocale = (floored > 0) ? ((unsigned int)floored - 1) : 0;
        } else {
            binlocale = (unsigned int)floored;
        }

        return binlocale;
    }

    class Score {
    public:
        Score(){};
        virtual double operator() (std::vector<unsigned int> & bins) = 0;
    };

    Score * score;// = nullptr;

    class ReciSpaceScore : public Score {
    private:
        unsigned int maxBin;
        unsigned int totalInWorkingSet;
        double invTotal;
        std::vector<float> icalc;
        std::vector<Datum> & pWorkingSet;
        std::vector<float> & pSincQrValues;
        std::vector<float> & pInvVar;

    public:
        ReciSpaceScore(unsigned int maxBin, std::vector<Datum> & workingSet, std::vector<float> & qrvalues, std::vector<float> & invVarValues) :  maxBin(maxBin), pWorkingSet(workingSet), pSincQrValues(qrvalues), pInvVar(invVarValues){
            totalInWorkingSet=workingSet.size();
            invTotal = 1.0d/(double)totalInWorkingSet;
            icalc.resize(totalInWorkingSet);
        };

        virtual double operator() (std::vector<unsigned int> & modelPR) {

            // calculate chi2
            float value, diff, top=0, bottom=0;

            std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());
            double totalCounts = 0.0d;
            for (auto & pr : modelPR_float) {
                totalCounts += pr;
            }
            double invTotalCounts = 1.0d/totalCounts;

            // for each q-value, calculate I_calc from P(r)
            for (unsigned int qr = 0; qr<totalInWorkingSet; qr++) {
                float iofq=0;
                unsigned int count=0;
                for (auto & pr : modelPR_float) { // iterate over r-values
                    iofq += pSincQrValues[qr*maxBin + count]*pr*invTotalCounts;
                    count++;
                }

                icalc[qr] = iofq;
                value = pInvVar[qr];
                top += iofq*pWorkingSet[qr].getI()*value;
                bottom += iofq*iofq*value;
            }

            FILE * pFile;
            const char *outputFileName;
            std::string nameOf = "icalc.txt";
            outputFileName = nameOf.c_str() ;
            pFile = fopen(outputFileName, "w");

            std::cout << " ICALC" << std::endl;
            float chi2=0;
            float scale = top/bottom;
            for (unsigned int qr = 0; qr<totalInWorkingSet; qr++) {
                diff = pWorkingSet[qr].getI() - scale*icalc[qr];
                fprintf(pFile, "%.6E %.5E  %.5E \n", pWorkingSet[qr].getQ(), scale*icalc[qr], pWorkingSet[qr].getI());
                chi2 += diff*diff*pInvVar[qr];
            }
            fclose(pFile);

            return chi2*invTotal;
        }
    };


    /*
     * Functor to calculate score based on Pr distributin
     * Functor is set in Annealing program
     */
    class RealSpaceScore : public Score {

    private:
        unsigned int zeroBin;
        double invBin;
        std::vector<float> & pWorkingProbPerBin;

    public:
        RealSpaceScore(unsigned int zeroBin, double invBinSize, std::vector<float> & it) : zeroBin(zeroBin), invBin(invBinSize), pWorkingProbPerBin(it){};

        virtual double operator() (std::vector<unsigned int> & modelPR) {

            double totalCounts = 0.0;
            double kl=0.0;
            unsigned int emptyBins = 0;
            double prob, *value;

            std::vector<double> modelPR_float(modelPR.begin(), modelPR.end());

            for (auto & pr : modelPR_float) {
                totalCounts += pr;
            }
            //double invTotal = 1.0/totalCounts;
            for (unsigned int i=0; i < zeroBin; i++){
                // i know every value in working_probability up to zeroBin is nonzero
                value = &modelPR_float[i];
                if (*value > 0) {
                    prob = pWorkingProbPerBin[i];  // bounded by experimental Shannon Number
                    kl += prob * std::log(prob / (*value) * totalCounts);
                } else if (*value == 0 ){ // if any bins are zero before dmax bin
//                    kl += 1.10/invBin; // must penalize empty bins
                    emptyBins += 1;
                }
            }

            // zeroBin should have at least one value
//            if (modelPR_float[zeroBin] == 0){
//                kl += 0.1;
//            }
            double penalty = 0.0d;
            if (emptyBins == 1){
                penalty = 1.1;
            } else if (emptyBins > 1){
                penalty = std::pow(10,emptyBins);
            }

            return kl*invBin + penalty;
        }
    };

    float getQmax(){ return qmax;}


};


#endif //SUGATAMA_DATA_H
