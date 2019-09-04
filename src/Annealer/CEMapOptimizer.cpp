//
// Created by xos81802 on 09/02/2019.
//


#include "../Anneal.h"

void Anneal::ceMapOptimization(Model *pModel, Data *pData, unsigned int topN, const unsigned int minN, const unsigned int maxN, std::vector<KDE::ProbabilityBead> & lattice){

    std::clock_t startTime;
    double runtime;
    std::random_device rd;
    std::mt19937 gen(rd());
    unsigned int currentN, indexOf;

    std::map<unsigned int, unsigned int> indicesInUse; // use shared pointer? make instance on heap
    for(auto & point : lattice){
        indicesInUse.emplace(point.index, 0.0d);
    }


    maxbin = pModel->getMaxBin(pData);
    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");

    maxbin += 1;  // maximum bin index, so vector size must be +1

    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
    contactCutOff = interconnectivityCutOff;

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    std::uniform_int_distribution<unsigned int> randomIndex(minN,maxN-1); // guaranteed unbiased
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    std::vector<Trial> topTrials;
    topTrials.reserve(topN);
    unsigned int totalInLattice = lattice.size();
    const unsigned int last = topN-1;
    unsigned int halfOf = totalInLattice/2;
    unsigned int diff, pickMe, neighbor;
    std::sort(lattice.begin(), lattice.end());

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<unsigned int> binCount(maxbin);
    std::vector<unsigned int> bead_indices(totalInLattice);
    std::vector<KDE::ProbabilityBead> latticeInUse(totalInLattice);

    std::vector<unsigned int> trialX;
    std::set<unsigned int> trialInUseSet;
    std::set<unsigned int> hull;
    trialX.reserve(halfOf);

    std::vector<double> running_average(3);
    std::fill(running_average.begin(), running_average.end(), 0.0d);

    double prior = 0;
    std::cout << "*******************             MAP REFINEMENT             *******************" << std::endl;
    logger("MAX BIN", std::to_string(maxbin));
    logger("BOUNDS LOWER", std::to_string(minN));
    logger("BOUNDS UPPER", std::to_string(maxN));
    logger("TOTAL ROUNDS", std::to_string(highTempRounds));
    logger("TRIALS PER ROUND", std::to_string(ccmultiple));
    logger("TOPN", std::to_string(topN));
    logger("CONTACT CUTOFF (Angstroms)", formatNumber(contactCutOff, 2));
    std::cout << "*******************       CROSS ENTROPY OPTIMIZATION       *******************" << std::endl;
    for(unsigned int round=0; round < highTempRounds; round++){

        unsigned int topAdded=0;
        startTime = std::clock();
        for(unsigned int trial=0; trial<ccmultiple; trial++){
            trialX.clear();
            trialInUseSet.clear();
            std::copy(lattice.begin(), lattice.end(), latticeInUse.begin()); // 2to3x faster than new instance
            //std::vector<KDE::ProbabilityBead> latticeInUse(lattice);
            unsigned int * const ptr = bead_indices.data();
            /*
             * binary method
             */
            // should the order of LatticeInUse be randomized?
            std::shuffle(latticeInUse.begin(), latticeInUse.end(), gen);
            currentN=randomIndex(gen);
            indexOf=0;

            for(auto & pLattice : latticeInUse){
                if (distribution(gen) < pLattice.prob){
                    trialX.emplace_back(pLattice.index);
                    ptr[indexOf] = pLattice.index;
                    indexOf++;
                    if (indexOf == currentN){
                        break;
                    }
                }
            }
            currentN = indexOf;

            /*
             * random start and connected random walk (Euler tour of 1)
             */
//            for(auto & pLattice : latticeInUse){
//                if (distribution(gen) < pLattice.prob){
//                    trialX.emplace_back(pLattice.index);
//                    ptr[indexOf] = pLattice.index;
//                    indexOf++;
//                    trialInUseSet.clear();
//                    hull.clear();
//                    trialInUseSet.insert(pLattice.index);
//                    populateLayeredDeadlimitUsingSet(trialInUseSet, hull, pModel);
//                    break;
//                }
//            }
//
//            std::uniform_int_distribution<unsigned int> randomHull(0,hull.size()-1); // guaranteed unbiased
//            std::vector<unsigned int> hullIndices(hull.begin(), hull.end());
//            // now add based on neighborhood
//            unsigned int quitNow = 0;
//            while (indexOf < currentN || quitNow < totalInLattice){
//                // get neighbor to add
//                for(auto & neighbor : hullIndices){
//
//                    auto fit = std::find_if(latticeInUse.begin(), latticeInUse.end(), find_pt_by_key(neighbor));
//                    if (distribution(gen) < (*fit).prob && trialInUseSet.find(neighbor) == trialInUseSet.end() ){
//                        trialX.emplace_back((*fit).index);
//                        ptr[indexOf] = (*fit).index;
//                        indexOf++;
//                        trialInUseSet.insert((*fit).index);
//                        populateLayeredDeadlimitUsingSet(trialInUseSet, hull, pModel);
//                        quitNow = 0;
//                        randomHull = std::uniform_int_distribution<unsigned int>(0,hull.size()-1);
//                        hullIndices = std::vector<unsigned int>(hull.begin(), hull.end());
//                        std::shuffle(hullIndices.begin(), hullIndices.end(), gen);
//                        break;
//                    }
//                    quitNow++;
//                }
//            }
//            currentN = indexOf;

            //score it D_kl, contacts distribution, connectivity?
            calculateModelPrDistributionDirect(&bead_indices, &binCount, currentN, pModel, pData);
            double currentKL = pData->getScore(binCount);

            std::set<unsigned int>beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + currentN);
//            EulerTour tour(beads_in_use_tree, pModel);
//            diff = tour.getNumberOfComponents() - 1;

            populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
            double currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

            //auto current_energy = (currentKL + 0.0001d*currentKLDivContacts);
            auto current_energy = currentKL;

            //update best list
            if (topAdded < topN){
                topTrials.emplace_back(Trial(current_energy, trialX));
                topAdded++;
                std::sort(topTrials.begin(), topTrials.end());
            } else {
                if (current_energy < topTrials[last].value){
                    /*
                     * replace last entry and sort
                     */
                    Trial * pTrial = &topTrials[last];
                    pTrial->value = current_energy;
                    pTrial->indices = std::vector<unsigned int>(trialX);
                    std::sort(topTrials.begin(), topTrials.end());
                }
            }
        }

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
        logger("TRIAL TIME", formatNumber(runtime,6));

        double diff = topTrials[last].value - topTrials[0].value;
        double ave = (topTrials[0].value + topTrials[last].value)*0.5;

        double average = (running_average[0] + running_average[1] + running_average[2])/3.0d;
        double epsilon = std::abs(average - topTrials[last].value)/average;

        std::printf(" RND: %2i 1st: %.1E LAST: %.1E GAP: %.1E EPSLON: %.2f \n", round, topTrials[0].value, topTrials[last].value, diff, epsilon);

        running_average[0] = running_average[1];
        running_average[1] = running_average[2];
        running_average[2] = topTrials[last].value;

        // update CDF
        for(auto & trial : topTrials){
            for(auto & selectedIndex : trial.indices){
                (*indicesInUse.find(selectedIndex)).second++;
            }
        }

        // update lattice probabilities
        //double oldprob, invTotal = 1.0d/totalModels;
        double oldprob, invTotal = 1.0d/(double)topN;
        double updateAlpha=0.67;

        //probabilities should sum to 1 over all beads
        double maxprob = 0, minprob = 1;
        for(auto & point : lattice){
            auto mit = indicesInUse.find(point.index);
            oldprob = (1.0 - updateAlpha)*point.prob;

            if (mit != indicesInUse.end()){
                point.prob = updateAlpha*(*mit).second*invTotal + oldprob;
                if (point.prob > maxprob){
                    maxprob = point.prob;
                }
            } else {
                point.prob = oldprob;
            }

            if (point.prob < minprob){
                minprob = point.prob;
            }
            point.prior = point.prob;
        }

        if (epsilon < 0.007 && round > 5){
            std::string line = "EPSILON " + formatNumber(epsilon,3) + " < 0.007";
            logger("TERMINATING",line);
            break;
        }


        logger("MAX PROBABILITY", formatNumber(maxprob,3));
        logger("MIN PROBABILITY", formatNumber(minprob,3));

        // logger("MAX DELTA PROBABILITY", formatNumber(maxdiff,5));
        //std::sort(lattice.begin(), lattice.end());

        for(auto & mit : indicesInUse){
            mit.second = 0;
        }
        topTrials.clear();
        topTrials.reserve(topN);
    }
}



void Anneal::ceMapOptimizationSym(Model *pModel, Data *pData, unsigned int topN, const unsigned int minN, const unsigned int maxN, std::vector<KDE::ProbabilityBead> & lattice){

    std::clock_t startTime;
    double runtime, currentKL;
    std::random_device rd;
    std::mt19937 gen(rd());
    unsigned int currentN, indexOf;


    contactCutOff = pModel->getNieghborCutOffLimit();
    violation_limit = sqrt(3)*pModel->getBeadRadius();

    std::cout << " cutoffs " << contactCutOff << " " << interconnectivityCutOff<< std::endl;
    //contactCutOff = interconnectivityCutOff;
    violation_limit = sqrt(3)*pModel->getBeadRadius();

    std::map<unsigned int, unsigned int> indicesInUse; // use shared pointer? make instance on heap
    std::vector<vector3> allSubunits(pModel->getNumberOfSubUnits()*lattice.size());
    std::map<unsigned int, std::vector<vector3> > allSubunitsMap;

    // pre-calculate the symmetry space since we are using a defined set of lattice points
    for(auto & point : lattice){
        indicesInUse.emplace(point.index, 0);

        const auto pVec = pModel->getBead(point.index)->getVec();
        auto mit = allSubunitsMap.emplace(point.index, std::vector<vector3>(pModel->getNumberOfSubUnits()));
        (*mit.first).second[0] = vector3(pVec);
        vector3 * const pVecArray = (*mit.first).second.data();

        for(unsigned int j=1; j<pModel->getNumberOfSubUnits(); j++){
            pVecArray[j] = pModel->transformVectorBySymmetry(j, pVec);
        }
    }

    maxbin = pModel->getMaxBin(pData);
    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");

    maxbin += 1;  // maximum bin index, so vector size must be +1

    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two


    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    unsigned int totalViolations, last = topN-1;
    std::uniform_int_distribution<unsigned int> randomIndex(minN,maxN-1); // guaranteed unbiased
    std::uniform_real_distribution<double> distribution(0.0,1.0);


    std::vector<Trial> topTrials;
    topTrials.reserve(topN);
    unsigned int totalInLattice = lattice.size();
    unsigned int halfOf = totalInLattice/2;
    std::sort(lattice.begin(), lattice.end());

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<unsigned int> binCount(maxbin);
    std::vector<unsigned int> bead_indices(totalInLattice);
    std::set<unsigned int>beads_in_use_tree;
    std::vector<KDE::ProbabilityBead> latticeInUse(totalInLattice);

    std::set<unsigned int> inUse;
    std::vector<unsigned int> trialX;
    trialX.reserve(halfOf);

    double prior_first=0, prior_last=0;
    std::vector<double> running_average(3);
    std::fill(running_average.begin(), running_average.end(), 0.0d);

    std::cout << "*******************             MAP REFINEMENT             *******************" << std::endl;
    logger("SYMMETRY", pModel->getSymmetryString());
    logger("MAX BIN", std::to_string(maxbin));
    logger("BOUNDS LOWER", std::to_string(minN));
    logger("BOUNDS UPPER", std::to_string(maxN));
    logger("TOTAL ROUNDS", std::to_string(highTempRounds));
    logger("TRIALS PER ROUND", std::to_string(ccmultiple));
    logger("TOPN", std::to_string(topN));
    logger("CONTACT CUTOFF (Angstroms)", formatNumber(contactCutOff, 2));
    logger("VIOLATION LIMIT", formatNumber(violation_limit, 2));
    std::cout << "*******************       CROSS ENTROPY OPTIMIZATION       *******************" << std::endl;

    for(unsigned int round=0; round < highTempRounds; round++){

        unsigned int topAdded=0;
        startTime = std::clock();
        std::clock_t  startTime2 = std::clock();

        for(unsigned int trial=0; trial<ccmultiple; trial++){
            trialX.clear();
            std::copy(lattice.begin(), lattice.end(), latticeInUse.begin()); // 2to3x faster than new instance

            std::shuffle(latticeInUse.begin(), latticeInUse.end(), gen);

            currentN = randomIndex(gen);

            indexOf=0;
            unsigned int * const ptr = bead_indices.data();

            for(auto & pLattice : latticeInUse){ // iterate through the whole set, accept with probability
                if (distribution(gen) < pLattice.prob){
                    trialX.emplace_back(pLattice.index);
                    ptr[indexOf] = pLattice.index;
                    indexOf++;
                    if (indexOf == currentN){
                        break;
                    }
                }
            }
            currentN = indexOf;

//            unsigned int stopAt = lattice.size(),  startAt=0, lastPlace = stopAt-1;
////            currentN = randomIndex(gen);
//            inUse.clear();
//            while( indexOf < maxN ){
//                unsigned int l=startAt;
//                for(; l<stopAt; l++){
//                    auto * pLattice = &latticeInUse[l];
//
//                    if (distribution(gen) < pLattice->prob && inUse.find(pLattice->index) == inUse.end()){
//                        trialX.emplace_back(pLattice->index);
//                        ptr[indexOf] = pLattice->index;
//                        inUse.insert(pLattice->index);
//                        indexOf++;
//                        std::iter_swap(latticeInUse.begin()+l, latticeInUse.begin()+lastPlace);
//                        lastPlace--;
//                        stopAt--;
//                        startAt=l;
//                        break;
//                    }
//                }
//                startAt = (startAt == stopAt-1) ? 0 : startAt;
//                if (indexOf >= minN && distribution(gen) > 0.57){
//                    break;
//                }
//            }
//            currentN = indexOf;

//            unsigned int * const ptr = bead_indices.data();
//            indexOf = 0;
//            std::shuffle(latticeInUse.begin(), latticeInUse.end(), gen);
//            for(auto & pLattice : latticeInUse){
//                if (distribution(gen) < pLattice.prob){
//                    trialX.emplace_back(pLattice.index);
//                    ptr[indexOf] = pLattice.index;
//                    indexOf++;
//                }
//            }
//            currentN = indexOf;

//            for(unsigned int index=0; index<currentN; index++){
//                std::uniform_real_distribution<double> distribution(0.0,upperProbLimit);
//                double locale = distribution(gen);
//                double sumIt = 0.0d;
//
//                for(auto it = latticeInUse.begin(); it != latticeInUse.end(); ++it){
//                    sumIt += (*it).prob;
//                    if (sumIt >= locale){
//                        trialX.emplace_back((*it).index);
//                        ptr[index] = (*it).index;
//                        upperProbLimit -= (*it).prob;
//                        std::iter_swap(it, latticeInUse.begin()+lastPlace);
//                        lastPlace--;
//                        break;
//                    }
//                };
//            }

            //score it D_kl, contacts distribution, connectivity?
            beads_in_use_tree.clear();
            beads_in_use_tree = std::set<unsigned int>(bead_indices.begin(), bead_indices.begin() + currentN);

//            startTime = std::clock();
            calculateModelPrDistributionSymCE(&bead_indices, beads_in_use_tree, allSubunitsMap, contactsDistributionOfModel, &binCount, currentN, totalViolations, pModel, pData );
//            runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//            std::cout << " SECOND " << runtime << std::endl;
            currentKL = pData->getScore(binCount);
//            currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

            auto current_energy = (currentKL + beta*totalViolations);

            //update best list
            if (topAdded < topN){
                topTrials.emplace_back(Trial(current_energy, trialX));
                topAdded++;
                std::sort(topTrials.begin(), topTrials.end());
            } else if (current_energy < topTrials[last].value){
                    /*
                     * replace last entry and sort
                     */
                    Trial * pTrial = &topTrials[last];
                    pTrial->value = current_energy;
                    pTrial->indices = std::vector<unsigned int>(trialX);
                    std::sort(topTrials.begin(), topTrials.end());
            }

            if (trial%1999==0 && trial > 0){
                std::string line = std::to_string(trial) + " [ last " + formatNumber(topTrials[last].value,8) +"]";
                logger("TRIAL", line);
                runtime = (std::clock() - startTime2)/(double) CLOCKS_PER_SEC;
                logger("TIME (SEC)", formatNumber(runtime, 3));
                startTime2 = std::clock();
            }
        }



        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
        double diff = topTrials[last].value - topTrials[0].value;
        double ave = (topTrials[0].value+ topTrials[last].value)*0.5;
        double epsilon = diff/ave;

        double average = (running_average[0] + running_average[1] + running_average[2])/3.0d;
        double diffave = std::abs(average - topTrials[last].value)/average;

        printf(" RND: %2i 1st: %.1E LAST: %.1E GAP: %.1E EPSLON: %.2f TIME(SEC): %.1f \n", round, topTrials[0].value, topTrials[last].value, diff, diffave, runtime);

        running_average[0] = running_average[1];
        running_average[1] = running_average[2];
        running_average[2] = topTrials[last].value;

        // update CDF
        for(auto & trial : topTrials){
            for(auto & selectedIndex : trial.indices){
                (*indicesInUse.find(selectedIndex)).second++;
            }
        }

        // update lattice probabilities
        double oldprob, invTotal = 1.0d/(double)topN;
        double updateAlpha=0.61;

        /*
         * update probabilities using a smoothing function
         */
        double maxprob = 0, minprob = 1;
        for(auto & point : lattice){
            auto mit = indicesInUse.find(point.index);
            oldprob = (1.0 - updateAlpha)*point.prob;

            if (mit != indicesInUse.end()){
                point.prob = (updateAlpha*(double)(*mit).second)*invTotal + oldprob;
                if (point.prob > maxprob){
                    maxprob = point.prob;
                }
            } else {
                point.prob = oldprob;
            }

            if (point.prob < minprob){
                minprob = point.prob;
            }
            point.prior = point.prob;
        }


        if (diffave < 0.05 && round > 5){
            std::string line = "EPSILON " + formatNumber(diffave,3) + " < 0.05";
            logger("TERMINATING",line);
            break;
        }
        logger("MAX PROBABILITY", formatNumber(maxprob,3));
        logger("MIN PROBABILITY", formatNumber(minprob,3));

        //std::sort(lattice.begin(), lattice.end());
        // reset count to zero for next round
        for(auto & mit : indicesInUse){
            mit.second = 0;
        }
        topTrials.clear();
        topTrials.reserve(topN);

    }

}

