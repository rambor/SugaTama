//
// Created by xos81802 on 07/08/2018.
//
#include <random>
#include "../Anneal.h"
#include "../EulerTour/EulerTour.h"

std::string Anneal::refineHomogenousBodyASAHybridEx(Model *pModel, Data *pData, std::string outputname) {
    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;

    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.000001 : (double)highTempStartForCooling;
    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    const char * pOutputFileName = outputname.c_str() ;
    std::string status = "RAMPING TO TEMP";

    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");

    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> active_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    std::vector<unsigned int>::iterator itIndex, pSwap2;

    // copy Starting_Set from initial model
    pModel->copyStartingModelIntoVector(bead_indices);
    unsigned int workingLimit = pModel->getStartingWorkingLimit();

    // reset iterators to internal bead_indices
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    // set deadLimit of the selected set
    // convert distances in Search Sphere to ShannonBin membership
    this->fillPrBinsAndAssignTotalBin(pModel, pData);
    unsigned short int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    /*
     * create set of beads within CONVEX HULL
     */
    std::set<unsigned int> hull;
//    recalculateDeadLimit(workingLimit, bead_indices, tempHull, pModel, totalBeadsInSphere);
    populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
    std::string hullName = "hull";
    pModel->writeSetToFile(hull, hullName);
    std::set<unsigned int> original_hull(beads_in_use_tree.begin(), beads_in_use_tree.end());

//    for(auto sit = beads_in_use_tree.begin(); sit != beads_in_use_tree.end(); ++sit){
//        original_hull.insert(*sit);
//    }
    const unsigned int targetCount = 0.8*(original_hull.size());

//    std::set<unsigned int> universe(bead_indices.begin(), bead_indices.end());
//    hullName = "uni";
//    pModel->writeSetToFile(universe, hullName);

    /*
     * create cross-validation set or working observed probability distribution that encompasses search space
     */
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);
    //std::vector<double> testContactsDistributionOfModel(13);

    std::cout << "      TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "                BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "             BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    double testKL, currentKL = pData->getScore(binCount);
    double lowestKL = currentKL, startingKL = currentKL;

    // CVX hull calculation
    char flags[] = "qhull FA";
    float starting_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
    /*
     * model for initial search is too dense, so deadlimit layer should be sufficient for a confined search
     */
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;

    // coupon collector's problem
    unsigned int couponbase = hull.size();// + workingLimit;
    auto coupons = (unsigned int)(couponbase*std::log((double)couponbase) + 0.5772156649*couponbase + 0.5d);
    unsigned int updateCount = ccmultiple*coupons;

   // float step_limit = updateCount/0.65f;///0.1;
    float step_limit = updateCount;
    std::vector<float> acceptanceRateDuringRun((unsigned long int)step_limit);
    std::vector<double> tempDuringRun((unsigned long int)step_limit);
    std::vector<double> divergenceDuringRun((unsigned long int)step_limit);
    std::vector<unsigned int> workingLimitDuringRun((unsigned long int)step_limit);

    double * pTempDuringRun = &tempDuringRun.front();
    float * pAcceptanceRateDuringRun = &acceptanceRateDuringRun.front();
    double * pDivergenceDuringRun = &divergenceDuringRun.front();
    unsigned int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents = eulerTour.getNumberOfComponents();
    unsigned int tempNumberOfComponents = currentNumberOfComponents;

    std::cout << "             EULER TOURS : " << currentNumberOfComponents << std::endl;

    bool isUpdated = false, once = true;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    double runtime, inv_kb_temp = 1.0/lowTempStop;
    //int diffContacts = pModel->getSizeOfNeighborhood() - this->contactsPerBead;
    char addRemoveText[75];

    // contacts potential
//    double tempContactSum, currentContactsSum = contactPotential(&beads_in_use_tree, pModel);
//    double tempKLDivContacts, currentKLDivContacts = currentContactsSum*invTotal;
    //float target = 9.0f;

    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
    std::copy(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), tempContactsDistributionOfModel.begin());
    double tempKLDivContacts, currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

    const auto fifteenPercent = (unsigned int)(0.15*step_limit);
    const auto sixtyfivePercent = (unsigned int)(0.65*step_limit);
    const auto eightyfivePercent = (unsigned int)(0.85*step_limit);

    // slowly increase weight of total Contact energy over D_kl
    double currentCDW = alpha/1000;
    double alphaConstant = alpha > 0 ? std::pow(alpha/currentCDW, 1.0/(double)fifteenPercent) : 0;
    double contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
    double contactsAvg=3.0d;
    double currentLambda = lambda;

    double current_energy = (currentKL + currentLambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts);
    double this_energy, lowest_energy = current_energy;

    std::clock_t startTime;
    int attempts=0, failures=0, updated=0;

    std::vector<unsigned int> selections(totalBeadsInSphere);
    for(unsigned int i=0; i < workingLimit; i++){
        selections[i] = i; // some distances will exceed dmax
    }

    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomHull(0,hull.size()-1); // guaranteed unbiased

    unsigned int compCOunt=0, numberOfCoolingTempSteps=0;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        startTime = std::clock();
        auto set_it = hull.begin();

        if ( distribution(gen) < percentAddRemove) { //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
            // build a list of indices within defined region of convex hull
            if (distribution(gen) < 0.531 ){ // ADD BEAD?
                std::cout << "*******************                  ADD                   *******************" << std::endl;
                std::advance(set_it, randomHull(gen));
                unsigned int addMe = *set_it;
                itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), addMe);

                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);  // alters backUpState
                addToPr(addMe, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = pData->getScore(binCount);

                beads_in_use_tree.insert(addMe);
//                tempContactSum = addToContactsPotential(addMe, currentContactsSum, &beads_in_use_tree, pModel);
//                tempKLDivContacts=tempContactSum/(double)beads_in_use_tree.size();

                addToContactsDistribution(addMe, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                //populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

                tempNumberOfComponents = eulerTour.addNode(addMe, pModel);

                this_energy = (float)(testKL + currentLambda*connectivityPotential(tempNumberOfComponents) + contactDistributionWeight*tempKLDivContacts);

                if ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "   ADDED => %i", addMe);
                } else { // undo changes (rejecting)
                    beads_in_use_tree.erase(addMe);
                    eulerTour.removeNode(addMe);
                    std::sprintf(addRemoveText, "     ADD => %i FAILED", addMe);
                    auto beginBinCount = binCount.begin();
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }

            } else { // REMOVE BEADS?
                std::cout << "*******************                 REMOVE                 *******************" << std::endl;
                std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
                unsigned int original = bead_indices[selections[0]];

                for(const auto & select : selections){
                    tempNumberOfComponents = eulerTour.removeNode(original);
                    if (tempNumberOfComponents <= currentNumberOfComponents){
                        break;
                    }
                    eulerTour.addNode(original, pModel);
                    original = bead_indices[select];
                }

                removeLatticePositionToModel(bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);

                // still need to sort, swap changes the order
                testKL = pData->getScore(binCount);
                //tempContactSum = removeFromContactsPotential(original, currentContactsSum, &beads_in_use_tree, pModel);

                removeFromContactsDistribution(original, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                beads_in_use_tree.erase(original);
                //populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);
                //tempKLDivContacts=tempContactSum/(double)beads_in_use_tree.size();
                this_energy = testKL + currentLambda*connectivityPotential(tempNumberOfComponents) + contactDistributionWeight*tempKLDivContacts;

                if ( this_energy < current_energy || ( exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, " REMOVED => %i", original);
                } else { // undo changes and move to next bead (rejecting)
                    beads_in_use_tree.insert(original);
                    auto beginIt = bead_indices.begin();
                    auto beginBinCount = binCount.begin();
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);

                    eulerTour.addNode(original, pModel);
                    std::sprintf(addRemoveText, "  REMOVE => %i FAILED", original);
                }
            }

        } else { // positional refinement
            // only search within deadLimit, no need to recalculate at end
            attempts +=1;
            // shuffling should not change the location of the iterator
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);

            unsigned int position = selections[0];
            unsigned int swap1 = bead_indices[position];

            bool retry = true;

            if (numberOfCoolingTempSteps > eightyfivePercent){
                for(const auto & select : selections){
                    position = select;
                    swap1 = bead_indices[position];
                    if (numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1) < 2){
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            retry = false;
                            compCOunt++;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                }
            }

            if (retry){ // favors densly packed areas
                for(const auto & select : selections){
                    position =select;
                    swap1 = bead_indices[position];
                    if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                        break;
                    }
                    eulerTour.addNode(swap1, pModel);
                }
            }


            // remove contribution of swap1
            itIndex = bead_indices.begin() + position;
            // remove selected index from P(r)
            //std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
            removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

            // remove from beads_in_use_tree
//            tempContactSum = removeFromContactsPotential(swap1, currentContactsSum, &beads_in_use_tree, pModel);
            removeFromContactsDistribution(swap1, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            beads_in_use_tree.erase(swap1);

            unsigned int originalSwap2Value = swap1;

            // find a point in the neighborhood of swap1
//            if (distribution(gen) < 1.0 - acceptRate) { // try to only move beads that are singly connected
            if (false) { // try to only move beads that are singly connected
                bool retry=true;
                unsigned int neighborOfSwap1 = getConnectedNeighborFromSet(&beads_in_use_tree, pModel, swap1);

                for(unsigned int i=0; i<2; i++){
                    originalSwap2Value = getUseableNeighborFromSet(&beads_in_use_tree, pModel, neighborOfSwap1);
                    if (originalSwap2Value != noNeigborIndex && originalSwap2Value != swap1){
                        retry = false;
                        break;
                    }
                    neighborOfSwap1 = getConnectedNeighborFromSet(&beads_in_use_tree, pModel, swap1);
                }

                if (retry){
                    std::advance(set_it, randomHull(gen));
                    originalSwap2Value = *set_it;
                    while (originalSwap2Value == swap1){
                        set_it = hull.begin();
                        std::advance(set_it, randomHull(gen));
                        originalSwap2Value = *set_it;
                    }
                }
            } else {
                /*
                 * pick random point in HULL that doesn't equal swap1
                 */
                std::advance(set_it, randomHull(gen));
                originalSwap2Value = *set_it;
                while (originalSwap2Value == swap1){
                    set_it = hull.begin();
                    std::advance(set_it, randomHull(gen));
                    originalSwap2Value = *set_it;
                }
            }

            // make the swap, sort and update P(r)
            pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);
            std::iter_swap(itIndex, pSwap2);
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->getScore(binCount);

            beads_in_use_tree.insert(originalSwap2Value); // add new lattice
            addToContactsDistribution(originalSwap2Value, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
//            tempContactSum = addToContactsPotential(originalSwap2Value, tempContactSum, &beads_in_use_tree, pModel);
//            tempKLDivContacts=tempContactSum/(double)beads_in_use_tree.size();
//            populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

            this_energy = testKL + contactDistributionWeight*tempKLDivContacts;
            tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
            this_energy += currentLambda*connectivityPotential(tempNumberOfComponents);

            if (this_energy < current_energy || (std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                isUpdated = true;
                std::sprintf(addRemoveText, " SWAPPED => %i to %i ", swap1, originalSwap2Value);
            } else {
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                beads_in_use_tree.insert(swap1); // add back swap one to tree
                beads_in_use_tree.erase(originalSwap2Value);
                std::sprintf(addRemoveText, "  FAILED => %i", swap1);
                eulerTour.removeNode(originalSwap2Value);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
            }
        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        // Adaptive simulated annealing part
        if (isUpdated){

            currentKL = testKL;
            current_energy = this_energy;
            currentNumberOfComponents = tempNumberOfComponents;
            currentKLDivContacts = tempKLDivContacts;

            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

            populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
            randomHull = std::uniform_int_distribution<unsigned int>(0,(unsigned int)hull.size()-1);
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());

            unsigned int * const pSelections = &selections[0]; // initialized as empty in Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i;
            }

            if (numberOfCoolingTempSteps > sixtyfivePercent){
                status = once ? "ANNEALING" : status;
                if (once && numberOfCoolingTempSteps > eightyfivePercent && contactDistributionWeight*currentKLDivContacts/currentKL > alpha){
                    contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
                    current_energy = currentKL + currentLambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts;
                    lowest_energy = current_energy;
                    once = false;
                    status = "FINAL LOW TEMP EQ";
                }

                if (current_energy < lowest_energy){
                    lowestKL = currentKL;
                    lowest_energy = current_energy;
                }
            }

            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), contactsDistributionOfModel.begin());
            contactsAvg = 2.0*binCount[0]/(float)workingLimit;
        } else {
            std::copy(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), tempContactsDistributionOfModel.begin());
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        // rescale contactsWeight in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent){
            currentCDW *= alphaConstant;
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = (float)(currentKL + currentLambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts);
            lowest_energy = current_energy;
            status = "RAMPING";
        }

        double tempDKLValue = contactDistributionWeight*currentKLDivContacts;

        printf("      TEMP : %-.3E MAXSTEPS => %.0f (%i) %s \n", lowTempStop, step_limit, numberOfCoolingTempSteps, status.c_str());
        printf("    ACCEPT : %9.5f FAILURES => %i  TIME : %.4f\n", acceptRate, failures, runtime);
        printf("     LIMIT : %9i  COMP => %3i (%i) LAM : %.2E\n", workingLimit, currentNumberOfComponents, compCOunt, currentLambda);
        printf("  CONTACTS : %.3E (%.1E | %.1E) ALPHA : %.1E AVG %.2f\n", tempDKLValue, currentKLDivContacts, (tempDKLValue/currentKL), currentCDW, contactsAvg);
        printf("           : %s %s\n", pOutputFileName, addRemoveText);
        printf("   %s => %-4.3E ( %.4E ) ENRGY : %.4E ( %.4E )\n", score, currentKL, lowestKL, current_energy, lowest_energy);

        // random access
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pTempDuringRun[numberOfCoolingTempSteps] = contactsAvg;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKLDivContacts;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        //calculateModelPrDistributionDirect(&bead_indices, &testBinCount, workingLimit, pModel, pData);
        //calculateModelPrDistribution(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        float testKL1 = pData->getScore(testBinCount);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            std::cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << std::endl;
//            return "stopped";
//        }

        /*
         * check distribution
         */
//        populateContactsDistribution(testContactsDistributionOfModel, &beads_in_use_tree, pModel);
//        for(int c=0; c<13; c++){
//            if (testContactsDistributionOfModel[c] != contactsDistributionOfModel[c]){
//                std::cout << c << " contacts " << testContactsDistributionOfModel[c] << " != " << contactsDistributionOfModel[c] << std::endl;
//                exit(0);
//            }
//        }
        //checkSetAndVector(workingLimit, &bead_indices, &beads_in_use_tree);
        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        // update for running average
        // correlate the updating based on D_kl, which means only update if we get lower values
        // if we correlated with overall energy, then improvements in volume or contacts will cause constants to be re-updated
    } // end of steps

    //printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);

    unsigned int originalCount=0;
    for(auto sit = beads_in_use_tree.begin(); sit != beads_in_use_tree.end(); ++sit){
        if (original_hull.find(*sit) != original_hull.end()){
            originalCount++;
        }
    }

    std::cout << " ORIGINAL FROM HULL TOTAL : " << originalCount << " target count " << targetCount << std::endl;
//    std::string name = "ar_" + outputname;
//    pModel->writeModelToFile2(
//            currentKL,
//            workingLimit,
//            bead_indices,
//            binCount,
//            name,
//            this,
//            pData,
//            numberOfCoolingTempSteps,
//            10000,
//            0);
//    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
//    currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

    int totalCounts=0;
    for (unsigned int c=1; c<13; c++){
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalCounts++;
            }
        }
    }

    std::cout << " DISTRIBUTION OF CONTACTS" << std::endl;
    for (int c=1; c<13; c++){
        int totalContactsAt = 0;
        for (unsigned int i=0; i<workingLimit; i++){
            //contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, contactCutOff, pModel, pDistance);
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/(double)totalCounts);
    }

//    std::cout << "KL DIVERGENCE " << std::endl;
//    pData->printKLDivergence(binCount);

    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    pModel->setCVXHullVolume(current_volume);
    // reset temp for low temp annealing (strictly positional)
    float tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }
//
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    pModel->setAverageNumberOfContactsInModel(average_number_of_contacts);

    // CONSTANT TEMP REFINEMENT
    // make Contact Potential a percentage of final energy
//    // At end of each temp, update a probability model for volume?  Use this to select
//    // string tempName = "rough";
//    // pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, tempName, this, pData);
//
//    pModel->writeModelToFile(workingLimit, bead_indices, "before_final", numberOfCoolingTempSteps);
//    cout << "------------------------------------------------------------------------------" << endl;
//    printf(" NUMBER OF STEPS %i\n", numberOfCoolingTempSteps);
//    printf(" LATTCE AVG => %.0f        STDEV => %.0f\n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
//    cout << "------------------------------------------------------------------------------" << endl;
//

    std::cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << std::endl;
    std::cout << "------------------------------------------------------------------------------" << std::endl;
    std::cout << " FINAL LOW TEMP REFINEMENT" << std::endl;
    std::cout << "              KL DIVERGENCE   ENERGY " << std::endl;
    std::printf("      START =>    %.4E           \n", startingKL);
    std::printf("      FINAL =>    %.4E           %.4E  \n", currentKL, current_energy);
    std::cout << "    CONTACTS" << std::endl;
    std::printf("      AVERAGE =>    %.2f TOTAL => %.0f \n", average_number_of_contacts, tempAverageContacts);
    std::printf("    BIN [ONE] => %d AVG => %.3f \n", binCount[0], binCount[0]/(float)workingLimit);
    std::printf("    VOLUME  START => %.4E FINAL => %.4E \n", starting_volume, current_volume);
    std::cout << "------------------------------------------------------------------------------" << std::endl;
//
//   this->printContactList(bead_indices, &beads_in_use_tree, workingLimit, pModel);

//   //pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
//
//    // if multithread, must put a lock on these two step
//    //CVX HULL STUFF

    /*
     * print distances
     */
//    for(unsigned int i=0; i<workingLimit; i++){
//
//        unsigned next = i+1;
//        const vector3 * vec1 = &pModel->getBead(bead_indices[i])->getVec();
//
//        for(; next<workingLimit; next++){
//            auto vec3 = *vec1 - pModel->getBead(bead_indices[next])->getVec();
//            std::cout << i << " " << next << " " << vec3.length() << " " << pData->convertToBin(vec3.length()) << std::endl;
//        }
//    }


    //float average_number_of_contacts = 4.1;
    std::string nameOfModel = pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            outputname,
            this,
            pData,
            numberOfCoolingTempSteps,
            current_volume,
            average_number_of_contacts);

    return nameOfModel;
}


/**
 * Refinement only uses the lattice specified in the input seed file, no other lattice points are considered
 * Input should be an aligned set of models like from damaver or kde
 *
 * @param pModel
 * @param pData
 * @param outputname
 * @return
 */
std::string Anneal::refineHomogenousBodyInputEx(Model *pModel, Data *pData, std::string outputname) {
    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;

    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.00001 : (double)highTempStartForCooling;
    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    const char * pOutputFileName = outputname.c_str() ;
    std::string status = "REFINEMENT";

    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");

    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> active_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    std::vector<unsigned int>::iterator itIndex, pSwap2;

    /*
     * copy model into bead_indices
     */
    pModel->copyStartingModelIntoVector(bead_indices);

    unsigned int workingLimit = pModel->getStartingWorkingLimit();
    const auto addLimit = (unsigned int)(workingLimit*0.23f);
    pModel->writeSubModelToFile(0, totalBeadsInSphere, bead_indices, "universe");
    // reset iterators to internal bead_indices
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> restricted_set(bead_indices.begin(), bead_indices.begin() + workingLimit); // only these beads are allowable
    std::set<unsigned int> hull; // only these beads are allowable
    // set deadLimit of the selected set
    // convert distances in Search Sphere to ShannonBin membership
    this->fillPrBinsAndAssignTotalBin(pModel, pData);
    unsigned short int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    /*
     * create cross-validation set or working observed probability distribution that encompasses search space
     */
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);
    std::vector<double> testContactsDistributionOfModel(13);

    std::cout << "      TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "                BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "             BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;
    std::cout << "            WORKINGLIMIT : " << workingLimit << std::endl;

    calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    double testKL, currentKL = pData->getScore(binCount);
    double startingKL = currentKL;
    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;
    // CVX hull calculation
    char flags[] = "qhull FA";
    float starting_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    /*
     * model for initial search is too dense, so deadlimit layer should be sufficient for a confined search
     */
    // coupon collector's problem
    auto coupons = (unsigned int)(workingLimit*std::log((double)workingLimit) + 0.5772156649*workingLimit + 0.5d);

    unsigned int updateCount = ccmultiple*coupons;

    float step_limit = updateCount/0.65f;

    // if renormalization is too frequent, then a arrangemetn that may be too out of norm will be reset to low
    // so, need to span few steps?
    //
    // if 0.2, essentially divides step_limit into 5 steps for recalibrating weight
    // if 0.1 => 10 steps
    // if 0.05 => 20 steps
    std::vector<float> acceptanceRateDuringRun((unsigned long int)step_limit);
    std::vector<double> tempDuringRun((unsigned long int)step_limit);
    std::vector<float> divergenceDuringRun((unsigned long int)step_limit);
    std::vector<unsigned int> workingLimitDuringRun((unsigned long int)step_limit);

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents = eulerTour.getNumberOfComponents();
    unsigned int tempNumberOfComponents = currentNumberOfComponents;

    std::cout << "             EULER TOURS : " << currentNumberOfComponents << std::endl;

    if (currentNumberOfComponents > 1){
        eulerTour.printTourInfo();
        std::cout << " TOO MANY TOURS, RERUN" << std::endl;
        exit(0);
    }

    bool isUpdated = false, once = true;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    double runtime, inv_kb_temp = 1.0/lowTempStop;
    char addRemoveText[75];

    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
    std::copy(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), tempContactsDistributionOfModel.begin());
    //double tempKLDivContacts, currentKLDivContacts = contactPotential(&beads_in_use_tree, pModel);
    double tempKLDivContacts = 0, currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);
    double contactsAvg=3.0d;

    auto fifteenPercent = (unsigned int)(0.15*step_limit);
    auto eightyfivePercent = (unsigned int)(0.85*step_limit);

    // slowly increase weight of total Contact energy over D_kl
    double currentCDW = alpha/1000;
    double alphaConstant = alpha > 0 ? std::pow(alpha/currentCDW, 1.0/(double)fifteenPercent) : 0;
    double contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);

    double this_energy, current_energy = (currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts);
    double lowest_energy = current_energy;

    std::clock_t startTime;
    int attempts=0, failures=0, animateCount=0, updated=0, successess=0;

    std::vector<unsigned int> selections(totalBeadsInSphere);
    for(unsigned int i=0; i < workingLimit; i++){
        selections[i] = i; // some distances will exceed dmax
    }

    unsigned int hullsize = hull.size();
    unsigned int hullsizeLimit = hull.size()*0.1;
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomHull(0,hull.size()-1);   // guaranteed unbiased

    unsigned int numberOfCoolingTempSteps=0, compCOunt=0;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++){

        startTime = std::clock();

        // if hullsize is less than 10, need to remove beads
        if (numberOfCoolingTempSteps < addLimit || distribution(gen) < percentAddRemove ) { //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
            // build a list of indices within defined region of convex hull
            if (hullsize > hullsizeLimit && numberOfCoolingTempSteps > addLimit && distribution(gen) < 0.5317){ // ADD BEAD?
                std::cout << "*******************                  ADD                   *******************" << std::endl;

                auto set_it = hull.begin();
                std::advance(set_it, randomHull(gen));
                unsigned int original = *set_it;

                if (currentNumberOfComponents > 1) { // only add to existing bead
                    original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    while (original == noNeigborIndex) { // find a new position
                        original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    }
                }

                itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);

                // make the swap at the border (workingLimit)
                addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);  // alters backUpState
                addToPr(original, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                testKL = pData->getScore(binCount);

                beads_in_use_tree.insert(original);
                //populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                addToContactsDistribution(original, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                //tempKLDivContacts=contactPotential(&beads_in_use_tree, pModel);
                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

                tempNumberOfComponents=eulerTour.addNode(original, pModel);

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + contactDistributionWeight*tempKLDivContacts;

                if ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     ADD => %i", 1);
                    hull.erase(original);
                } else { // undo changes (rejecting)
                    beads_in_use_tree.erase(original);
                    eulerTour.removeNode(original);
                    std::sprintf(addRemoveText, "     ADD => %i", 0);
                    auto beginBinCount = binCount.begin();
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }

            } else { // REMOVE BEADS?
                std::cout << "*******************                 REMOVE                 *******************" << std::endl;
                // test for deletion
                std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
                unsigned int original = bead_indices[selections[0]];

                for(const auto & select : selections){
                    tempNumberOfComponents = eulerTour.removeNode(original);
                    if (tempNumberOfComponents <= currentNumberOfComponents){
                        break;
                    }
                    eulerTour.addNode(original, pModel);
                    original = bead_indices[select];
                }

                removeLatticePositionToModel(bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &original);

                // still need to sort, swap changes the order
                testKL = pData->getScore(binCount);

                removeFromContactsDistribution(original, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                beads_in_use_tree.erase(original);

                //populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);
                //tempKLDivContacts=contactPotential(&beads_in_use_tree, pModel);

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + contactDistributionWeight*tempKLDivContacts;

                if (this_energy < current_energy || (exp((current_energy - this_energy)*inv_kb_temp) > distribution(gen)) ) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     REMOVE => %i", original);
                    hull.insert(original);
                } else { // undo changes and move to next bead (rejecting)
                    std::sprintf(addRemoveText, "     REMOVE => %i FAILED", original);
                    beads_in_use_tree.insert(original);
                    auto beginIt = bead_indices.begin();
                    auto beginBinCount = binCount.begin();
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);

                    eulerTour.addNode(original, pModel);
                }
            }

//            calculateModelPrDistribution(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//            float testKL1 = pData->getScore(testBinCount);
//            if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//                std::cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << std::endl;
//                return "stopped";
//            }

        } else { // positional refinement
            // only search within deadLimit, no need to recalculate at end
            attempts +=1;
            // shuffling should not change the location of the iterator
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);

            unsigned int position = selections[0];
            unsigned int swap1 = bead_indices[position];

            bool retry = true;

            if (numberOfCoolingTempSteps > eightyfivePercent){
                for(const auto & select : selections){
                    position = select;
                    swap1 = bead_indices[position];
                    if (numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1) < 2){
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            compCOunt++;
                            retry = false;
                            break;
                        }
                        eulerTour.addNode(swap1, pModel);
                    }
                }
            }

            if (retry){ // favors densly packed areas
                for(const auto & select : selections){
                    position =select;
                    swap1 = bead_indices[position];
                    if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                        break;
                    }
                    eulerTour.addNode(swap1, pModel);
                }
            }

            // remove contribution of swap1
            itIndex = bead_indices.begin() + position;
            // remove selected index from P(r)
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
            removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // remove contribution to Contacts Distribution
            removeFromContactsDistribution(swap1, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            beads_in_use_tree.erase(swap1); // remove from in_use tree
            // find bead to move to
            /*
             * pick random point in HULL that doesn't equal swap1
             */
            auto set_it = hull.begin();
            std::advance(set_it, randomHull(gen));
            unsigned int originalSwap2Value = *set_it;
            while (originalSwap2Value == swap1){
                set_it = hull.begin();
                std::advance(set_it, randomHull(gen));
                originalSwap2Value = *set_it;
            }

            // make the swap, sort and update P(r)
            pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);
            std::iter_swap(itIndex, pSwap2);
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            addToPr(originalSwap2Value, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->getScore(binCount);

            beads_in_use_tree.insert(originalSwap2Value); // add new lattice
            addToContactsDistribution(originalSwap2Value, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

            this_energy = testKL + contactDistributionWeight*tempKLDivContacts;

            tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
            this_energy += lambda*connectivityPotential(tempNumberOfComponents);

            if (this_energy < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     SWAPPED => %i to %i ", swap1, originalSwap2Value);
                hull.erase(originalSwap2Value);
                hull.insert(swap1);
            }  else {
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                beads_in_use_tree.insert(swap1); // add back swap one to tree
                beads_in_use_tree.erase(originalSwap2Value);

                std::sprintf(addRemoveText, "      FAILED => %i", swap1);
                eulerTour.removeNode(originalSwap2Value);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
            }

        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        // Adaptive simulated annealing part
        if (isUpdated){

            currentKL = testKL;
            current_energy = this_energy;
            currentNumberOfComponents = tempNumberOfComponents;
            currentKLDivContacts = tempKLDivContacts;

            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased
            randomHull = std::uniform_int_distribution<unsigned int>(0,(unsigned int)hull.size()-1);
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
            hullsize = hull.size();

            unsigned int * const pSelections = &selections[0]; // initialized as empty in Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i;
            }

            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), contactsDistributionOfModel.begin());
            contactsAvg = 2.0*binCount[0]/(float)workingLimit;
        } else {
            std::copy(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), tempContactsDistributionOfModel.begin());
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }


        printf("      TEMP : %-.2E MAXSTEPS => %.0f (%4i : %i) %s \n", lowTempStop, step_limit, numberOfCoolingTempSteps, addLimit,  status.c_str());
        printf("    ACCEPT : %.5f FAILURES => %i  TIME : %.4f\n", acceptRate, failures, runtime);
        printf("     LIMIT : %5i  COMP => %3i (%i) HULL : %i \n", workingLimit, currentNumberOfComponents, compCOunt, hullsize);
        printf("  CONTACTS : %.3E (%.2E) ALPHA : %.2E AVG %.2f\n", contactDistributionWeight*currentKLDivContacts, currentKLDivContacts, currentCDW, contactsAvg);
        printf("           : %s %s\n", pOutputFileName, addRemoveText);
        printf("   %s => %-4.3E ENRGY : %.4E ( %.4E )\n", score, currentKL, current_energy, lowest_energy);

        // rescale etaFactor in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent ){
            currentCDW *= alphaConstant;
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts;
            lowest_energy = current_energy;// + totalContactEnergy ;
        }

        if (once && numberOfCoolingTempSteps > eightyfivePercent && contactDistributionWeight*currentKLDivContacts/currentKL > alpha){
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts;
            lowest_energy = current_energy;
            once = false;
            status = "FINAL LOW TEMP EQ";
        }
        /*
         * check distribution
         */
//        populateContactsDistribution(testContactsDistributionOfModel, &beads_in_use_tree, pModel);
//        for(int c=0; c<13; c++){
//            if (testContactsDistributionOfModel[c] != contactsDistributionOfModel[c]){
//                std::cout << c << " contacts " << testContactsDistributionOfModel[c] << " != " << contactsDistributionOfModel[c] << std::endl;
//                exit(0);
//            }
//        }
//        checkSetAndVector(workingLimit, &bead_indices, &beads_in_use_tree);

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        // update for running average
        // correlate the updating based on D_kl, which means only update if we get lower values
        // if we correlated with overall energy, then improvements in volume or contacts will cause constants to be re-updated
    } // end of steps

    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);

//    std::string name = "ar_" + outputname;
//    pModel->writeModelToFile2(
//            currentKL,
//            workingLimit,
//            bead_indices,
//            binCount,
//            name,
//            this,
//            pData,
//            numberOfCoolingTempSteps,
//            10000,
//            0);

//    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
//    currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);


    int totalCounts=0;
    for (unsigned int c=1; c<13; c++){
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalCounts++;
            }
        }
    }

    std::cout << " DISTRIBUTION OF CONTACTS" << std::endl;
    for (int c=1; c<13; c++){
        int totalContactsAt = 0;
        for (unsigned int i=0; i<workingLimit; i++){
            //contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, contactCutOff, pModel, pDistance);
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/(double)totalCounts);
    }

//    std::cout << "KL DIVERGENCE " << std::endl;
//    pData->printKLDivergence(binCount);

    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    pModel->setCVXHullVolume(current_volume);
    // reset temp for low temp annealing (strictly positional)
    float tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }
//
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    pModel->setAverageNumberOfContactsInModel(average_number_of_contacts);

    // CONSTANT TEMP REFINEMENT
    // make Contact Potential a percentage of final energy
//    // At end of each temp, update a probability model for volume?  Use this to select
//    // string tempName = "rough";
//    // pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, tempName, this, pData);
//
//    pModel->writeModelToFile(workingLimit, bead_indices, "before_final", numberOfCoolingTempSteps);
//    cout << "------------------------------------------------------------------------------" << endl;
//    printf(" NUMBER OF STEPS %i\n", numberOfCoolingTempSteps);
//    printf(" LATTCE AVG => %.0f        STDEV => %.0f\n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
//    cout << "------------------------------------------------------------------------------" << endl;
//

    std::cout << " AVERAGE CONTACTS FROM SET " << average_number_of_contacts << std::endl;
    std::cout << "------------------------------------------------------------------------------" << std::endl;
    std::cout << " FINAL LOW TEMP REFINEMENT" << std::endl;
    std::cout << "              KL DIVERGENCE   ENERGY " << std::endl;
    std::printf("      START =>    %.4E           \n", startingKL);
    std::printf("      FINAL =>    %.4E           %.4E  \n", currentKL, current_energy);
    std::cout << "    CONTACTS" << std::endl;
    std::printf("      AVERAGE =>    %.2f TOTAL => %.0f \n", average_number_of_contacts, tempAverageContacts);
    std::printf("    BIN [ONE] => %d AVG => %.3f \n", binCount[0], binCount[0]/(float)workingLimit);
    std::printf("    VOLUME  START => %.4E FINAL => %.4E \n", starting_volume, current_volume);
    std::cout << "------------------------------------------------------------------------------" << std::endl;
//
//   this->printContactList(bead_indices, &beads_in_use_tree, workingLimit, pModel);

//   //pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
//
//    // if multithread, must put a lock on these two step
//    //CVX HULL STUFF

    //float average_number_of_contacts = 4.1;
    std::string nameOfModel = pModel->writeModelToFile2(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            outputname,
            this,
            pData,
            numberOfCoolingTempSteps,
            current_volume,
            average_number_of_contacts);

    return nameOfModel;
}

