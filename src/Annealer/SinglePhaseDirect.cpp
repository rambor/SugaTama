//
// Created by xos81802 on 20/07/2018.
//

#include "../Anneal.h"

bool Anneal::createInitialModelCVXHullDirect(Model *pModel, Data *pData, std::string name) {

    std::cout << " --------------------------------------------------- " << std::endl;
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    contactCutOff = interconnectivityCutOff;

    maxbin = pModel->getMaxBin(pData);
    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");
    const char * pOutputFileName = name.c_str() ;

    maxbin += 1;  // maximum bin index, so vector size must be +1
    std::cout << "         INITIAL MAX BIN : " << maxbin << std::endl;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::cout << "      TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "                BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "             BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> lowest_bead_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);
    //std::clock_t start;
    // c-style is slightly faster for large vector sizes
    // start = std::clock();
    unsigned int * const ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    // cout << "C-STYLE " << (std::clock() - start)/(double) CLOCKS_PER_SEC << endl;
    std::random_device rd;
    std::mt19937 gen(rd());

    lowerV = (unsigned int)(pData->getVolume()  - pData->getVolume() *0.39);
    upperV = (unsigned int)(pData->getVolume()  + pData->getVolume() *0.05);

    float radius_larger = pModel->getBeadRadius()*std::sqrt(7.0f)/2.0f;
    auto lowerN = (unsigned int)std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger));
    auto upperN = (unsigned int)std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger));

    // pick random number between lowerN and upperN
    unsigned int workingLimit = lowerN;//number_of_beads_to_use(gen);
    unsigned int lowestWorkingLimit = workingLimit;
    // create initial model
    // randomize bead indices
    // sort to workingLimit
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);

    std::set<unsigned int> hull;  // sort hull points into bead_indices
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());

    std::cout << "        WORKINGLIMIT SET : " << workingLimit << std::endl;
    // setup parameters for hull
    char flags[] = "qhull FA";
    recalculateDeadLimit(workingLimit, bead_indices, hull, pModel, totalBeadsInSphere);
    std::string switched = "HULL";
    coordT points[3*(upperN+(unsigned int)(0.10*upperN))];

    float test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    float lowest_volume = current_volume;
    float initialVolume = current_volume;

    std::cout << "              CVX VOLUME : " << current_volume  << std::endl;

    EulerTour eulerTour(bead_indices.cbegin(), workingLimit, pModel);

    unsigned int lowestNumberOfComponents=0, tempNumberOfComponents, currentNumberOfComponents  = eulerTour.getNumberOfComponents();
    std::cout << "      CREATED EULER TOUR : " << currentNumberOfComponents  << std::endl;

    calculateModelPrDistributionDirect(&bead_indices, &binCount, workingLimit, pModel, pData);
    float currentKL = pData->getScore(binCount);
    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy

    float hlambda = pData->getIsIntensity() ? 10.0f : 0.0001f;
//    float muConstant = pData->getIsPr() ? 0.000001f : 0.1f; // chi tends to start in 100's whereas DKL is 0.1;
    double muPercent = 0.1;
    double muConstant = currentKL*muPercent/((1.0 - muPercent)*current_volume/(double)workingLimit);


    float targetVolume = upperV;
    float diff = std::abs(current_volume-targetVolume)/targetVolume;
    float temp_volume_energy, current_volume_energy = muConstant*current_volume/(double)workingLimit;
    float testKL, test_energy, current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + current_volume_energy;//(float)workingLimit;

    std::cout << "           INITIAL SCORE : " << currentKL << std::endl;
    float lowest_energy = current_energy;

    auto lowTempStop =  (double)0.1;//highT;
    double inv_kb_temp = 1.0f/lowTempStop;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;
    unsigned int failures=0;

    double sumIt = 0.0, counter=0.0;
    bool modulo = false, updateCVX = false;
    //int sizeOfNeighborhood = pModel->getSizeOfNeighborhood();  seems on average there are less than 12 available to look at
    std::cout << " STARTING CONSTANT TEMP SEARCH " << currentNumberOfComponents << std::endl;

    unsigned int high, hullsize=hull.size();
    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomHull(0,hullsize); // guaranteed unbiased

    double probabilityAddRemove = 0.361;

    for (high=0; high < highTempRounds; high++) { // iterations during the high temp search

        if (distribution(gen) < probabilityAddRemove){

            if (distribution(gen) < -0.1 && currentNumberOfComponents > 1){ // scramble
                recalculateDeadLimit(workingLimit, bead_indices, hull, pModel, totalBeadsInSphere);
                hullsize = hull.size();
                randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);

                unsigned int limit = workingLimit;
                for(auto bit : hull){ // sort all beads in hull and combine with workinglimit
                    auto findIt = std::find(bead_indices.begin() + limit, bead_indices.end(), bit);
                    if (findIt != bead_indices.end()){
                        std::iter_swap(findIt, bead_indices.begin() + limit);
                        limit++;
                    } else {
                        exit(0);
                    }
                }
                unsigned int tempWorkingLimit = workingLimit; //number_of_beads_to_use(gen);

                std::shuffle(bead_indices.begin(), bead_indices.begin() + limit, gen);
                std::sort(bead_indices.begin(), bead_indices.begin() + tempWorkingLimit);

                calculateModelPrDistributionDirect(&bead_indices, &binCount, tempWorkingLimit, pModel, pData);
//                calculateModelPrDistribution(&bead_indices, &binCount, tempWorkingLimit, totalBeadsInSphere, pModel, pData);
                testKL = pData->getScore(binCount);

                test_volume = calculateCVXHULLVolume(flags, &bead_indices, tempWorkingLimit, pModel);
                //temp_volume_energy = muConstant*std::abs(test_volume-targetVolume)/targetVolume;
                temp_volume_energy = muConstant*test_volume/(double)tempWorkingLimit;

                EulerTour tempTour(bead_indices.cbegin(), tempWorkingLimit, pModel);
                tempNumberOfComponents=tempTour.getNumberOfComponents();
                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy;

                if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) {
                    beads_in_use_tree = std::set<unsigned int>(bead_indices.begin(), bead_indices.begin()+tempWorkingLimit);
                    workingLimit = tempWorkingLimit;
                    eulerTour = EulerTour(tempTour); // copy constructor called
                    updateCVX = true;
                } else { // undo changes (rejecting)
                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());   // make backup copy
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // make backup copy
                }

            } else {
                if ( (distribution(gen) < 0.5 && workingLimit < upperN) || workingLimit < lowerN ){
                    std::cout << "*******************                  ADD                   ******************* "  << std::endl;

                    std::set<unsigned int>::const_iterator set_it=hull.begin();
                    std::advance(set_it, randomHull(gen));
                    unsigned int testIndex = *set_it;

                    while ( beads_in_use_tree.find(testIndex) != beads_in_use_tree.end()) { // find a new position
                        set_it=hull.begin();
                        std::advance(set_it, randomHull(gen));
                        testIndex = *set_it;
                    }

                    auto itIndex = std::find(bead_indices.begin(), bead_indices.end(), testIndex);
                    // make the swap at the border (workingLimit)
                    addLatticPositionToModelDirect(&bead_indices, &workingLimit, &itIndex);
                    addToPrDirect(testIndex, bead_indices, workingLimit, binCount, pModel, pData);

                    testKL = pData->getScore(binCount);
                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

                    tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                    //temp_volume_energy = muConstant*std::abs(test_volume-targetVolume)/targetVolume;
                    temp_volume_energy = muConstant*test_volume/(double)workingLimit;
                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy;

                    if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                        beads_in_use_tree.insert(testIndex);
                        updateCVX = true;
                    } else { // undo changes (rejecting)
                        restoreAddingFromBackUpDirect(&workingLimit, &binCountBackUp, &binCount);
                        eulerTour.removeNode(testIndex);
                    }


                } else {

                    std::cout << "*******************                 REMOVE                 ******************* "  << std::endl;

                    // grab from randomized active_indices list
                    unsigned int testIndex = bead_indices[randomIndex(gen)];

                    removeLatticePositionToModelDirect(bead_indices, binCount, &workingLimit, &testIndex, pModel, pData);

                    testKL = pData->getScore(binCount);
                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

                    tempNumberOfComponents = eulerTour.removeNode(testIndex);

                    //temp_volume_energy = muConstant*std::abs(test_volume-targetVolume)/targetVolume;
                    temp_volume_energy = muConstant*test_volume/(double)workingLimit;
                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy;

                    if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                        beads_in_use_tree.erase(testIndex);
                        updateCVX = true;
                    } else { // undo changes and move to next bead (rejecting)
                        eulerTour.addNode(testIndex, pModel);
                        restoreRemovingLatticePointFromBackUpDirect(&workingLimit, &binCountBackUp, &binCount);
                    }
                }
            }


        } else {
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            // recalculate
            // calculate convex hull and get hull points
            // can be threaded

            unsigned int swap1 = bead_indices[randomIndex(gen)];

            if (currentNumberOfComponents > 1 && modulo){ // only move CVX points if true
                const unsigned int * pBead = bead_indices.data();
                std::vector<unsigned int> active_indices(workingLimit);
                for (unsigned int i = 0; i < workingLimit; i++) {
                    beadToPoint(&points[i*3], pModel->getBead(pBead[i]));
                    active_indices[i] = pBead[i];
                }

                // needs to be optimized
                qh_new_qhull(3, workingLimit, points, 0, flags, nullptr, nullptr);
                vertexT * vertices = qh vertex_list;
                auto totalV = (unsigned int)qh num_vertices;

                // only move CVX hull points
                std::vector<unsigned int> indices_to_check(totalV);
                for (unsigned int v = 0; v < totalV; v++) { //
                    indices_to_check[v] = active_indices[qh_pointid( vertices->point)];
                    vertices = vertices->next;
                }

                qh_freeqhull(true);
                std::uniform_int_distribution<unsigned int> randomVertices(0, totalV-1);
                swap1 = indices_to_check[randomVertices(gen)];
            }


            // find bead to swap in active set
            auto itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
            // remove selected index from P(r)
            removeFromPrDirect(swap1, bead_indices, workingLimit, binCount, pModel, pData);
            // find better position
            eulerTour.removeNode(swap1);

            // swap to anywhere in the CVX HULL
            auto set_it=hull.begin();
            std::advance(set_it, randomHull(gen));
            unsigned int neighbor = *set_it;

            while (beads_in_use_tree.find(neighbor) != beads_in_use_tree.end()) { // find a new position
                set_it=hull.begin();
                std::advance(set_it, randomHull(gen));
                neighbor = *set_it;
            }

            auto pSwap2 = std::find(bead_indices.begin(), bead_indices.end(), neighbor); // potential problem
            std::iter_swap(itIndex, pSwap2);

            addToPrDirect(neighbor, bead_indices, workingLimit, binCount, pModel, pData);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->getScore(binCount);
            tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

            test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

            //temp_volume_energy = muConstant*std::abs(test_volume-targetVolume)/targetVolume;
            temp_volume_energy = muConstant*test_volume/(double)workingLimit;
            test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy;

            if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen)) ) {
                beads_in_use_tree.erase(swap1);
                beads_in_use_tree.insert(neighbor);
                updateCVX = true;
            } else {
                std::iter_swap(pSwap2, itIndex);
                eulerTour.removeNode(neighbor);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
            }
        }


        if (updateCVX){
            current_energy = test_energy;
            currentKL = testKL;
            current_volume = test_volume;
            current_volume_energy = temp_volume_energy;
            currentNumberOfComponents = tempNumberOfComponents;

            acceptRate = inv500slash499*acceptRate+inv500;
            updateCVX = false;
            failures=0;

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy

//            if (current_volume < targetVolume){
                populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
                switched = "NEAREST NEIGHBOR";
//            } else {
//                switched = "HULL";
//                recalculateDeadLimit(workingLimit, bead_indices, hull, pModel, totalBeadsInSphere);
//            }

            hullsize = hull.size();
            randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);
            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased

            if (current_volume_energy > 1.5*currentKL && currentNumberOfComponents == 1){
                muConstant = currentKL*muPercent/((1.0 - muPercent)*current_volume/(double)workingLimit);
                current_volume_energy = muConstant*current_volume/(double)workingLimit;
            }

        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASAConstantTemp(high, highTempRounds, acceptRate, lowTempStop, inv_kb_temp);
        std::cout << "*******************                                        *******************" << std::endl;
        printf("       ACCEPTRATE : %5.3f       STEP : %-7i (%7i) DIRECT\n", acceptRate, high, highTempRounds);
        printf("            GRAPH : %5i   FAILURES : %i :%s: %s\n", currentNumberOfComponents, failures, pOutputFileName, switched.c_str());
        printf("            LIMIT : %5i   UPPER >= %i LOWER <= %i POOL : %i\n", workingLimit, upperN, lowerN, hullsize);
        printf("           VOLUME : %.0f (%.0f) MU : %.2E \n", current_volume, targetVolume, muConstant);
        printf("         %s : %.4E E_VOL : %.2E ENRGY : %.4E \n", score, currentKL, current_volume_energy, current_energy);

        modulo = (high%50) ? true : modulo;

//        std::cout << "*******************                                        *******************" << std::endl;
//        printf("       ACCEPTRATE : %5.3f       STEP : %-7i (%7i) DIRECT\n", acceptRate, high, highTempRounds);
//        printf("            GRAPH : %5i   FAILURES : %i %s\n", currentNumberOfComponents, failures, switched.c_str());
//        printf("            LIMIT : %5i   UPPER >= %i LOWER <= %i POOL : %i\n", workingLimit, upperN, lowerN, hullsize);
//        printf("           VOLUME : %.0f (%.0f) MU : %.2E \n", current_volume, targetVolume, muConstant);
//        printf("         %s : %.4E E_VOL : %.2E ENRGY : %.4E \n", score, currentKL, current_volume_energy, current_energy);
//        std::cout << "*******************                                        *******************" << std::endl;


        // check Pr values
        // uncomment to check update of Pr and direct calc are equal
//        calculateModelPrDistributionDirect(&bead_indices, &testBinCount, workingLimit, pModel, pData);
//        float testKL1 = pData->getScore(testBinCount);
//
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            printf("               D_KL : %.8E : %.8E \n", currentKL, testKL1);
//            return "stopped";
//        }
//
        // test bead_indices and vector are in sync
//        for(int i=0; i<workingLimit; i++){
//            if (beads_in_use_tree.find(bead_indices[i]) == beads_in_use_tree.end()){
//                std::cout << "!!!!! [ " << i << " ] NOT FOUND IN TREE => " << bead_indices[i] << std::endl;
//                exit(0);
//            }
//        }
//
//        for(auto bit = beads_in_use_tree.begin(); bit != beads_in_use_tree.end(); ++bit){
//            if(std::find(bead_indices.begin(), bead_indices.begin()+workingLimit, *bit) == bead_indices.begin()+workingLimit){
//                std::cout << "!!!!! " << " NOT FOUND in BEAD_INDICES " << *bit << std::endl;
//                exit(0);
//            }
//        }
//
//        if(beads_in_use_tree.size() != workingLimit){
//            cout << "!!!!! " << " Incorrect size WORKINGLIMIT "  << endl;
//            exit(0);
//        }

        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            lowestNumberOfComponents=1;
            lowest_volume = current_volume;
            lowestWorkingLimit = workingLimit;
            sumIt += workingLimit;
            counter += 1.0;
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
//        if (checkForRepeats(lowest_bead_indices)){
//            cout << " STOPPED REPEATS found " << endl;
//            return "stopped";
//        }
            // update D_kl
//            muConstant = mu*currentKL/((1.0f-mu)*current_volume);
//            current_energy = currentKL + muConstant*current_volume;
            lowest_energy = current_energy;
            //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_HT_" +  std::to_string(high));
        }
    }

    EulerTour tempEulerTour(lowest_bead_indices.begin(), lowestWorkingLimit, pModel);
    std::cout << " Lowest Tour : " << tempEulerTour.getNumberOfComponents() << std::endl;
    pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name, high);

    std::set<unsigned int> centeringhullpts;
    std::set<unsigned int> lowestSet(lowest_bead_indices.begin(), lowest_bead_indices.begin()+lowestWorkingLimit);
    getHullPoints(centeringhullpts, lowestSet, pModel);

    pModel->centerLatticeModel(&lowestWorkingLimit, lowest_bead_indices, centeringhullpts);


    tempEulerTour = EulerTour(lowest_bead_indices.begin(), lowestWorkingLimit, pModel);
    std::cout << " Lowest Tour AFTER  : " << tempEulerTour.getNumberOfComponents() << std::endl;

    pModel->setStartingSet(lowest_bead_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
    pModel->setBeadAverageAndStdev(lowestWorkingLimit, lowestWorkingLimit*0.2f);

    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    printf("          AVERAGE => %.0f (%.0f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    printf("      TRK AVERAGE => %.2f (%i) \n", (sumIt/counter), (int)counter);
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************                 VOLUME                 *******************" << std::endl;
    printf("          INITIAL => %.0f \n", initialVolume);
    printf("           LOWEST => %.0f \n", lowest_volume);
    printf("            FINAL => %.0f\n", current_volume);
    std::cout << "*******************                                        *******************" << std::endl;


    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "    SHANNON BINS IN DATA : " << pData->getShannonBins() << std::endl;

    if (pData->getIsPr()){
        std::cout << "            ZERO BIN AT : " << pData->getZeroBin() << std::endl;
    }


    std::cout << " EXITING INITIAL MODEL BUILD => " << currentNumberOfComponents << std::endl;
    float tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }
//
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    std::cout << " AVERAGE NUMBER CONTACTS : " << average_number_of_contacts << std::endl;

    if (lowestNumberOfComponents == 1 ) {
        return true;
    } else {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return false;
    }

}


std::string Anneal::refineHomogenousBodyASAHybridDirect(Model *pModel, Data *pData, std::string outputname) {

    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "#        STARTING ASA REFINEMENT OF HOMOGENOUS BODY             # " << std::endl;
    std::cout << "#                                                               # " << std::endl;
    std::cout << "########################<<<<<>>>>>>############################## " << std::endl;

    double lowTempStop = (asaAcceptanceRate < 0.11) ? (double) highTempStartForCooling * 0.000000001
                                                    : (double) highTempStartForCooling;
    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    const char * pOutputFileName = outputname.c_str() ;
    std::string status = "RAMPING TO TEMP";

    std::random_device rd;
    std::mt19937 gen(rd());
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0, 1.0);

    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> active_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    // copy Starting_Set from initial model
    pModel->copyStartingModelIntoVector(bead_indices);

    unsigned int workingLimit = pModel->getStartingWorkingLimit();
    //const auto lowerLimit = (unsigned int) (workingLimit / 2.31f);

    // reset iterators to internal bead_indices
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> hull;
    populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
    // set deadLimit of the selected set
    // convert distances in Search Sphere to ShannonBin membership
    this->fillPrBinsAndAssignTotalBinDirect(pModel, pData);

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

    std::cout << "      TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "                BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "             BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    calculateModelPrDistributionDirect(&bead_indices, &binCount, workingLimit, pModel, pData);
    float testKL, currentKL = pData->getScore(binCount);
    float lowestKL = currentKL, startingKL = currentKL;

    // CVX hull calculation
    char flags[] = "qhull FA";
    float starting_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    // Surface area calculations
    std::vector<float> weights(workingLimit);
    //float deltar = pModel->getBeadRadius()*0.74;
    //float deltar = pModel->getBeadRadius()*sqrt(3.0)/3.0;
    std::fill(weights.begin(), weights.end(), pModel->getBeadRadius() + delta_r);
    std::vector<Eigen::Vector3f> coordinates(workingLimit);

    for (unsigned int i = 0; i < workingLimit; i++) {
        vector3 pBeadVec = pModel->getBead(bead_indices[i])->getVec();
        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
    }

    const auto sa = (float) (4.0f * M_PI * (pModel->getBeadRadius() + delta_r) * (pModel->getBeadRadius() + delta_r));
    float currentStoV = surfaceToVolume(workingLimit, weights, coordinates);
    // end surface area calculations

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
    /*
     * model for initial search is too dense, so deadlimit layer should be sufficient for a confined search
     */
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;

    // coupon collector's problem
    auto coupons = (unsigned int) (workingLimit * std::log((double) workingLimit) + 0.5772156649 * workingLimit + 0.5d);

    const unsigned int updateCount = ccmultiple * coupons;
    unsigned int original, updatedCount = 0;

    float step_limit = updateCount / 0.65f;///0.1;

    // if renormalization is too frequent, then a arrangemetn that may be too out of norm will be reset to low
    // so, need to span few steps?
    //
    // if 0.2, essentially divides step_limit into 5 steps for recalibrating weight
    // if 0.1 => 10 steps
    // if 0.05 => 20 steps
    std::vector<float> acceptanceRateDuringRun((unsigned long int) step_limit);
    std::vector<double> tempDuringRun((unsigned long int) step_limit);
    std::vector<float> divergenceDuringRun((unsigned long int) step_limit);
    std::vector<unsigned int> workingLimitDuringRun((unsigned long int) step_limit);

    double *pTempDuringRun = &tempDuringRun.front();
    float *pAcceptanceRateDuringRun = &acceptanceRateDuringRun.front();
    float *pDivergenceDuringRun = &divergenceDuringRun.front();
    unsigned int *pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents = eulerTour.getNumberOfComponents();
    unsigned int tempNumberOfComponents = currentNumberOfComponents;

    std::cout << "             EULER TOURS : " << currentNumberOfComponents << std::endl;

//    if (currentNumberOfComponents > 1) {
//        std::cout << " TOO MANY TOURS, RERUN" << std::endl;
//        exit(0);
//    }

    bool isUpdated = false, once=true;
    float acceptRate = 0.5f, inv500 = 1.0f / 500.0f;
    float inv500slash499 = 499.0f / 500.0f;

    double runtime, inv_kb_temp = 1.0 / lowTempStop;
    //int diffContacts = pModel->getSizeOfNeighborhood() - this->contactsPerBead;
    float this_energy;
    char addRemoveText[75];

    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
    std::copy(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), tempContactsDistributionOfModel.begin());
    double tempKLDivContacts = 0, currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

    const auto fifteenPercent = (unsigned int)(0.15*step_limit);
    const auto sixtyfivePercent = (unsigned int)(0.65*step_limit);
    const auto eightyfivePercent = (unsigned int)(0.85*step_limit);

    // slowly increase weight of total Contact energy over D_kl
    double currentCDW = alpha/1000;
    double alphaConstant = alpha > 0 ? std::pow(alpha/currentCDW, 1.0/(double)fifteenPercent) : 0;
    double contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
    double contactsAvg=3.0d;

    auto current_energy = (float) (currentKL + lambda * connectivityPotential(currentNumberOfComponents) + contactDistributionWeight * currentKLDivContacts);
    float lowest_energy = current_energy;

    std::clock_t startTime;
    int attempts = 0, failures = 0, animateCount = 0, updated = 0, successess = 0;

    std::vector<unsigned int> selections(totalBeadsInSphere);
    for (unsigned int i = 0; i < workingLimit; i++) {
        selections[i] = i; // some distances will exceed dmax
    }

    std::uniform_int_distribution<unsigned int> randomHull(0,hull.size()-1); // guaranteed unbiased
    unsigned int compCOunt=0, numberOfCoolingTempSteps=0;


    for (; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {

        auto set_it = hull.begin();
        startTime = std::clock();

        if (distribution(gen) < percentAddRemove) { //add or remove bead within working Set (exclude deadzone)
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            // add remove based on lower and upper limits
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
            // build a list of indices within defined region of convex hull
            if (distribution(gen) < 0.531) { // ADD BEAD?
                std::cout << "*******************                  ADD                   *******************" << std::endl;
                std::advance(set_it, randomHull(gen));
                original = *set_it;
                auto itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);

                // make the swap at the border (workingLimit)
                addLatticPositionToModelDirect(&bead_indices, &workingLimit, &itIndex);
                addToPrDirect(original, bead_indices, workingLimit, binCount, pModel, pData);

                testKL = pData->getScore(binCount);

                beads_in_use_tree.insert(original);
                addToContactsDistribution(original, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

                tempNumberOfComponents = eulerTour.addNode(original, pModel);

                this_energy = testKL + lambda * connectivityPotential(tempNumberOfComponents) + contactDistributionWeight * tempKLDivContacts;

                if (this_energy < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                    isUpdated = true;
                    sprintf(addRemoveText, "   ADDED => %i", original);
                } else { // undo changes (rejecting)
                    beads_in_use_tree.erase(original);
                    eulerTour.removeNode(original);
                    sprintf(addRemoveText, "     ADD => %i FAILED", original);
                    restoreAddingFromBackUpDirect(&workingLimit, &binCountBackUp, &binCount);
                }
            } else { // REMOVE BEADS?
                // test for deletion, only delete those that maintain or improve currentNumberOfComponents
                std::cout << "*******************                 REMOVE                 *******************" << std::endl;
                std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
                original = bead_indices[selections[0]];

                for(const auto & select : selections){
                    tempNumberOfComponents = eulerTour.removeNode(original);
                    if (tempNumberOfComponents <= currentNumberOfComponents){
                        break;
                    }
                    eulerTour.addNode(original, pModel);
                    original = bead_indices[select];
                }

                removeLatticePositionToModelDirect(bead_indices, binCount, &workingLimit, &original, pModel, pData);
                // still need to sort, swap changes the order
                testKL = pData->getScore(binCount);

                removeFromContactsDistribution(original, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                beads_in_use_tree.erase(original);

                //populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

                this_energy = testKL + lambda * connectivityPotential(tempNumberOfComponents) +
                        contactDistributionWeight * tempKLDivContacts;

                if (this_energy < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                    sprintf(addRemoveText, "    REMOVED => %i", original);
                    isUpdated = true;
                } else { // undo changes and move to next bead (rejecting)
                    sprintf(addRemoveText, "     REMOVE => %i FAILED", original);
                    beads_in_use_tree.insert(original);
                    restoreRemovingLatticePointFromBackUpDirect(&workingLimit, &binCountBackUp, &binCount);
                    eulerTour.addNode(original, pModel);
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
            auto itIndex = bead_indices.begin() + position;
            // remove selected index from P(r)
            //std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
            removeFromPrDirect(swap1, bead_indices, workingLimit, binCount, pModel, pData);
            // remove contribution to Contacts Distribution
            removeFromContactsDistribution(swap1, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            beads_in_use_tree.erase(swap1); // remove from in_use tree
            // find bead to move to
            std::advance(set_it, randomHull(gen));
            unsigned int originalSwap2Value = *set_it;
            while (originalSwap2Value == swap1){
                set_it = hull.begin();
                std::advance(set_it, randomHull(gen));
                originalSwap2Value = *set_it;
            }

            // make the swap, sort and update P(r)
            auto pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);
            std::iter_swap(itIndex, pSwap2);

//                std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            addToPrDirect(originalSwap2Value, bead_indices, workingLimit, binCount, pModel, pData);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->getScore(binCount);
            beads_in_use_tree.insert(originalSwap2Value); // add new lattice

            addToContactsDistribution(originalSwap2Value, tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            //populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

            this_energy = testKL  + contactDistributionWeight*tempKLDivContacts;
            tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
            this_energy += lambda*connectivityPotential(tempNumberOfComponents);

            if (this_energy < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen))) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     SWAPPED => %i to %i ", swap1, originalSwap2Value);
            } else {
                std::iter_swap(pSwap2, itIndex);
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                beads_in_use_tree.insert(swap1); // add back swap one to tree
                beads_in_use_tree.erase(originalSwap2Value);

                std::sprintf(addRemoveText, "      FAILED => %i", swap1);
                eulerTour.removeNode(originalSwap2Value);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
            }

        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

        printf("      TEMP : %-.3E MAXSTEPS => %.0f (%i) %s \n", lowTempStop, step_limit, numberOfCoolingTempSteps, status.c_str());
        printf("    ACCEPT : %9.5f FAILURES => %i  TIME : %.4f\n", acceptRate, failures, runtime);
        printf("     LIMIT : %9i  COMP => %3i (%i) \n", workingLimit, currentNumberOfComponents, compCOunt);
        printf("  CONTACTS : %.3E (%.2E) ALPHA : %.2E AVG %.2f\n", contactDistributionWeight*currentKLDivContacts, currentKLDivContacts, currentCDW, contactsAvg);
        printf("    DIRECT : %s %s\n", pOutputFileName, addRemoveText);
        printf("   %s => %-4.3E ( %.4E ) ENRGY : %.4E ( %.4E )\n", score, currentKL, lowestKL, current_energy, lowest_energy);

        // random access
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pTempDuringRun[numberOfCoolingTempSteps] = currentStoV;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

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
                    current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts;
                    lowest_energy = current_energy;
                    once = false;
                    status = "FINAL LOW TEMP EQ";
                }
            }

            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), contactsDistributionOfModel.begin());
            contactsAvg = 2.0*binCount[0]/(float)workingLimit;
        } else {
            std::copy(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), tempContactsDistributionOfModel.begin());
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        /*
         * check distribution
         */
//        populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
//        for(int c=0; c<13; c++){
//            if (tempContactsDistributionOfModel[c] != contactsDistributionOfModel[c]){
//                std::cout << c << " contacts " << tempContactsDistributionOfModel[c] << " != " << contactsDistributionOfModel[c] << std::endl;
//                exit(0);
//            }
//        }
//        checkSetAndVector(workingLimit, &bead_indices, &beads_in_use_tree);

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);
        // update for running average
        // correlate the updating based on D_kl, which means only update if we get lower values
        // if we correlated with overall energy, then improvements in volume or contacts will cause constants to be re-updated
        // rescale etaFactor in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent){
            currentCDW *= alphaConstant;
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = (float)(currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts);
            lowest_energy = current_energy;
            status = "RAMPING";
        }

    }

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

    weights.resize(workingLimit);
    std::fill(weights.begin(), weights.end(), pModel->getBeadRadius() + delta_r);
    coordinates.resize(workingLimit);

    for (unsigned int i = 0; i < workingLimit; i++) {
        vector3 pBeadVec = pModel->getBead(bead_indices[i])->getVec();
        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
    }
    float finalStoV = surfaceToVolume(workingLimit, weights, coordinates);


    int totalCounts=0;
    for (int c=1; c<13; c++){
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

    std::cout << "KL DIVERGENCE " << std::endl;
    pData->printKLDivergence(binCount);

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

    outputname += "_direct";
std::cout << " DIRECT " << std::endl;

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
 * Random add/remove to generate initial model for seeded modeling
 *
 * Try to find the minimial set of lattice points that agree with atomistic P(r)
 * Uses Number of Contacts per bead in objective function, set eta to 0 for none
 *
 */
bool Anneal::createSeedFromPDBDirect(Model *pModel, Data *pData, std::string name, std::string PDBFilename, unsigned int numberOfUniqueConnectedPhases){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    totalNumberOfPhasesForSeededModeling = numberOfUniqueConnectedPhases;
    contactCutOff = interconnectivityCutOff;
    auto lowTempStop = (double)highTempStartForCooling;

    // convert distances within the large search space to bins based on input P(R)-DATA file

    maxbin = pModel->getMaxBin(pData);
    maxbin += 1;  // maximum bin index, so vector size must be +1
    std::cout << "         INITIAL MAX BIN : " << maxbin << std::endl;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    //this->fillPrBinsAndAssignTotalBin( pModel, pData);
    //unsigned short int * const pBin = pModel->getPointerToBins(); // initialized as emptyin Model class

    std::vector<double> prPDB(maxbin);
    pModel->createSeedFromPDB(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size
    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);
    auto trueModelBeginIt = pModel->getSeedBegin();
    unsigned int workingLimit = pModel->getTotalInSeed();
    // copy trueModel into BeadIndices
    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    /*
     * Create set of Anchors if exists
     * Anchors must be included in the model and connected to the model
     * Can have more than one anchor
     *
     * totalComponents default is 0, set by setAnchorPoints
     */
    std::set<unsigned int> excludeAnchorsList;
    for(unsigned int i=0; i < totalComponents; i++) {
        // the selected set of beads that make up each component will be held by Component object
        // add randomly selected indices to Component
        Component * component = &(components[i]);
        if (!component->anchorsEmpty()){
            for(auto it = component->getAnchors()->begin(); it != component->getAnchors()->end(); ++it) {  // find the anchor
                excludeAnchorsList.insert(*it);
            }
            component->writeAnchorsToFile("anchors_"+std::to_string(i+1));
        }
    }

    std::cout << "    TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "    MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "              BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "           BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    auto minWorkingLimit = (unsigned int)(0.17*workingLimit);

    std::vector<unsigned int> bead_indices(workingLimit); // large vector ~1000's
    std::vector<unsigned int> lowest_bead_indices(workingLimit); // large vector ~1000's
    std::vector<unsigned int> backUpState(workingLimit);

    // c-style is slightly faster for large vector sizes
    const unsigned int num = workingLimit;
    unsigned int * ptr = (num != 0) ? &bead_indices.front() : nullptr;
    for(unsigned int i = 0; i < num; i++) {
        ptr[i] = *(trueModelBeginIt+i);
    }
    std::cout << "        4a: " << std::endl;
    // prepare bead_indices by copying in truemodel
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    // check for anchors, if anchors are present, move them out of the core model (anchors are only considered in the component)
//    std::cout << "                       ANCHORS ?  => " << excludeAnchorsList.size() << endl;
//    std::cout << "              INITIAL WORKINGLIMIT   " << workingLimit << endl;
//    if(excludeAnchorsList.size() > 0){
//        for(int i=0; i<workingLimit; i++){
//            if (excludeAnchorsList.find(bead_indices[i]) != excludeAnchorsList.end()){
//                // move to workingLimit and decrement
//                int decrement=1;
//                while( !( excludeAnchorsList.find(*(bead_indices.begin()+(workingLimit-decrement))) == excludeAnchorsList.end() ) ){
//                    decrement--;
//                }
//                cout << " FOUND ANCHOR AND SWAPPING TO END " << bead_indices[i] << endl;
//                std::iter_swap(bead_indices.begin() + i, bead_indices.begin() + (workingLimit-decrement));
//                workingLimit -= decrement;
//            }
//        }
//    }
//    std::cout << " (AFTER ANCHOR CHECK) WORKINGLIMIT  " << workingLimit << endl;
    const unsigned int deadLimit = workingLimit;; // as bead indices are discarded, set upper limit of vector
    //std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);

    // Surface area calculations
    std::vector<unsigned int> unsorted_bead_indices(workingLimit);   // large vector ~1000's
    std::vector<float> weights(workingLimit);
    std::fill(weights.begin(), weights.end(), pModel->getBeadRadius() + delta_r);
    std::vector<Eigen::Vector3f> coordinates(workingLimit);

    for(unsigned int i=0; i < workingLimit; i++){
        const vector3 & pBeadVec = pModel->getBead(bead_indices[i])->getVec();
        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
        unsorted_bead_indices[i] = bead_indices[i];
    }

    float currentStoV = surfaceToVolume(workingLimit, weights, coordinates);
    float startingStoV = currentStoV;

    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.end());
    double inv_kb_temp = 1.0/(double)highTempStartForCooling;

    auto lowerN = (unsigned int)(0.1*workingLimit);
    unsigned int upperN = workingLimit;
    std::cout << " TOTAL LATTICE IN SEED : " << workingLimit << std::endl;
    std::cout << "        LATTICE LIMITS : " << lowerN << " <= N <= " << upperN << std::endl;
    // randomize and take the workingLength as first set, shuffling takes about 10x longer than copy and sort

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isConnected = true;
    // bead_indices contains only the indices that relate to the input PDB

    std::cout << " CHECKING CONNECTIVITY  " << std::endl;
    std::cout << "        IS CONNECTED ? : " << isConnected << std::endl;
    std::cout << "  NUMBER OF COMPONENTS : " << currentNumberOfComponents << std::endl;

    //const int components = tempNumberOfComponents;
    // calculate Pr distribution 0.000865 so 10000*100 is 13 minutes
    auto beginBinCount = binCount.begin();
    auto endBinCount = binCount.end();

    // calculate KL against known
    // calculateKLEnergy populates binCount based on model in bead_indices
    // create custom D_KL function supplying Pr_of_PDB as target

    calculateModelPrDistributionDirect(&bead_indices, &binCount, workingLimit, pModel, pData);

    float invShannonNumber = 1.0f/(float)pData->getShannonBins();
    float testKL, currentKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB)*invShannonNumber;

    //float tempTotalContactEnergy = calculateTotalContactEnergy(&bead_indices, workingLimit, pModel, pDistance);
    std::copy(beginBinCount, endBinCount, binCountBackUp.begin());
    std::cout << "          INITIAL D_KL : " << currentKL << std::endl;
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    unsigned int lowestWorkingLimit = workingLimit;
    // int priorWorkingLimit;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    //output for plotting
    bool isUpdated = false;

    unsigned int seedHighTempRounds = 31*deadLimit;

    // calculate average number of contacts per bead
    float startContactSum=0.0f;
    for (unsigned int i=0; i<workingLimit; i++){
        startContactSum += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
    }
    startContactSum = startContactSum/(float)workingLimit;

    float this_energy, startKL = currentKL, lowestKL = currentKL;
    float current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);
    char addRemoveText[50];

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    unsigned int counter=1, original;
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased

    unsigned int high = 0;

    for (; high < seedHighTempRounds; high++){ // iterations during the high temp search

        std::cout << "******************************************************************************" << std::endl;

        if (distribution(gen) >= 0.41){ // ADD BEAD?

            if (workingLimit < (deadLimit-3)){
                std::cout << "*******************                  ADD                   *******************" << std::endl;
                std::sprintf(addRemoveText, " ");

                original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                auto itIndex = std::find(bead_indices.begin(), bead_indices.end(), original);
                auto distance = (unsigned int) std::distance(bead_indices.begin(), itIndex);

                while(original == noNeigborIndex || distance >= deadLimit){
                    original = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                    itIndex = std::find(bead_indices.begin(), bead_indices.end(), original);
                    distance = (unsigned int)std::distance(bead_indices.begin(), itIndex);
                }

                addLatticPositionToModelDirect(&bead_indices, &workingLimit, &itIndex);
                addToPrDirect(original, bead_indices, workingLimit, binCount, pModel, pData);

                testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB)*invShannonNumber;

                // I AM ONLY ADD POSITIONS THAT ARE IN CONTACT VIA NEIGHBORS LIST
                //this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents);
                this_energy = testKL;
                beads_in_use_tree.insert(original);

                if ( this_energy < current_energy ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = 1;
                    isUpdated = true;
                    eulerTour.addNode(original, pModel);
                    std::sprintf(addRemoveText, "     ADD => %i", 1);
                } else if ( exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = 1;
                    isUpdated = true;
                    eulerTour.addNode(original, pModel);
                    std::sprintf(addRemoveText, "     ADD => %i", 1);
                } else { // undo changes (rejecting)
                    beads_in_use_tree.erase(original);
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }
            }

        } else { // REMOVE BEADS?
            std::cout << "*******************                 REMOVE                 *******************" << std::endl;
            // test for deletion
            std::sprintf(addRemoveText, "     REMOVE => %i", 1);
            original = bead_indices[randomIndex(gen)];
            //tempNumberOfComponents = eulerTour.removeNode(original);

            bool tourtest = true;
            while(tourtest){
                if (eulerTour.removeNode(original) == 1){
                    tourtest = false;
                } else {
                    eulerTour.addNode(original, pModel);
                    original = bead_indices[ randomIndex(gen) ]; // potential to select the same thing twice
                }
            }

            if (!tourtest){
                // grab from randomized active_indices list
                removeLatticePositionToModelDirect(bead_indices, binCount, &workingLimit, &original, pModel, pData);
                testKL = calculateKLDivergenceAgainstPDBPR(binCount, prPDB)*invShannonNumber;
                //this_energy = testKL + lambda*connectivityPotential(1);
                this_energy = testKL;

                beads_in_use_tree.erase(original);

                if (this_energy < current_energy ) {
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = 1;
                    isUpdated = true;
                } else if ((testKL > 0 ) && exp((current_energy - this_energy)*inv_kb_temp) > distribution(gen) ){
                    currentKL = testKL;
                    current_energy = this_energy;
                    currentNumberOfComponents = 1;
                    isUpdated = true;
                } else { // undo changes and move to next bead (rejecting)
                    auto beginIt = bead_indices.begin();
                    beads_in_use_tree.insert(original);
                    restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.addNode(original, pModel);
                }
            } else {
                eulerTour.addNode(original, pModel);
            }
        }

        std::cout << "*******************                                        *******************" << std::endl;
        printf("       TEMP => %-.3E ( %.3f )\n     INVKBT => %.2f\n", lowTempStop, acceptRate, inv_kb_temp);
        printf("   MAXSTEPS => %i ( %4i ) \n", seedHighTempRounds, high);
        printf("      GRAPH => %i \n", currentNumberOfComponents);
        printf("LIMIT: %5i ( >= MIN: %i )  \n", workingLimit, minWorkingLimit);
        printf("D_KL: %.4E SCORE : %.4E \n", currentKL, current_energy);

        if (isUpdated){
            acceptRate = inv500*(499*acceptRate+1);
            isUpdated = false;

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::copy(beginBinCount, endBinCount, binCountBackUp.begin()); // make backup copy
            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1); // guaranteed unbiased

        } else {
            acceptRate = inv500*(499*acceptRate);
        }

        updateASATemp(high, seedHighTempRounds, acceptRate, lowTempStop, inv_kb_temp);

        if (currentKL < lowestKL){
            lowestKL = currentKL;
        }

        counter++;
    } // end of HIGH TEMP EQUILIBRATION

    // average number of contacts per bead
    float contactSum=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        contactSum += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
    }

    // determine the distribution of contacts
    // distribution matching
    unsigned int totalCounts=0;
    for (unsigned int c=1; c<13; c++){
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalCounts++;
            }
        }
    }

    std::cout << " DISTRIBUTION OF CONTACTS" << std::endl;
    for (unsigned int c=1; c<13; c++){
        unsigned int totalContactsAt = 0;
        for (unsigned int i=0; i<workingLimit; i++){
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::cout << c << " " << totalContactsAt/(double)totalCounts << std::endl;
    }


    float average_number_of_contacts = contactSum/(float)workingLimit;

    std::cout << "KL DIVERGENCE : " << std::endl;
    std::cout << "                   INITIAL D_KL => " << startKL << std::endl;
    std::cout << "                    LOWEST D_KL => " << lowestKL << std::endl;
    std::cout << "AVERAGE CONTACTS : (per lattice point)" << std::endl;
    std::cout << "                        INITIAL => " << startContactSum << std::endl;
    std::cout << "                       bin[one] => " << ((2.0*binCount[0])/(double)workingLimit) << std::endl;
    // this is fixed model for initial high temp search?
    // set this as the seed
    // remap to bead universe for storage and printing
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
    bead_indices.resize(totalBeadsInSphere);
    ptr = &bead_indices.front();
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    std::vector<unsigned int>::iterator it;
    for (unsigned int i=0; i<workingLimit; i++){
        it = std::find(bead_indices.begin(), bead_indices.end(), lowest_bead_indices[i]); // if itTrueIndex == endTrue, it means point is not within the set
        std::iter_swap(bead_indices.begin()+i, it);
    }

    pModel->setStartingSet(bead_indices);
    pModel->setStartingWorkingLimit(workingLimit);
    // set seed model to be invariant during reconstruction

    char flags[25];
    sprintf(flags, "qhull s FA");
    float cvx = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    /*
     * REMOVE ALL ANCHORS FROM REDUCED MODEL
     * Anchors are not part of the core reduced model but will be part of the model that is built ab initio
     */
    std::cout << "                       ANCHORS ?  => " << excludeAnchorsList.size() << std::endl;
    std::cout << "              INITIAL WORKINGLIMIT   " << workingLimit << std::endl;
    if(!excludeAnchorsList.empty()){
        for (auto eit : excludeAnchorsList){
            auto found = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, eit);
            auto dis = (unsigned int)std::distance(bead_indices.begin(), found);
            if (dis < workingLimit){ // move it
                std::cout << *it <<  " FOUND ANCHOR AND SWAPPING TO END " << " " << workingLimit << std::endl;
                std::iter_swap(found, bead_indices.begin() + workingLimit-1);
                workingLimit--;
            }
        }
    }

    currentStoV = surfaceToVolume(workingLimit, weights, coordinates);
    std::cout << "                      WORKINGLIMIT   " << workingLimit << std::endl;
    std::cout << "                      INITIAL StoV   " << startingStoV << std::endl;
    std::cout << "                        FINAL StoV   " << currentStoV << std::endl;
    std::sort(bead_indices.begin(), bead_indices.begin()+workingLimit);
    pModel->setReducedSeed(workingLimit, bead_indices);

    // Surface area calculations
    unsorted_bead_indices.resize(workingLimit);   // large vector ~1000's
    weights.resize(workingLimit);
    std::fill(weights.begin(), weights.end(), pModel->getBeadRadius() + delta_r);
    coordinates.resize(workingLimit);

    for(unsigned int i=0; i < workingLimit; i++){
        const vector3 & pBeadVec = pModel->getBead(bead_indices[i])->getVec();
        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
        unsorted_bead_indices[i] = bead_indices[i];
    }

    currentStoV = surfaceToVolume(workingLimit, weights, coordinates);

    std::string nameOfModel = pModel->writeModelToFileBare(
            currentKL,
            workingLimit,
            bead_indices,
            binCount,
            name,
            this,
            high,
            cvx,
            average_number_of_contacts);


    //create a distribution from the reduced model
    std::vector<unsigned int> contactsDistributionSeed(13);
    std::fill(contactsDistributionSeed.begin(), contactsDistributionSeed.end(), 0);

    for(unsigned int i = 0;i<workingLimit; i++){
        contactsDistributionSeed[numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i])]++;
    }

    unsigned int totalContactsInDistribution=0;
    for(auto in : contactsDistributionSeed){
        totalContactsInDistribution += in;
    }

    contactsDistribution.resize(13); // set contactDistribution using input PDB model instead of default
    for(unsigned int i = 0;i < 13; i++){
        contactsDistribution[i] = (float)contactsDistributionSeed[i]/(float)totalContactsInDistribution ;
    }

    return true;
}


inline void Anneal::addLatticPositionToModelDirect(std::vector<unsigned int> * pIndices,
                                                   unsigned int * pWorkingLimit,
                                                   std::vector<unsigned int>::iterator * pItIndex){

    // make the swap at the border (workingLimit)
    std::iter_swap(pIndices->begin() + *pWorkingLimit, *pItIndex); // this swaps to working position and changes backup
    // increment workingLimit to include new position
    *pWorkingLimit += 1;
}



/**
 *
 * @param addMe
 * @param beadsInUse - unsorted vector of indices in use
 * @param upperLimit
 * @param prBins
 * @param pModel
 * @param pData
 */
inline void Anneal::addToPrDirect(const unsigned int addMe, std::vector<unsigned int> & beadsInUse, const unsigned int upperLimit, std::vector<unsigned int> & prBins, Model * pModel, Data * pData ){

    const vector3 * pVecAddMe = &pModel->getBead(addMe)->getVec();

    const unsigned int * const ptr = beadsInUse.data();
    for(unsigned int i = 0; i < upperLimit; i++) {
        const vector3 * pVec = &pModel->getBead(ptr[i])->getVec();
        prBins[pData->convertToBinUnsignedInt( (*pVecAddMe - *pVec).length() )]++; // some distances will exceed dmax
    }

    prBins[0]--; // remove the double count from addMe counting itself
}


inline void Anneal::removeLatticePositionToModelDirect(
        std::vector<unsigned int> & bead_indices,
        std::vector<unsigned int> & pBinCount,
        unsigned int * pWorkingLimit,
        const unsigned int * pLatticePointToRemove,
        Model * pModel, Data * pData){

    auto pBeginIt = bead_indices.begin();
    auto itIndex = std::find(pBeginIt, pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
    // remove original from P(r)
    // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
    removeFromPrDirect(*pLatticePointToRemove, bead_indices, *pWorkingLimit, pBinCount, pModel, pData);
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
}


/**
 * beadsInUse is Unsorted
 * lattice point to be removed (removeMe) must be in the vector
 * removeMe is the index of the bead from origination
 *
 */
inline void Anneal::removeFromPrDirect(unsigned int removeMe, std::vector<unsigned int> & beadsInUse, const unsigned int upperLimit, std::vector<unsigned int> & prBins, Model * pModel, Data * pData){

    const vector3 * pVecRemoveMe = &pModel->getBead(removeMe)->getVec();
    for(unsigned int i = 0; i<upperLimit; i++){
        const vector3 * pVec = &pModel->getBead(beadsInUse[i])->getVec();
        prBins[ pData->convertToBinUnsignedInt( (*pVecRemoveMe - *pVec).length() ) ]--;
    }

    prBins[0]++; // correct for double count
}


inline void Anneal::restoreRemovingLatticePointFromBackUpDirect(
        unsigned int * pWorkingLimit,
        std::vector<unsigned int> * pBinCountBackUp,
        std::vector<unsigned int> * pBinCount){
    // value we want to recover is at wl = 9
    // wl += 1 => 10
    // 0 1 2 3 9 5 6 7 8 4 10
    // sort to wl
    // 0 1 2 3 4 5 6 7 8 9 10
    //
    *pWorkingLimit += 1;
    // if I have 5000 lattice points copy from backup is O(n) versus sorting to number of lattice points in working limit
    // sorting is n*log(n) with n = 200?  should be much smaller
    std::copy(pBinCountBackUp->begin(), pBinCountBackUp->end(), pBinCount->begin()); //copy to bin count
}



inline void Anneal::restoreAddingFromBackUpDirect(
                                                  unsigned int * pWorkingLimit,
                                                  std::vector<unsigned int> * pBinCountBackUp,
                                                  std::vector<unsigned int> * pBinCount){

    *pWorkingLimit -= 1;
    std::copy(pBinCountBackUp->begin(), pBinCountBackUp->end(), pBinCount->begin()); //copy to bin count
}


inline void Anneal::fillPrBinsAndAssignTotalBinDirect(Model * pModel, Data * pData){

    maxbin = pModel->getMaxBin(pData) + 1;
    totalBins = pData->getShannonBins(); // set by experimental input Data file
    // maxbin and totalBins will be differenct
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
}