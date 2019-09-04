//
// Created by xos81802 on 10/07/2018.
//

#include <random>
#include "../Anneal.h"
#include "../EulerTour/EulerTour.h"

bool Anneal::createInitialModelCVXHull(Model *pModel, Data *pData, std::string name) {
    std::cout << " --------------------------------------------------- " << std::endl;
    srand(time(0));
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    contactCutOff = interconnectivityCutOff;
    // convert distances to ShannonBin membership
    maxbin = pModel->populateBins(pData);

    char score[9];
    (pData->getIsPr()) ? sprintf(score, "    D_KL") : sprintf(score, "CHI_FREE");
    /*
     * pBin holds a vector of small ints that are used to determine which bins to add to
     */
    unsigned short int * pBin = pModel->getPointerToBins(); // initialized as empty in Model class

    maxbin += 1;  // maximum bin index, so vector size must be +1
    std::cout << "         INITIAL MAX BIN : " << maxbin << std::endl;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
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
    //std::clock_t start = std::clock();
    unsigned int * const ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    pModel->writeSubModelToFile(0, totalBeadsInSphere, bead_indices, "universe");

    // cout << "C-STYLE " << (std::clock() - start)/(double) CLOCKS_PER_SEC << endl;
    std::random_device rd;
    std::mt19937 gen(rd());

//    float invBeadVolume = 1.0f/pModel->getBeadVolume();
    lowerV = (unsigned int)(pData->getVolume()  - pData->getVolume() *0.39);
    upperV = (unsigned int)(pData->getVolume()  + pData->getVolume() *0.05);

    float radius_larger = pModel->getBeadRadius()*std::sqrt(7.0f)/2.0f;
    auto lowerN = (unsigned int)std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger));
    auto upperN = (unsigned int)std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger));
    // pick random number between lowerN and upperN
    unsigned int workingLimit = upperN;//number_of_beads_to_use(gen);
    /*
     * create initial model
     * randomize bead indices
     * sort to workingLimit
     */
    std::shuffle(bead_indices.begin(), bead_indices.end(), gen);
    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);

    char flags[] = "qhull FA";
    float test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    std::set<unsigned int> hull;  // sort hull points into bead_indices
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    //populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
    recalculateDeadLimit(workingLimit, bead_indices, hull, pModel, totalBeadsInSphere);
    std::string switched = "HULL";
    // make copy as initial best model
    unsigned int lowestWorkingLimit = workingLimit;
    std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());

    std::cout << "        WORKINGLIMIT SET : " << workingLimit << std::endl;
    // setup parameters for hull
    coordT points[3*(upperN+(int)(0.30*upperN))];
    //char flags[] = "qhull FA";
    //float test_volume, current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    float lowest_volume = current_volume;
    float initialVolume = current_volume;

    std::cout << "              CVX VOLUME : " << current_volume  << std::endl;

    EulerTour eulerTour(bead_indices, workingLimit, pModel);
    unsigned int tempNumberOfComponents, currentNumberOfComponents  = eulerTour.getNumberOfComponents();
    std::cout << "      CREATED EULER TOUR : " << currentNumberOfComponents  << std::endl;

    calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);

    float currentKL = pData->getScore(binCount);
    float lowest_kl = currentKL;
    float hlambda = pData->getIsIntensity() ? 10 : 0.0001f;
//    float muConstant = pData->getIsPr() ? 0.000001f : 0.1f; // chi tends to start in 100's whereas DKL is 0.1;
    double muPercent = 0.31;
    double muConstant = currentKL*muPercent/((1.0 - muPercent)*current_volume/(double)workingLimit);

    float targetVolume = upperV;
    float diff = std::abs(current_volume-targetVolume)/targetVolume;
    //float temp_volume_energy, current_volume_energy = muConstant*diff;
    float temp_volume_energy, current_volume_energy = muConstant*current_volume/(double)workingLimit;

    float testKL, test_energy, current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + current_volume_energy;

    std::cout << "           INITIAL SCORE : " << currentKL << std::endl;
    float lowest_energy = current_energy;

    auto lowTempStop =  (double)0.001;//highT;
    double inv_kb_temp = 1.0f/lowTempStop;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;
    unsigned int failures=0;

    double sumIt = 0.0, counter=0.0;
    bool modulo = false, updateCVX = false;
    //int sizeOfNeighborhood = pModel->getSizeOfNeighborhood();  seems on average there are less than 12 available to look at
    std::cout << " STARTING CONSTANT TEMP SEARCH " << currentNumberOfComponents << std::endl;
    unsigned int high, hullsize=hull.size();
    unsigned int testIndex, rI, tempWorkingLimit;
    // std::clock_t startTime;
    //exit(0);
    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomHull(0,hullsize); // guaranteed unbiased

    double probabilityAddRemove = 0.361;

    for (high=0; high < highTempRounds; high++) { // iterations during the high temp search

        if (distribution(gen) < probabilityAddRemove ){
            /*
             * reshuffle within hull
             */
            if (distribution(gen) < -0.031 && currentNumberOfComponents > 1){ //scrambling mechanism

//                recalculateDeadLimit(workingLimit, bead_indices, hull, pModel, totalBeadsInSphere);
//                hullsize = hull.size();
//                randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);

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
                tempWorkingLimit = workingLimit; //number_of_beads_to_use(gen);

                std::shuffle(bead_indices.begin(), bead_indices.begin() + limit, gen);
                std::sort(bead_indices.begin(), bead_indices.begin() + tempWorkingLimit);

                calculateModelPrDistribution(&bead_indices, &binCount, tempWorkingLimit, totalBeadsInSphere, pModel, pData);
                testKL = pData->getScore(binCount);

                test_volume = calculateCVXHULLVolume(flags, &bead_indices, tempWorkingLimit, pModel);

                //temp_volume_energy = muConstant*std::abs(test_volume-targetVolume)/targetVolume;
                temp_volume_energy = muConstant*test_volume/(double)tempWorkingLimit;

                EulerTour tempTour(bead_indices.cbegin(), tempWorkingLimit, pModel);
                tempNumberOfComponents=tempTour.getNumberOfComponents();
                test_energy = (float)(testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);

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

                    auto set_it=hull.begin();
                    std::advance(set_it, randomHull(gen));
                    testIndex = *set_it;

                    while ( beads_in_use_tree.find(testIndex) != beads_in_use_tree.end()) { // find a new position
                        //testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                        set_it=hull.begin();
                        std::advance(set_it, randomHull(gen));
                        testIndex = *set_it;
                    }

                    auto itIndex = std::find(bead_indices.begin()+workingLimit, bead_indices.end(), testIndex);
                    //make the swap at the border (workingLimit)

                    addLatticPositionToModel(&bead_indices, &workingLimit, &itIndex);
                    addToPr(testIndex, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);

                    testKL = pData->getScore(binCount);
                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                    // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                    //startTime = std::clock();
                    tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);
                    //std::cout << " runtime " << (std::clock() - startTime)/(double) CLOCKS_PER_SEC << std::endl;
//                    temp_volume_energy = muConstant*std::abs(test_volume-targetVolume)/targetVolume;
                    temp_volume_energy = muConstant*test_volume/(double)workingLimit;

                    test_energy = (testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);
                    if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) {
                        beads_in_use_tree.insert(testIndex);
                        updateCVX = true;
                    } else { // undo changes (rejecting)
                        auto beginBinCount = binCount.begin();
                        restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                        eulerTour.removeNode(testIndex);
                    }

                } else {

                    std::cout << "*******************                 REMOVE                 ******************* "  << std::endl;
                    // grab from randomized active_indices list
                    rI = randomIndex(gen);
                    testIndex = bead_indices[rI];

                    removeLatticePositionToModelByIndex(bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, rI);
//                    removeLatticePositionToModel(bead_indices, binCount, pBin, &workingLimit, totalBeadsInSphere, &testIndex);
                    testKL = pData->getScore(binCount);

                    test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
                    //test_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);
                    //startTime = std::clock();
                    tempNumberOfComponents = eulerTour.removeNode(testIndex);
                    //std::cout << " runtime " << (std::clock() - startTime)/(double) CLOCKS_PER_SEC << std::endl;
//                    temp_volume_energy = muConstant*std::abs(test_volume-targetVolume)/targetVolume;
                    temp_volume_energy = muConstant*test_volume/(double)workingLimit;

                    test_energy = (testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);
                    if (test_energy < current_energy || (exp((current_energy - test_energy)*inv_kb_temp) > distribution(gen))) {
                        beads_in_use_tree.erase(testIndex);
                        updateCVX = true;
                    } else { // undo changes and move to next bead (rejecting)
                        eulerTour.addNode(testIndex, pModel);
                        auto beginIt = bead_indices.begin();
                        auto beginBinCount = binCount.begin();
                        restoreRemovingLatticePointFromBackUp(&beginIt, &workingLimit, &binCountBackUp, &beginBinCount);
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
                std::vector<unsigned int> active_indices(workingLimit);
                unsigned int * const pSelections = bead_indices.data(); // initialized as empty in Model class
                for (unsigned int i = 0; i < workingLimit; i++) {
                    beadToPoint(&points[i*3], pModel->getBead(pSelections[i]));
                    active_indices[i] = pSelections[i];
                }

                // needs to be optimized
                qh_new_qhull(3, workingLimit, points, 0, flags, nullptr, nullptr);
                vertexT * vertices = qh vertex_list;
                auto totalV = (unsigned int)qh num_vertices;

                // only move CVX hull points
                std::vector<unsigned int> indices_to_check(totalV); // large vector ~1000's
                for (unsigned int v = 0; v < totalV; v++) { //
                    indices_to_check[v] = active_indices[qh_pointid( vertices->point)];
                    vertices = vertices->next;
                }

                qh_freeqhull(true);
                std::uniform_int_distribution<unsigned int> randomVertices(0, totalV-1);
                swap1 = indices_to_check[randomVertices(gen)];
                modulo=false;
            }

            // find bead to swap in active set
            auto itIndex = std::find(bead_indices.begin(), bead_indices.begin() + workingLimit, swap1);
            // remove selected index from P(r)
            removeFromPr(swap1, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // find better position
            beads_in_use_tree.erase(swap1);
            eulerTour.removeNode(swap1);
            // find a bead in use, and see if it has an empty neighbor to use
            //unsigned int neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]); // biased towards denser regions
            //while (neighbor == noNeigborIndex || neighbor == swap1){
            //    neighbor = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
            //}
            // if neighbor is not in use, then find it in bead_indices and assign to pSwap2
            // make the swap, sort and update P(r)

            auto set_it=hull.begin();
            std::advance(set_it, randomHull(gen));
            unsigned int neighbor = *set_it;

            while (beads_in_use_tree.find(neighbor) != beads_in_use_tree.end() || neighbor == swap1) { // find a new position
                //testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, bead_indices[randomIndex(gen)]);
                set_it=hull.begin();
                std::advance(set_it, randomHull(gen));
                neighbor = *set_it;
            }

            auto pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), neighbor); // potential problem

            std::iter_swap(itIndex, pSwap2);
            std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted

            addToPr(neighbor, bead_indices, workingLimit, pBin, totalBeadsInSphere, binCount);
            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->getScore(binCount);
            tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

            test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

//            temp_volume_energy = muConstant*std::abs(test_volume-targetVolume)/targetVolume;
            temp_volume_energy = muConstant*test_volume/(double)workingLimit;

            test_energy = (testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + temp_volume_energy);

            if (test_energy < current_energy || (exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))) { //
                beads_in_use_tree.insert(neighbor);
                updateCVX = true;
            } else {
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                eulerTour.removeNode(neighbor);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
                beads_in_use_tree.insert(swap1);
            }
        }

        // Adaptive simulated annealing part
        if (updateCVX){
            currentNumberOfComponents = tempNumberOfComponents;
            currentKL = testKL;
            current_volume = test_volume;
            current_volume_energy = temp_volume_energy;
            current_energy = test_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            updateCVX = false;
            failures=0;

            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());   // make backup copy
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy


            if (currentNumberOfComponents < 13){
                populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
                switched = "NEAREST NEIGHBOR";
            } else {
                switched = "HULL";
                recalculateDeadLimit(workingLimit, bead_indices, hull, pModel, totalBeadsInSphere);
            }

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
        printf("       ACCEPTRATE : %5.3f       STEP : %-7i (%7i)\n", acceptRate, high, highTempRounds);
        printf("            GRAPH : %5i   FAILURES : %i %s\n", currentNumberOfComponents, failures, switched.c_str());
        printf("            LIMIT : %5i   UPPER >= %i LOWER <= %i POOL : %i\n", workingLimit, upperN, lowerN, hullsize);
        printf("           VOLUME : %.0f (%.0f) MU : %.2E \n", current_volume, targetVolume, muConstant);
        printf("         %s : %.4E E_VOL : %.2E ENRGY : %.4E \n", score, currentKL, current_volume_energy, current_energy);
        std::cout << "*******************                                        *******************" << std::endl;


        modulo = (high%171) ? true : modulo;
        // check Pr values
        // uncomment to check update of Pr and direct calc are equal
//        float testKL1 = calculateKLEnergy(&bead_indices, &testBinCount, workingLimit, totalBeadsInSphere, pModel, pData);
//        if (currentKL != testKL1 || checkForRepeats(bead_indices)){
//            cout << " STOPPED POSITIONAL " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << endl;
//            return "stopped";
//        }
//
        // test bead_indices and vector are in sync
//        for(int i=0; i<workingLimit; i++){
//            if (beads_in_use_tree.find(bead_indices[i]) == beads_in_use_tree.end()){
//                cout << "!!!!! [ " << i << " ] NOT FOUND IN TREE => " << bead_indices[i] << endl;
//                exit(0);
//            }
//        }
//
//        for(auto bit : beads_in_use_tree){
//            if(std::find(bead_indices.begin(), bead_indices.begin()+workingLimit, bit) == bead_indices.begin()+workingLimit){
//                std::cout << "!!!!! " << " NOT FOUND in BEAD_INDICES " << bit << std::endl;
//                exit(0);
//            }
//        }
//
//        if(beads_in_use_tree.size() != workingLimit){
//            cout << "!!!!! " << " Incorrect size WORKINGLIMIT "  << endl;
//            exit(0);
//        }

        if (currentNumberOfComponents == 1 && current_energy < lowest_energy){
            lowestWorkingLimit = workingLimit;
            lowest_volume = current_volume;
            sumIt += workingLimit;
            counter += 1.0;
            std::copy(bead_indices.begin(), bead_indices.end(), lowest_bead_indices.begin());
            lowest_energy = current_energy;
            lowest_kl = currentKL;
            //pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_HT_" +  std::to_string(high));
        }

    }

    highTempStartForCooling = (float)lowTempStop;

    EulerTour tempEulerTour(lowest_bead_indices.begin(), lowestWorkingLimit, pModel);
    std::cout << " Lowest Tour : " << tempEulerTour.getNumberOfComponents() << std::endl;

    std::set<unsigned int> centeringhullpts;
    std::set<unsigned int> lowestSet(lowest_bead_indices.begin(), lowest_bead_indices.begin()+lowestWorkingLimit);
    getHullPoints(centeringhullpts, lowestSet, pModel);

    pModel->centerLatticeModel(&lowestWorkingLimit, lowest_bead_indices, centeringhullpts);

    pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, name, high);

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
    std::cout << "*******************                    KL                  *******************" << std::endl;
    printf("           LOWEST => %.4E \n", lowest_kl);

    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "    SHANNON BINS IN DATA : " << pData->getShannonBins() << std::endl;

    if (pData->getIsPr()){
        std::cout << "            ZERO BIN AT : " << pData->getZeroBin() << std::endl;
    }

//    pData->printICalc(binCount);
    std::cout << " EXITING INITIAL MODEL BUILD => " << currentNumberOfComponents << std::endl;
    float tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }
//
    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    std::cout << " AVERAGE NUMBER CONTACTS : " << average_number_of_contacts << std::endl;

    if (tempEulerTour.getNumberOfComponents() == 1 ) {
        return true;
    } else {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, -g (current value => " << highTempRounds << " )" << std::endl;
        return true;
    }
}

bool Anneal::initializeModelToRefine(Model *pModel, Data *pData, std::string name, std::string PDBFilename) {

    srand(time(0));
    contactCutOff = interconnectivityCutOff;

//    unsigned long int totalDistancesInSphere = pModel->getTotalDistances();
//    float * pDistance = pModel->getPointerToDistance();

    // convert distances within the large search space to bins based on input P(R)-DATA file
    this->fillPrBinsAndAssignTotalBin( pModel, pData);

    // initialize Universe and fill indices
    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's

    unsigned int * ptr = (totalBeadsInSphere != 0) ? &bead_indices.front() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    // returns only seed, needs to be mapped to entire Universe
    std::vector<double> prPDB(maxbin);
    pModel->createSeedFromPDB(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size
    unsigned int workingLimit = pModel->getTotalInSeed();


    // sort reassigned bead model into working universe
    unsigned int locale=0;
    for(auto it = pModel->getSeedBegin(); it != pModel->getSeedEnd(); ++it) {
        auto fit = std::find(bead_indices.begin(), bead_indices.end(), *it); // if itTrueIndex == endTrue, it means point is not within the set
        std::iter_swap(bead_indices.begin()+locale, fit);
        locale++;
    }

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);
    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50

    std::cout << "      TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "                BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "             BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    calculateModelPrDistribution(&bead_indices, &binCount, workingLimit, totalBeadsInSphere, pModel, pData);
    float currentKL = pData->getScore(binCount);
    //
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::cout << " WORKINGLIMIT SET => " << workingLimit << std::endl;
    // setup parameters for hull
    char flags[] = "qhull FA";

    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    // populate deadLimit
    // layer of beads within interconnectivity cutOff
    // randomize points within deadlimit
    // changing working limit to add/remove
    // for each randomized arrangement, does volume get smaller?, connectivity improve?  D_KL?
    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents  = eulerTour.getNumberOfComponents();

//    float lowest_volume, test_volume, current_volume = pData->calculateVolumeFromDistribution(maxbin, &binCount);//calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//    float initialVolume = current_volume;

    std::cout << "   STARTING SCORE => " << currentKL << std::endl;
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    std::cout << " STARTING CONSTANT TEMP SEARCH " << currentNumberOfComponents<< std::endl;

    pModel->writeModelToFile(workingLimit, bead_indices, name, 0);
    pModel->setStartingSet(bead_indices);
    pModel->setStartingWorkingLimit(workingLimit);
    pModel->setBeadAverageAndStdev(workingLimit, workingLimit*0.1f);

    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    std::printf("          AVERAGE => %.0f (%.0f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************                 VOLUME                 *******************" << std::endl;
    std::printf("            FINAL => %.0f\n", current_volume);
    std::cout << "*******************                                        *******************" << std::endl;

    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "    SHANNON BINS IN DATA : " << pData->getShannonBins() << std::endl;
    std::cout << "             ZERO BIN AT : " << pData->getZeroBin() << std::endl;

    std::cout << " EXITING INITIAL MODEL BUILD => " << currentNumberOfComponents << std::endl;
    unsigned int tempAverageContacts=0;
    for (unsigned int i=0; i<workingLimit; i++){
        tempAverageContacts += numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
    }

    float average_number_of_contacts = (float)tempAverageContacts/(float)workingLimit;
    std::cout << " AVERAGE NUMBER CONTACTS : " << average_number_of_contacts << std::endl;

    return true;
}

