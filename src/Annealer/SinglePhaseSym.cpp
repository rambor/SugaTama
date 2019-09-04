//
// Created by xos81802 on 13/07/2018.
//
#include "../Anneal.h"
#include "../EulerTour/EulerTour.h"

/**
 * create initial model of complete object
 * does not require pre-computed distances, this is a direct method meaning all pairwise distances are calcualted on the fly
 */
bool Anneal::createInitialModelSymmetry(Model *pModel, Data *pData) {

    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    auto lowTempStop =  (double)0.0001;
    double inv_kb_temp = 1.0f/lowTempStop;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    contactCutOff = pModel->getNieghborCutOffLimit();
    violation_limit = pModel->getBeadRadius();
    //violation_limit = sqrt(3)*pModel->getBeadRadius();
    //violation_limit = (float)(2.0/3.0*(sqrt(3)*pModel->getBeadRadius()));

    maxbin = pModel->getMaxBin(pData);
    // adjust this part
    maxbin += 2;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> lowestbinCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);

    std::cout << "      TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "                BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "             BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // as bead indices are discarded, set upper limit of vector
    std::vector<unsigned int> subUnit_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<unsigned int> test_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<unsigned int> lowest_subUnit_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> active_indices(totalBeadsInSphere);         // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    const unsigned int num = totalBeadsInSphere;
    unsigned int * const ptr = (num != 0) ? subUnit_indices.data() : nullptr;
    for(unsigned int i = 0; i < num; i++) {
        ptr[i] = i;
    }

//    printSymModel(&subUnit_indices, totalBeadsInSphere, pModel, pData);
//    std::set<unsigned int> universe(subUnit_indices.begin(), subUnit_indices.end());
//    std::string universename = "uni";
//    pModel->writeSetToFile(universe, universename);
//    universe.clear();

    lowerV = (unsigned int)(pData->getVolume()  - pData->getVolume() *0.331);
    upperV = (unsigned int)(pData->getVolume()  + pData->getVolume() *0.05);

    const double targetSubUnitVolume = (2*lowerV+upperV)/3/(float)pModel->getNumberOfSubUnits();
    const auto radius_larger = (float)(pModel->getBeadRadius()*std::sqrt(7.0)/2.0);
    const auto lowerN = (unsigned int)(std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()))/3;
    const auto upperN = (unsigned int)(std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()));

//    const auto lowerN = (unsigned int)(std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()))/3;
//    const auto upperN = (unsigned int)(std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()))/2;

    unsigned int subUnitWorkingLimit = (upperN + lowerN)/2;
    const unsigned int multipleWL = subUnitWorkingLimit*7;

    std::cout << "  Bead Search Limited to : " << lowerN << " <= N <= " << upperN << std::endl;
    std::cout << "        INITIAL MODEL WL : " << subUnitWorkingLimit << std::endl;
    std::cout << "                SYMMETRY : " << pModel->getSymmetryString() << std::endl;
    std::cout << "          TOTAL SUBUNITS : " << pModel->getNumberOfSubUnits() << std::endl;

    // randomize and take the workingLength as first set
    // shuffle beads in asymmetric unit
    std::cout << "        CREATING INITIAL RANDOM MODEL " << std::endl;
    std::shuffle(subUnit_indices.begin(), subUnit_indices.end(), gen);
    std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());

    std::set<unsigned int> beads_in_use_tree(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
    std::set<unsigned int> hull;  // sort hull points into bead_indices
    recalculateDeadLimit(subUnitWorkingLimit, subUnit_indices, hull, pModel, totalBeadsInSphere);
    std::string switched = "HULL";
    // randomly pick from each selected index to create monomer
    // each index maps back to beads in universe and also the symmetry grouping
    // calculate volume subunit
    coordT points[3*(2*upperN)];
    char flags[25];
    std::sprintf(flags, "qhull s FA");

    float subunit_test_volume, subunit_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);
//    float test_volume, current_volume = calculateCVXVolumeSymmetry(&subUnit_indices, subUnitWorkingLimit, pModel);

    // calculate starting energy
    unsigned int totalViolations, temp_violations;
    std::clock_t start = std::clock();
    calculateModelPrDistributionSym(&subUnit_indices, &binCount, subUnitWorkingLimit, totalViolations, pModel, pData );

    // fill binCount for first time
    double currentKL = pData->getScore(binCount);

    double hlambda = pData->getIsIntensity() ? 10.0d : 0.01d;
//    double muConstant = pData->getIsPr() ? 0.00001d : 0.1d; // chi tends to start in 100's whereas DKL is 0.1;
    double muPercent = 0.731;
    double muConstant = currentKL*muPercent/((1.0 - muPercent)*subunit_volume/(double)subUnitWorkingLimit);
//    const float targetVolume = (0.5f*(lowerV+upperV));///(float)pModel->getNumberOfSubUnits();
//    const auto lowerSubUnitV = (unsigned int)((pData->getVolume()  - pData->getVolume() *0.21)/(float)pModel->getNumberOfSubUnits());
//    const auto upperSubUnitV = (unsigned int)((pData->getVolume()  + pData->getVolume() *0.05)/(float)pModel->getNumberOfSubUnits());

    //float diff = biCameralVolumePotential(lowerSubUnitV, upperSubUnitV, targetSubUnitVolume);

    //float invTarVolmuConstant = muConstant/(targetVolume*targetVolume);
    const double invTarVolSubUnitmuConstant = muConstant/targetSubUnitVolume;
    unsigned int failures=0;

    bool updateCVX = false;
    std::cout << " THREAD TIME "  << (std::clock() - start)/(double) CLOCKS_PER_SEC << " seconds " << std::endl;
    // high temp search
    // pModel->writeModelToFile(groupWorkingLimit, subUnitIndices, "symStart");

    double testKL, test_energy; // sets alpha as a constant during high temp eq

    unsigned int lowestWorkingLimit = subUnitWorkingLimit;
    unsigned int tempNumberOfComponents, currentNumberOfComponents;

    EulerTour eulerTour(subUnit_indices.begin(), subUnitWorkingLimit, pModel);
    currentNumberOfComponents = eulerTour.getNumberOfComponents();

    double current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) ;
    //double current_energy = currentKL + hlambda*connectivityPotential(currentNumberOfComponents);
    std::cout << "       INITIAL => ENERGY :  " << current_energy << std::endl;
    std::cout << "       INITIAL =>   D_KL : " <<  currentKL << std::endl;
    std::cout << "       INITIAL => VOLUME : " << subunit_test_volume << std::endl;
    std::cout << "              VIOLATIONS : " << totalViolations << std::endl;
    std::cout << "       INITIAL        WL : " << subUnitWorkingLimit << std::endl;

    char addRemoveText[50];
    float sum_x_squared=0;
    unsigned int volumeCount=0, testIndex;
    float volumeSum=0, workingLimitSum=0;

    // what is tolerable number of violations per subsunit?
    current_energy += beta * totalViolations;

    double lowest_energy = current_energy;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
    std::copy(binCount.begin(), binCount.end(), lowestbinCount.begin());  // unaltered P(r)
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup cop

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    unsigned int high, hullsize=hull.size();
    std::uniform_int_distribution<unsigned int> randomHull(0,hullsize-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomIndex(0,subUnitWorkingLimit-1); // guaranteed unbiased

    //double temp_subunit_volume_energy, current_subunit_volume_energy = std::abs(subunit_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
    double temp_subunit_volume_energy, current_subunit_volume_energy = muConstant*subunit_volume/(double)subUnitWorkingLimit;

    unsigned int addFailed = 0, removeFailed =0;
    for (high=0; high < highTempRounds; high++){ // iterations during the high temp search

        if (distribution(gen) < 0.17 || currentNumberOfComponents > 5){
            // create points from workSet to determine HULL
            std::sprintf(addRemoveText, "POSITIONAL");
            unsigned int position = randomIndex(gen);
            unsigned int swap1 = subUnit_indices[position];
            auto itIndex = subUnit_indices.begin() + position;

            // setup parameters for hull
            if ( (distribution(gen) < 0.37129) || currentNumberOfComponents > 2 ){ // only move CVX points if true

                unsigned int * const pSelections = subUnit_indices.data(); // initialized as empty in Model class
                for (unsigned int i = 0; i < subUnitWorkingLimit; i++) {
                    beadToPoint(&points[i*3], pModel->getBead(pSelections[i]));
                    active_indices[i] = pSelections[i];
                }

                // needs to be optimized
                qh_new_qhull(3, subUnitWorkingLimit, points, 0, flags, nullptr, nullptr);
                vertexT * vertices = qh vertex_list;
                auto totalV = (unsigned int)(qh num_vertices);

                // only move CVX hull points
                std::vector<unsigned int> indices_to_check(totalV); // large vector ~1000's
                for (unsigned int v=0; v < totalV; v++) { //
                    indices_to_check[v] = active_indices[qh_pointid(vertices->point)];
                    vertices = vertices->next;
                }

                qh_freeqhull(true);

                std::uniform_int_distribution<unsigned int> randomVertices(0, totalV-1);
                swap1 = indices_to_check[randomVertices(gen)];
                itIndex = std::find(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit, swap1);
            } // recalculate

            // find bead to swap in active set
            removeFromPrSym(swap1, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData); // remove selected index from P(r)
            //beads_in_use_tree.erase(swap1);
            eulerTour.removeNode(swap1);
            // find position to swap to
            auto set_it=hull.begin();
            std::advance(set_it, randomHull(gen));
            unsigned int neighbor = *set_it;
            while (neighbor == swap1 || beads_in_use_tree.find(neighbor) != beads_in_use_tree.end()) { // find a new position
                set_it=hull.begin();
                std::advance(set_it, randomHull(gen));
                neighbor = *set_it;
            }

            // make the swap, sort and update P(r)
            auto pSwap2 = std::find(subUnit_indices.begin()+subUnitWorkingLimit, subUnit_indices.end(), neighbor);

            std::iter_swap(itIndex, pSwap2);

            temp_violations = addToPrSym(neighbor, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);
            testKL = pData->getScore(binCount); // 100x faster than calculateKLEnergySymmetry
            tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

            subunit_test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel); // calculate volume of the subUnit
            //temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
            temp_subunit_volume_energy = muConstant*subunit_test_volume/(double)subUnitWorkingLimit;

            //test_energy = testKL + hlambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations;
            test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) +  beta*temp_violations;//diff*invTarVolmuConstant;

            if (( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy-temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                beads_in_use_tree.erase(swap1);
                beads_in_use_tree.insert(neighbor);
                updateCVX = true;
            } else {
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                std::copy(backUpState.begin(), backUpState.end(), subUnit_indices.begin());
                //beads_in_use_tree.insert(swap1);
                eulerTour.removeNode(neighbor);
                eulerTour.addNode(swap1, pModel);
            }

        } else { // add or remove

            // hard limit on number of beads
            if ( (distribution(gen) < 0.51 && subUnitWorkingLimit > lowerN) ||  subUnitWorkingLimit > upperN  ) { // REMOVE beads from sorted list into useable range < deadLimit
                // randomly swap positions with end of workingLength, could remove CVX Hull Point
                std::sprintf(addRemoveText, "  REMOVE  "); // must only consider removing beads that do not increase number of components

                testIndex = subUnit_indices[randomIndex(gen)];
                while(true){
                    tempNumberOfComponents = eulerTour.removeNode(testIndex);
                    if (tempNumberOfComponents <= currentNumberOfComponents){
                        break;
                    }
                    eulerTour.addNode(testIndex, pModel);
                    testIndex = subUnit_indices[randomIndex(gen)];
                }


                removeLatticePositionToModelSym(subUnit_indices, binCount, &subUnitWorkingLimit, &testIndex, pModel, pData);
                temp_violations = getViolations(&subUnit_indices, subUnitWorkingLimit, pModel);

                testKL = pData->getScore(binCount);

                subunit_test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel); // calculate volume of the subUnit
                //temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
                temp_subunit_volume_energy = muConstant*subunit_test_volume/(double)subUnitWorkingLimit;
//                tempNumberOfComponents = eulerTour.removeNode(testIndex);

                    test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations;
//                test_energy = testKL + hlambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations;

                if ( ( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy - temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                    beads_in_use_tree.erase(testIndex);
                    updateCVX=true;
                } else { // undo changes and move to next bead (rejecting)
                    eulerTour.addNode(testIndex, pModel);
                    std::copy(backUpState.begin(), backUpState.end(), subUnit_indices.begin());
                    subUnitWorkingLimit += 1;
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); //copy to bin count
                    removeFailed += 1;
                }

            } else { // ADD beads
                std::sprintf(addRemoveText, "    ADD   ");

                std::set<unsigned int>::const_iterator set_it=hull.begin();
                std::advance(set_it, randomHull(gen));
                testIndex = *set_it;

                while ( beads_in_use_tree.find(testIndex) != beads_in_use_tree.end()) { // find a new position
                    set_it=hull.begin();
                    std::advance(set_it, randomHull(gen));
                    testIndex = *set_it;
                }

                auto itIndex = std::find(subUnit_indices.begin()+ subUnitWorkingLimit, subUnit_indices.end(), testIndex);
                std::iter_swap(subUnit_indices.begin() + subUnitWorkingLimit, itIndex); // this swaps to working position and changes backup
                // increment workingLimit to include new position

                if (itIndex == subUnit_indices.end()){
                    std::cout << "Not found";
                    exit(0);
                }

                subUnitWorkingLimit += 1;

                temp_violations = addToPrSym(testIndex, subUnit_indices, subUnitWorkingLimit, binCount, pModel, pData);

                testKL = pData->getScore(binCount);

                subunit_test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel); // calculate volume of the subUnit
                //temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
                temp_subunit_volume_energy = muConstant*subunit_test_volume/(double)subUnitWorkingLimit;
                // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                test_energy = testKL + hlambda*(tempNumberOfComponents-1)*(tempNumberOfComponents-1) + beta*temp_violations;// + diff*invTarVolmuConstant + beta*temp_violations;
                //test_energy = testKL + hlambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations;

                if (( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy-temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                    beads_in_use_tree.insert(testIndex);
                    updateCVX=true;
                } else { // undo changes (rejecting)
                    auto beginBinCount = binCount.begin();
                    restoreAddingFromBackUp(&subUnit_indices, &backUpState, &subUnitWorkingLimit, &binCountBackUp, &beginBinCount);
                    eulerTour.removeNode(testIndex);
                    addFailed += 1;
                }
            }
        }



        if (updateCVX){

            current_energy = test_energy;
            currentKL = testKL;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            subunit_volume = subunit_test_volume;
            current_subunit_volume_energy = temp_subunit_volume_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            updateCVX = false;
            failures=0;

            if ( hullsize < subUnitWorkingLimit/2 || currentNumberOfComponents == 1){
                populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
                switched = "NEAREST NEIGHBOR";
            } else {
                switched = "HULL";
                recalculateDeadLimit(subUnitWorkingLimit, subUnit_indices, hull, pModel, totalBeadsInSphere);
            }

            hullsize = hull.size();
            randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);

            std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy
            randomIndex = std::uniform_int_distribution<unsigned int>(0,subUnitWorkingLimit-1); // guaranteed unbiased

            if (current_subunit_volume_energy > 1.5*currentKL && currentNumberOfComponents == 1){
                muConstant = currentKL*muPercent/((1.0 - muPercent)*subunit_volume/(double)subUnitWorkingLimit);
                current_subunit_volume_energy = muConstant*subunit_volume/(double)subUnitWorkingLimit;
            }
        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASAConstantTemp(high, highTempRounds, acceptRate, lowTempStop, inv_kb_temp);

        // if no positional refinement using populate, the deadlimit space is fine
        // if using positional refinement with populate, get artefacts

        std::printf("*******************             %s                 ******************* \n", addRemoveText);
        std::printf("      MAXSTEPS : %7i (%7i) ACCEPTRATE : %.3f TEMP : %.2E\n", highTempRounds, high, acceptRate, lowTempStop);
        std::printf("         GRAPH : %7i  MODE => %s \n", currentNumberOfComponents, switched.c_str());
        std::printf("         LIMIT : %7i UPPER <= %i  LOWER >= %i (%i)\n", subUnitWorkingLimit, upperN, lowerN, hullsize);
        std::printf("        VOLUME : %.0f (%.0f)  MU => %.1E  MU*VOL => %.2E\n", subunit_volume, targetSubUnitVolume, muConstant, current_subunit_volume_energy);
        std::printf("    VIOLATIONS : %7d  BETA => %.1E (%i) %i | %i\n", totalViolations, beta, volumeCount, addFailed, removeFailed);
        std::printf("          D_KL : %.4E ENRGY: %.4E (%.4E) \n", currentKL, current_energy, lowest_energy);
        std::cout << "*******************                                        *******************" << std::endl;

        // write to file to animate search
        if (currentNumberOfComponents == 1 && current_energy < lowest_energy ){
            std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());
            std::copy(binCount.begin(), binCount.end(), lowestbinCount.begin());  // unaltered P(r)
            lowestWorkingLimit = subUnitWorkingLimit;

            lowest_energy = current_energy;
            workingLimitSum += subUnitWorkingLimit;
            sum_x_squared += subUnitWorkingLimit*subUnitWorkingLimit;
            volumeSum += subunit_volume;
            volumeCount++;
        }

        if (currentNumberOfComponents > 61 && high%multipleWL == 0){
            hlambda = currentKL*0.971/((1.0 - 0.971)*(currentNumberOfComponents-1)*(currentNumberOfComponents-1));
            current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) + beta*totalViolations;
        }



//        checkSetAndVector(subUnitWorkingLimit, &subUnit_indices, &beads_in_use_tree);
//
//        if (checkForRepeats(subUnit_indices)){
//            std::cout << " STOPPED POSITIONAL " << " WL: " << subUnitWorkingLimit << " D_KL " << currentKL << " <=> " << std::endl;
//            exit(0);
//        }

//        violations = getViolations(&subUnit_indices, subUnitWorkingLimit, pModel, pData);
//        if (violations != totalViolations){
//            std::cout << " totalViolations " << totalViolations << " actual " << violations << std::endl;
//            testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, violations, pModel, pData );
//            if (currentKL != testKL){
//                printf("*******************             %s                 ******************* \n", addRemoveText);
//                std::cout << high  << "           Failed " << pData->getScore(binCount) << " != " << testKL << std::endl;
//                std::cout << high  << " Failed currentKL " << currentKL << " != " << testKL << std::endl;
//                printf("CKL %.8E  TKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//            }
//            exit(0);
//        }
//
//        calculateModelPrDistributionSym(&subUnit_indices, &testBinCount, subUnitWorkingLimit, temp_violations, pModel, pData );
//        testKL = pData->getScore(testBinCount); // 100x faster than calculateKLEnergySymmetry
//        if (currentKL != testKL){
//            printf("*******************             %s                 ******************* \n", addRemoveText);
//            std::cout << high  << "           Failed " << pData->getScore(binCount) << " != " << testKL << std::endl;
//            std::cout << high  << " Failed currentKL " << currentKL << " != " << testKL << std::endl;
//            printf("CKL %.8E  TKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//            exit(0);
//        }
//
//        testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, temp_violations, pModel, pData );
//        if (currentKL != testKL){
//            // if this fails, update of P(r) is wrong or subUnit_indices is corrupted
//            printf("*******************             %s                 ******************* \n", addRemoveText);
//            std::cout << high  << "           Failed " << pData->getScore(binCount) << " != " << testKL << std::endl;
//            std::cout << high  << " Failed currentKL " << currentKL << " != " << testKL << std::endl;
//            printf("CurrentKL %.8E  TestKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//
//            int bins = testBinCount.size();
//            std::cout << " BIN : " << std::endl;
//            for(int i=0; i<bins; i++){
//                if (testBinCount[i] != binCount[i]){
//                    std::cout << i << " : " << testBinCount[i] << " == " << binCount[i] << std::endl;
//                }
//            }
//
//            exit(0);
//        }

    } // end of HIGH TEMP EQUILIBRATION


    switched = "HULLFinal";
    recalculateDeadLimit(lowestWorkingLimit, lowest_subUnit_indices, hull, pModel, totalBeadsInSphere);
    pModel->writeSetToFile(hull, switched);


    highTempStartForCooling = (float)lowTempStop;
    // pModel->updateBeadIndices(lowestWorkingLimit, lowestDeadLimit, lowest_subUnit_indices);
    // calculate average volume and standard deviation
    float volumeAverage = workingLimitSum/(float)volumeCount;
    float volumeStdev = std::sqrt(sum_x_squared/(float)volumeCount - volumeAverage*volumeAverage);

    // remove points close to hull
    std::string tempName = "initial_CVX_subunit_" + filenameprefix;

    if (currentNumberOfComponents > 1){
        pModel->setStartingSet(lowest_subUnit_indices);
        pModel->setStartingWorkingLimit(lowestWorkingLimit);
        pModel->writeModelToFile(lowestWorkingLimit, lowest_subUnit_indices, tempName , high);
        pModel->writeSymModelToFile(currentKL, lowestWorkingLimit, lowest_subUnit_indices, lowestbinCount, "initial_CVX_sym_" + filenameprefix, this, pData, high, subunit_volume, 4.1);
    } else {
        pModel->setStartingSet(subUnit_indices);
        pModel->setStartingWorkingLimit(subUnitWorkingLimit);
        pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, tempName , high);
        pModel->writeSymModelToFile(currentKL, subUnitWorkingLimit, subUnit_indices, lowestbinCount, "initial_CVX_sym_" + filenameprefix, this, pData, high, subunit_volume, 4.1);
    }

    pModel->setBeadAverageAndStdev(volumeAverage, volumeStdev);

    temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, subUnitWorkingLimit, pModel);
    logger("VIOLATIONS (FINAL)", std::to_string(temp_violations));
    pModel->writeSymModelToFile(currentKL, subUnitWorkingLimit, subUnit_indices, binCount, "final_CVX_sym_" + filenameprefix, this, pData, high, subunit_volume, 4.1);


    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    logger("AVERAGE # BEADS HIGH TEMP SELECTION", std::to_string((int)volumeAverage));
    logger("SIGMA", std::to_string((int)volumeStdev));
    std::cout << "*******************                                        *******************" << std::endl;

    EulerTour finalEulerTour(lowest_subUnit_indices.begin(), lowestWorkingLimit, pModel);
    currentNumberOfComponents = finalEulerTour.getNumberOfComponents();

    if (currentNumberOfComponents == 1 ) {
        return true;
    } else {
        tempName = "failed_subunit_" + filenameprefix;
        pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, tempName , high);
        pModel->writeSymModelToFile(currentKL, subUnitWorkingLimit, subUnit_indices, binCount, "failed_sym" + filenameprefix, this, pData, high, subunit_volume, 4.1);

        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return false;
    }
}


/**
 * Input model must be aligned to a subunit
 *
 * @param pModel
 * @param pData
 * @param nameTo
 * @return
 */
std::string Anneal::refineSymModel(Model *pModel, Data *pData, std::string nameTo){
    const std::string sym = pModel->getSymmetryString();
    bool once = true;
    const char * pnameTo = nameTo.c_str() ;
    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.000000001 : (double)highTempStartForCooling;
    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    std::cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY WITH SYMMETRY" << std::endl;

    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    pModel->copyStartingModelIntoVector(bead_indices);
    unsigned int workingLimit = pModel->getStartingWorkingLimit();

    std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);
    /*
     * create set of beads to select from based on all those that are in contact with selected set
     */
    std::set<unsigned int> hull;
    populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
    std::set<unsigned int> original_hull(beads_in_use_tree.begin(), beads_in_use_tree.end());

//    this->updateContactsDistribution(&beads_in_use_tree, pModel);

    if (maxbin < 2){
        throw std::invalid_argument(" INCORRECT MAXBIN NOT INITIALIZED PROPERLY");
    }

    totalBins = pData->getShannonBins(); //
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin);   // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);      // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);    // smallish vector, typically < 50

    logger("TOTAL EXP N_S BINS", std::to_string(totalBins));
    logger("MAX MODEL N_S BINS", std::to_string(totalBins));
    logger("BINWIDTH (Angstroms)", formatNumber(pData->getBinWidth()));

    unsigned int temp_violations, totalViolations;
    calculateModelPrDistributionSym(&bead_indices, &binCount, workingLimit, totalViolations, pModel, pData );
    // fill binCount for first time
    double testKL, currentKL = pData->getScore(binCount);

    char flags[] = "qhull FA"; // CVX HULL STUFF

    // Surface area calculations
//    std::vector<unsigned int> unsorted_bead_indices(workingLimit);   // large vector ~1000's
//    std::vector<float> weights(workingLimit);
//    probe_radius = pModel->getBeadRadius() + delta_r;
//    std::fill(weights.begin(), weights.end(), probe_radius);
//    std::vector<Eigen::Vector3f> coordinates(workingLimit);
//
//    for(unsigned int i=0; i<workingLimit; i++){
//        vector3 pBeadVec = pModel->getBead(bead_indices[i])->getVec();
//        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
//        unsorted_bead_indices[i] = bead_indices[i];
//    }

    // based on the subunit estimated from initial modeling, use as a upper limit of the volume
    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    float subunit_test_volume;
    const float target_volume = 0.8f*current_volume;
    double muConstant = 0.000001;
    const float invTarVolSubUnitmuConstant = 0;//muConstant/target_volume;
    double temp_subunit_volume_energy=0, current_volume_energy = 0;

    const float targetSubUnitVolume = (upperV)/(float)pModel->getNumberOfSubUnits();
    logger("TARGET VOLUME (A^3)", formatNumber(target_volume, 3));
    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;

    float tempAverageContacts;
    // coupon collector's problem
    const unsigned int couponbase = 3*workingLimit;
    const auto coupons = (unsigned int)(couponbase*std::log((double)couponbase) + 0.5772156649*couponbase + 0.5d);
    const unsigned int updateCount = ccmultiple*coupons;
    float step_limit = (updateCount < 51357) ? 51357 : (float)updateCount;

    std::vector<float> acceptanceRateDuringRun((unsigned long int)step_limit);
    std::vector<double> tempDuringRun((unsigned long int)step_limit);
    std::vector<float> divergenceDuringRun((unsigned long int)step_limit);
    std::vector<unsigned int> workingLimitDuringRun((unsigned long int)step_limit);

    double * pTempDuringRun = &tempDuringRun.front();
    float * pAcceptanceRateDuringRun = &acceptanceRateDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    unsigned int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int tempNumberOfComponents = 0, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool isUpdated = false;
    float acceptRate = 0.5;
    const float inv500slash499 = 499.0f/500.0f, inv500 = 1.0f/500.0f;
    unsigned int original, counter=1, failures=0, updated=0;
    double inv_kb_temp = 1.0/lowTempStop ;

    char addRemoveText[50];

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);

//    double tempContactSum, currentContactsSum = contactPotential(&beads_in_use_tree, pModel);
//    double tempKLDivContacts, currentKLDivContacts = currentContactsSum/(double)beads_in_use_tree.size();
    totalViolations = getViolationsFromSet(&beads_in_use_tree, contactsDistributionOfModel, workingLimit, pModel);
//    populateContactsDistribution(contactsDistributionOfModel, &beads_in_use_tree, pModel);
    double tempKLDivContacts = 0, currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

    auto fifteenPercent = (unsigned int)(0.15*step_limit);
    auto sixtyFivePercent = (unsigned int)(0.65*step_limit);
    auto eightyfivePercent = (unsigned int)(0.85*step_limit);

    double currentCDW = alpha/1000;
    double alphaConstant = alpha > 0 ? std::pow(alpha/currentCDW, 1.0/(double)fifteenPercent) : 0;
    double contactDistributionWeight = (isRefine) ? 0 : currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
    double contactsAvg=3.0d;

    double this_energy, current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts;
    // Add violations
    current_energy += beta*totalViolations;

    double lowest_energy = current_energy ;

    logger("alpha", formatNumber(alpha, 3));
    logger("CDW(alpha)", formatNumber(contactDistributionWeight, 3));
    logger("lambda", formatNumber(lambda, 3));
    logger("beta", formatNumber(beta, 3));

    logger("NUMBER OF COMP", std::to_string(currentNumberOfComponents));
    logger("VIOLATIONS", std::to_string(totalViolations));
    logger("D_KL", formatNumber(currentKL, 7));
    logger("DIV_CON", formatNumber(currentKLDivContacts, 4));
    logger("TOTAL ENERGY", formatNumber(current_energy, 5));

    std::clock_t startTime;
    std::uniform_int_distribution<unsigned int> randomHull(0,hull.size()-1); // guaranteed unbiased

    std::vector<unsigned int> selections(workingLimit);
    for(unsigned int i=0; i < workingLimit; i++){
        selections[i] = i; // some distances will exceed dmax
    }

    const auto totalSubunits = (double)pModel->getNumberOfSubUnits();
    unsigned int numberOfCoolingTempSteps=0;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {

        auto set_it = hull.begin();
        auto beginBinCount = binCount.begin();
        startTime = std::clock();

        if ( distribution(gen) < percentAddRemove ){ //add or remove bead within working Set (exclude deadzone)
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
            // additional points to expand deadlimit will occur via enlarging CVX Hull
            if (distribution(gen) < 0.51){ // ADD BEAD?
                std::cout << "*******************                  ADD                   *******************" << std::endl;

                std::advance(set_it, randomHull(gen));
                original  = *set_it;
                auto itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);

                std::iter_swap(bead_indices.begin() + workingLimit, itIndex); // this swaps to working position and changes backup
                // increment workingLimit to include new position
                workingLimit += 1;

                beads_in_use_tree.insert(original);
                addToPrSym(original, bead_indices, workingLimit, binCount, pModel, pData);
                testKL = pData->getScore(binCount);

                tempNumberOfComponents = eulerTour.addNode(original, pModel);

                temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

//                subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//                temp_subunit_volume_energy = (subunit_test_volume > target_volume) ? (subunit_test_volume-target_volume)*invTarVolSubUnitmuConstant : 0 ;

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

                if ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     ADD => %i", 1);
                } else { // undo changes (rejecting)
                    eulerTour.removeNode(original);
                    beads_in_use_tree.erase(original);
                    std::sprintf(addRemoveText, "     ADD => %i", 0);
                    restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
                }

            } else { // REMOVE BEADS?

                std::cout << "*******************                 REMOVE                 *******************" << std::endl;
                std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
                original = bead_indices[selections[0]];
                // test for deletion

                if (numberOfCoolingTempSteps > eightyfivePercent){
                    for(const auto & select : selections){
                        original = bead_indices[select];
                        if (numberOfContactsFromSet(&beads_in_use_tree, pModel, original) == 1){
                            tempNumberOfComponents = eulerTour.removeNode(original);
                            break;
                        }
                    }
                } else {
                    for(const auto & select : selections){
                        tempNumberOfComponents = eulerTour.removeNode(original);
                        if (tempNumberOfComponents <= currentNumberOfComponents){
                            break;
                        }
                        eulerTour.addNode(original, pModel);
                        original = bead_indices[select];
                    }
                }


                removeLatticePositionToModelSym(bead_indices, binCount, &workingLimit, &original, pModel, pData);
                beads_in_use_tree.erase(original); // remove bead for next calculation

                temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
                testKL = pData->getScore(binCount);

                tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

//                subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//                temp_subunit_volume_energy = (subunit_test_volume > target_volume) ? (subunit_test_volume-target_volume)*invTarVolSubUnitmuConstant : 0 ;

                this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

                if (this_energy < current_energy || (std::exp((current_energy - this_energy)*inv_kb_temp) > distribution(gen) )) {
                    isUpdated = true;
                    std::sprintf(addRemoveText, "     REM => %i", original);
                } else { // undo changes and move to next bead (rejecting)
                    beads_in_use_tree.insert(original);
//                    auto pBeginIt = bead_indices.begin();
//                    restoreRemovingLatticePointFromBackUp(&pBeginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                    workingLimit += 1;
                    // if I have 5000 lattice points copy from backup is O(n) versus sorting to number of lattice points in working limit
                    // sorting is n*log(n) with n = 200?  should be much smaller
                    std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                    std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                    eulerTour.addNode(original, pModel);
                }
            }

        } else { // positional refinement
            std::cout << "______________________________________________________________________________" << std::endl;
            std::cout << "*******************               POSITIONAL               *******************" << std::endl;
            std::cout << "*******************                                        *******************" << std::endl;
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);

            unsigned int position = selections[0];
            unsigned int swap1 = bead_indices[position];

            bool retry = true;
            std::string foundIt = "false";

            if (numberOfCoolingTempSteps > eightyfivePercent){
                for(const auto & select : selections){
                    position = select;
                    swap1 = bead_indices[position];
                    if (numberOfContactsFromSet(&beads_in_use_tree, pModel, swap1) < 2){
                        if ( eulerTour.removeNode(swap1) <= currentNumberOfComponents){
                            retry = false;
                            foundIt = "true";
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
                        foundIt = "true";
                        break;
                    }
                    eulerTour.addNode(swap1, pModel);
                }
            }
            //std::cout << "selected " << foundIt << std::endl;
            // remove contribution of swap1
            auto itIndex = bead_indices.begin() + position;
            // remove selected index from P(r)
            removeFromPrSym(swap1, bead_indices, workingLimit, binCount, pModel, pData);
//            tempContactSum = removeFromContactsPotential(swap1, currentContactsSum, &beads_in_use_tree, pModel);
            beads_in_use_tree.erase(swap1);
            // find better position
            std::advance(set_it, randomHull(gen));
            unsigned int originalSwap2Value = *set_it;
            while (originalSwap2Value == swap1){
                set_it = hull.begin();
                std::advance(set_it, randomHull(gen));
                originalSwap2Value = *set_it;
            }
           // std::cout << "swapped " << *itIndex << " pos " << position << " swap1 " << swap1 <<  " <=> " << originalSwap2Value << std::endl;
            // get available neighbor position
            auto pSwap2 = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), originalSwap2Value);
            std::iter_swap(itIndex, pSwap2);
            //std::sort(bead_indices.begin(), bead_indices.begin() + workingLimit); // bead_indices needs to be sorted
            addToPrSym(originalSwap2Value, bead_indices, workingLimit, binCount, pModel, pData);

            // calculate energy as KL divergence, testKL is ~10x faster than numberOfContacts calculation
            testKL = pData->getScore(binCount);
            beads_in_use_tree.insert(originalSwap2Value);

//            tempContactSum = addToContactsPotential(originalSwap2Value, tempContactSum, &beads_in_use_tree, pModel);
//            tempKLDivContacts=tempContactSum/(double)beads_in_use_tree.size();
            temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
//            populateContactsDistribution(tempContactsDistributionOfModel, &beads_in_use_tree, pModel);
            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);
//            subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//            temp_subunit_volume_energy = (subunit_test_volume > target_volume) ? (subunit_test_volume-target_volume)*invTarVolSubUnitmuConstant : 0 ;
            tempNumberOfComponents = eulerTour.addNode(originalSwap2Value, pModel);
            this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

            if (this_energy < current_energy || (exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen)) ) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     SWAPPED => %i to %i", swap1, originalSwap2Value);
            } else {
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                beads_in_use_tree.insert(swap1); // add back swap one to tree
                beads_in_use_tree.erase(originalSwap2Value);
                std::sprintf(addRemoveText, "      FAILED => %i", swap1);
                eulerTour.removeNode(originalSwap2Value);
                currentNumberOfComponents = eulerTour.addNode(swap1, pModel);
            }
        } // end of positional refinement or add/remove if statement pModel->getVolumeAverage(), pModel->getVolumeStdev()

        printf("      TEMP :  %-.3E MAXSTEPS => %.0f (%4i) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        printf("    ACCEPT : %10.5f FAILURES => %6i  TIME : %.5f\n", acceptRate, failures, ((std::clock() - startTime)/(double) CLOCKS_PER_SEC));
        printf("     LIMIT : %10i     COMP => %3i VIOL => %d EVOL %.2E \n", workingLimit, currentNumberOfComponents, totalViolations, current_volume_energy);
        printf("  CONTACTS :   %.2E (%.2E ALPHA : %.1E) <AVE> : %.2f  \n", (currentKLDivContacts*contactDistributionWeight), currentKLDivContacts, currentCDW, contactsAvg);
        printf("       SYM : %10s  OUTFILE => %s\n", sym.c_str(), pnameTo);
        printf("      D_KL :  %-5.3E ( %.3E ) TOTAL ENRGY : %.4E\n", currentKL, lowest_energy, current_energy);

        // update run time parameters
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pTempDuringRun[numberOfCoolingTempSteps] = totalViolations;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){

            currentKL = testKL;
            current_energy = this_energy;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            currentKLDivContacts = tempKLDivContacts;
            current_volume_energy = temp_subunit_volume_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

//            if (numberOfCoolingTempSteps > sixtyFivePercent){
                populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
//            } else {
//                recalculateDeadLimit(workingLimit, bead_indices, hull, pModel, totalBeadsInSphere);
//            }
            //
            randomHull = std::uniform_int_distribution<unsigned int>(0,(unsigned int)hull.size()-1);
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

            selections.resize(workingLimit);
            unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i; // some distances will exceed dmax
            }

//            currentContactsSum = tempContactSum;
            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), contactsDistributionOfModel.begin());
            contactsAvg = 2.0*binCount[0]/((float)workingLimit*totalSubunits);

        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        //cout << "______________________________________________________________________________" << endl;
        //cout << "*******************                 TEST                   *******************" << endl;
        //cout << "*******************              -----------               *******************" << endl;

//        calculateModelPrDistributionSym(&bead_indices, &testBinCount, workingLimit, temp_violations, pModel, pData );
//        float testKL1 = pData->getScore(binCount);
//        if ((currentKL != testKL1)  || checkForRepeats(bead_indices)){
//            testKL = pData->getScore(binCount);
//            std::cout << "MAIN LOOP " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << " -- " << testKL << std::endl;
//            exit(0);
//        }

        // rescale etaFactor in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent ){
            currentCDW *= alphaConstant;
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
            lowest_energy = current_energy;// + totalContactEnergy ;
        }

        if (once && numberOfCoolingTempSteps > eightyfivePercent && contactDistributionWeight*currentKLDivContacts/currentKL > alpha){
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
            lowest_energy = current_energy;
            once = false;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);

        counter++;
    } // end of steps


    calculateModelParametersSymmetry(&beads_in_use_tree, pModel);

    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (unsigned int i=0; i<workingLimit; i++){
        unsigned int temp = bead_indices[i];
        tempAverageContacts += numberOfContactsFromSet(&beads_in_use_tree, pModel, temp);
    }

    tempAverageContacts = (float)(tempAverageContacts/(double)workingLimit);
    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);

    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);

    // calculate Rg
    std::cout << " CALCULATED RG : " << calculateRgSym(bead_indices, workingLimit, pModel) << std::endl;

    if (currentNumberOfComponents > 1 ) {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return "failed";
    }

    /*
     * calculate distribution of contacts
     */
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
            //contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, contactCutOff, pModel, pDistance);
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/(double)totalCounts);
    }


    // At end of each temp, update a probability model for volume?  Use this to select
    // perform positional refinement until delta E stabilizes?
    // final round to move any points close to body
    float volume_test = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    // pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
    // pModel->setBeadAverageAndStdev(oldN, oldStdev);
    pData->printKLDivergence(binCount);

    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");
    //pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "sphere");

//    string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_final_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);
//    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);

    //pModel->writeModelToFile(workingLimit, bead_indices, "refined_");
    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");

    //  pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_" + filenameprefix);
    //  pModel->writeSymModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_sym_");

    return nameOfModel;
}

/**
 * Input model must be aligned to a subunit
 * Input model is the masked Universe
 * Starting from a filled space, sample down
 * Starting Model has a connected Euler tour
 *
 * @param pModel
 * @param pData
 * @param nameTo
 * @return
 */
std::string Anneal::refineSymModelRefine(Model *pModel, Data *pData, std::string nameTo){
    const std::string sym = pModel->getSymmetryString();
    bool once = true;
    const char * pnameTo = nameTo.c_str() ;
    double lowTempStop =  (asaAcceptanceRate < 0.11) ? (double)highTempStartForCooling*0.0000001 : (double)highTempStartForCooling;
    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    std::cout << "STARTING SA REFINEMENT OF HOMOGENOUS BODY WITH SYMMETRY" << std::endl;

    //unsigned int workingLimit = pModel->getStartingWorkingLimit(); // restrict universe to size of the input model that is going to be refined
    //const unsigned int totalBeadsInSphere = workingLimit;
    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // make copy of bead_indices
    std::vector<unsigned int> bead_indices(totalBeadsInSphere);   // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    /*
     * copy model into bead_indices
     */
    pModel->copyStartingModelIntoVector(bead_indices);
    unsigned int workingLimit = pModel->getStartingWorkingLimit(); // restrict universe to size of the input model that is going to be refined
    printSymModel(&bead_indices, workingLimit, pModel, pData);
    //workingLimit = minimizeViolations(&bead_indices, workingLimit, pModel);
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin()+workingLimit);

    /*
     * create set of beads to select from based on all those that are in contact with selected set
     */
    std::set<unsigned int> hull;
    populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);

    if (maxbin < 2){
        throw std::invalid_argument(" INCORRECT MAXBIN NOT INITIALIZED PROPERLY");
    }

    totalBins = pData->getShannonBins(); //
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);          // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin);   // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);      // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);    // smallish vector, typically < 50

    logger("TOTAL EXP N_S BINS", std::to_string(totalBins));
    logger("MAX MODEL N_S BINS", std::to_string(totalBins));
    logger("BINWIDTH (Angstroms)", formatNumber(pData->getBinWidth()));

    unsigned int temp_violations, totalViolations;
    calculateModelPrDistributionSym(&bead_indices, &binCount, workingLimit, totalViolations, pModel, pData );
    // fill binCount for first time
    double testKL, currentKL = pData->getScore(binCount);

    char flags[] = "qhull FA"; // CVX HULL STUFF

    // Surface area calculations
//    std::vector<unsigned int> unsorted_bead_indices(workingLimit);   // large vector ~1000's
//    std::vector<float> weights(workingLimit);
//    probe_radius = pModel->getBeadRadius() + delta_r;
//    std::fill(weights.begin(), weights.end(), probe_radius);
//    std::vector<Eigen::Vector3f> coordinates(workingLimit);
//
//    for(unsigned int i=0; i<workingLimit; i++){
//        vector3 pBeadVec = pModel->getBead(bead_indices[i])->getVec();
//        coordinates[i] = Eigen::Vector3f(pBeadVec.x, pBeadVec.y, pBeadVec.z);
//        unsorted_bead_indices[i] = bead_indices[i];
//    }

    // based on the subunit estimated from initial modeling, use as a upper limit of the volume
    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    float subunit_test_volume;
    const float target_volume = 0.8f*current_volume;
    double muConstant = 0.000001;
    const float invTarVolSubUnitmuConstant = 0;//muConstant/target_volume;
    double temp_subunit_volume_energy=0, current_volume_energy = 0;

    const float targetSubUnitVolume = (upperV)/(float)pModel->getNumberOfSubUnits();
    logger("TARGET VOLUME (A^3)", formatNumber(target_volume, 3));
    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
    std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

    std::cout << "STARTING ADAPTIVE SIMULATED ANNEALING SEARCH\n STARTING D_KL => " << currentKL << std::endl;

    float tempAverageContacts;
    // coupon collector's problem
    const auto couponbase = workingLimit;
    const auto coupons = (unsigned int)(couponbase*std::log((double)couponbase) + 0.5772156649*couponbase + 0.5d);
    const unsigned int updateCount = ccmultiple*coupons;
    float step_limit = (updateCount < 51357) ? 51357 : (float)updateCount;

    std::vector<float> acceptanceRateDuringRun((unsigned long int)step_limit);
    std::vector<double> tempDuringRun((unsigned long int)step_limit);
    std::vector<float> divergenceDuringRun((unsigned long int)step_limit);
    std::vector<unsigned int> workingLimitDuringRun((unsigned long int)step_limit);

    double * pTempDuringRun = &tempDuringRun.front();
    float * pAcceptanceRateDuringRun = &acceptanceRateDuringRun.front();
    float * pDivergenceDuringRun = &divergenceDuringRun.front();
    unsigned int * pWorkingLimitDuringRun = &workingLimitDuringRun.front();

    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int tempNumberOfComponents = 0, currentNumberOfComponents = eulerTour.getNumberOfComponents();

    bool readyToAdd = true;
    auto readyToAddLimit = (unsigned int)(totalBeadsInSphere*0.71);
    bool isUpdated = false;
    float acceptRate = 0.5;
    const float inv500slash499 = 499.0f/500.0f, inv500 = 1.0f/500.0f;
    unsigned int original, counter=1, failures=0, updated=0;
    double inv_kb_temp = 1.0/lowTempStop ;

    char addRemoveText[50];

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);

//    double tempContactSum, currentContactsSum = contactPotential(&beads_in_use_tree, pModel);
//    double tempKLDivContacts, currentKLDivContacts = currentContactsSum/(double)beads_in_use_tree.size();
    totalViolations = getViolationsFromSet(&beads_in_use_tree, contactsDistributionOfModel, workingLimit, pModel);
    double tempKLDivContacts = 0, currentKLDivContacts = calculateKLDivergenceContactsDistribution(contactsDistributionOfModel);

    auto fifteenPercent = (unsigned int)(0.15*step_limit);
    auto sixtyFivePercent = (unsigned int)(0.65*step_limit);
    auto eightyfivePercent = (unsigned int)(0.85*step_limit);

    double currentCDW = alpha/1000;
    double alphaConstant = alpha > 0 ? std::pow(alpha/currentCDW, 1.0/(double)fifteenPercent) : 0;
    double contactDistributionWeight = (isRefine) ? 0 : currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
    double contactsAvg=3.0d;

    double this_energy, current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + contactDistributionWeight*currentKLDivContacts;
//    double this_energy, current_energy = currentKL + contactDistributionWeight*currentKLDivContacts + beta*totalViolations;
    // Add violations

    double lowest_energy = current_energy ;

    logger("alpha", formatNumber(alpha, 3));
    logger("CDW(alpha)", formatNumber(contactDistributionWeight, 3));
    logger("lambda", formatNumber(lambda, 3));
    logger("beta", formatNumber(beta, 3));

    logger("NUMBER OF COMP", std::to_string(currentNumberOfComponents));
    logger("VIOLATIONS", std::to_string(totalViolations));
    logger("D_KL", formatNumber(currentKL, 7));
    logger("DIV_CON", formatNumber(currentKLDivContacts, 4));
    logger("TOTAL ENERGY", formatNumber(current_energy, 5));

    std::clock_t startTime;
    unsigned int hullsize = hull.size();
    std::uniform_int_distribution<unsigned int> randomHull(0,hull.size()-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomIndex(0,workingLimit-1); // guaranteed unbiased

    std::vector<unsigned int> selections(workingLimit);
    unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
    for(unsigned int i=0; i < workingLimit; i++){
        *(pSelections+i) = i; // some distances will exceed dmax
    }

    const auto totalSubunits = (double)pModel->getNumberOfSubUnits();
    unsigned int numberOfCoolingTempSteps=0;
    bool doIt=false;

    for(; numberOfCoolingTempSteps < step_limit; numberOfCoolingTempSteps++) {

        startTime = std::clock();

        std::cout << "______________________________________________________________________________" << std::endl;
        std::cout << "*******************               ADD?REMOVE               *******************" << std::endl;
        // additional points to expand deadlimit will occur via enlarging CVX Hull
        if (readyToAdd && distribution(gen) < 0.49){ // ADD BEAD?
            std::cout << "*******************                  ADD                   *******************" << std::endl;
            auto set_it = hull.begin();
            std::advance(set_it, randomHull(gen));
            original  = *set_it;
            auto itIndex = std::find(bead_indices.begin() + workingLimit, bead_indices.end(), original);

            std::iter_swap(bead_indices.begin() + workingLimit, itIndex); // this swaps to working position and changes backup
            // increment workingLimit to include new position
            workingLimit += 1;

            beads_in_use_tree.insert(original);
            addToPrSym(original, bead_indices, workingLimit, binCount, pModel, pData);
            testKL = pData->getScore(binCount);

            tempNumberOfComponents = eulerTour.addNode(original, pModel);

            temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

//                subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//                temp_subunit_volume_energy = (subunit_test_volume > target_volume) ? (subunit_test_volume-target_volume)*invTarVolSubUnitmuConstant : 0 ;

            this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;
//                this_energy = testKL + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

            if ( this_energy < current_energy || ( std::exp((current_energy - this_energy) * inv_kb_temp) > distribution(gen) )) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     ADD => %i", 1);
                hull.erase(original);
            } else { // undo changes (rejecting)
                eulerTour.removeNode(original);
                beads_in_use_tree.erase(original);
                std::sprintf(addRemoveText, "     ADD => %i", 0);
                auto beginBinCount = binCount.begin();
                restoreAddingFromBackUp(&bead_indices, &backUpState, &workingLimit, &binCountBackUp, &beginBinCount);
            }

        } else { // REMOVE BEADS?

            std::cout << "*******************                 REMOVE                 *******************" << std::endl;
            std::shuffle(selections.begin(), selections.begin() + workingLimit, gen);
            //original = bead_indices[randomIndex(gen)];
            original = bead_indices[selections[0]];
            // test for deletion

            if (numberOfCoolingTempSteps > sixtyFivePercent){
                for(const auto & select : selections){
                    original = bead_indices[select];
                    if (numberOfContactsFromSet(&beads_in_use_tree, pModel, original) < contactsAvg){

                        tempNumberOfComponents = eulerTour.removeNode(original);
                        if (tempNumberOfComponents <= currentNumberOfComponents){
                            break;
                        }
                        eulerTour.addNode(original, pModel);
                    }
                }
            } else {
                for(const auto & select : selections){
                    tempNumberOfComponents = eulerTour.removeNode(original);
                    if (tempNumberOfComponents <= currentNumberOfComponents){
                        break;
                    }
                    eulerTour.addNode(original, pModel);
                    original = bead_indices[select];
                }
            }

            removeLatticePositionToModelSym(bead_indices, binCount, &workingLimit, &original, pModel, pData);
            beads_in_use_tree.erase(original); // remove bead for next calculation

            temp_violations = getViolationsFromSet(&beads_in_use_tree, tempContactsDistributionOfModel, workingLimit, pModel);
            testKL = pData->getScore(binCount);

            tempKLDivContacts = calculateKLDivergenceContactsDistribution(tempContactsDistributionOfModel);

//                subunit_test_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
//                temp_subunit_volume_energy = (subunit_test_volume > target_volume) ? (subunit_test_volume-target_volume)*invTarVolSubUnitmuConstant : 0 ;

            this_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;
//                this_energy = testKL + beta*temp_violations + contactDistributionWeight*tempKLDivContacts;

            if (this_energy < current_energy || (std::exp((current_energy - this_energy)*inv_kb_temp) > distribution(gen) )) {
                isUpdated = true;
                std::sprintf(addRemoveText, "     REM => %i", original);
                hull.insert(original);
            } else { // undo changes and move to next bead (rejecting)
                beads_in_use_tree.insert(original);
//                    auto pBeginIt = bead_indices.begin();
//                    restoreRemovingLatticePointFromBackUp(&pBeginIt, &workingLimit, &binCountBackUp, &beginBinCount);
                workingLimit += 1;
                // if I have 5000 lattice points copy from backup is O(n) versus sorting to number of lattice points in working limit
                // sorting is n*log(n) with n = 200?  should be much smaller
                std::copy(backUpState.begin(), backUpState.end(), bead_indices.begin());
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); // add back swap 1 from backup
                eulerTour.addNode(original, pModel);
            }
        }

        printf("      TEMP :  %-.3E MAXSTEPS => %.0f (%4i) \n", lowTempStop, step_limit, numberOfCoolingTempSteps);
        printf("    ACCEPT : %10.5f FAILURES => %6i  TIME : %.5f\n", acceptRate, failures, ((std::clock() - startTime)/(double) CLOCKS_PER_SEC));
        printf("       SYM : %10s  OUTFILE => %s\n", sym.c_str(), pnameTo);
        printf("    E_VIOL :   %.2E VIOL => %d LIMIT : %5i COMP %4i HULL %i\n", (beta*totalViolations), totalViolations, workingLimit, currentNumberOfComponents, hullsize);
        printf("  E_CNTCTS :   %.2E (%.2E ALPHA : %.1E) <AVE> : %.2f  \n", (currentKLDivContacts*contactDistributionWeight), currentKLDivContacts, currentCDW, contactsAvg);
        printf("      D_KL :  %-5.3E ( %.3E ) TOTAL ENRGY : %.4E\n", currentKL, lowest_energy, current_energy);

        // update run time parameters
        pAcceptanceRateDuringRun[numberOfCoolingTempSteps] = acceptRate;
        pTempDuringRun[numberOfCoolingTempSteps] = totalViolations;
        pDivergenceDuringRun[numberOfCoolingTempSteps] = currentKL;
        pWorkingLimitDuringRun[numberOfCoolingTempSteps] = workingLimit;

        // Adaptive simulated annealing part
        if (isUpdated){

            currentKL = testKL;
            current_energy = this_energy;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            currentKLDivContacts = tempKLDivContacts;
            current_volume_energy = temp_subunit_volume_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            isUpdated = false;
            failures=0;
            updated++;

            /*
             * do not update hull until target acceptance rate is reached
             */
            if(doIt){
                populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
            }
            hullsize = hull.size();
            randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);
//            randomIndex = std::uniform_int_distribution<unsigned int>(0,workingLimit-1);
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());
            std::copy(bead_indices.begin(), bead_indices.end(), backUpState.begin());

            selections.resize(workingLimit);
            unsigned int * const pSelections = &selections[0]; // initialized as emptyin Model class
            for(unsigned int i=0; i < workingLimit; i++){
                *(pSelections+i) = i; // some distances will exceed dmax
            }

            if (workingLimit <= readyToAddLimit){
                readyToAdd=true;
            }

//            currentContactsSum = tempContactSum;
            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), contactsDistributionOfModel.begin());
            contactsAvg = 2.0*binCount[0]/((float)workingLimit*totalSubunits);

        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        if (!doIt && numberOfCoolingTempSteps > fifteenPercent){
            doIt = true;
        }

        //cout << "______________________________________________________________________________" << endl;
        //cout << "*******************                 TEST                   *******************" << endl;
        //cout << "*******************              -----------               *******************" << endl;
//        calculateModelPrDistributionSym(&bead_indices, &testBinCount, workingLimit, temp_violations, pModel, pData );
//        double testKL1 = pData->getScore(binCount);
//        if ((currentKL != testKL1)  || checkForRepeats(bead_indices)){
//            testKL = pData->getScore(binCount);
//            std::cout << "MAIN LOOP " << " WL: " << workingLimit << " D_KL " << currentKL << " <=> " << testKL1 << " -- " << testKL << std::endl;
//            exit(0);
//        }

        // rescale etaFactor in small steps until final value
        if (numberOfCoolingTempSteps < fifteenPercent ){
            currentCDW *= alphaConstant;
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
//            current_energy = currentKL + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
            lowest_energy = current_energy;// + totalContactEnergy ;
        }

        if (once && numberOfCoolingTempSteps > eightyfivePercent && contactDistributionWeight*currentKLDivContacts/currentKL > alpha){
            contactDistributionWeight = currentKL*currentCDW/((1.0 - currentCDW)*currentKLDivContacts);
            current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
//            current_energy = currentKL + beta*totalViolations + contactDistributionWeight*currentKLDivContacts + current_volume_energy;
            lowest_energy = current_energy;
            once = false;
        }

        updateASATemp(numberOfCoolingTempSteps, step_limit, acceptRate, lowTempStop, inv_kb_temp);

        counter++;
    } // end of steps


    calculateModelParametersSymmetry(&beads_in_use_tree, pModel);

    tempAverageContacts=0.0;
    // calculate average_number_of_contacts
    for (unsigned int i=0; i<workingLimit; i++){
        unsigned int temp = bead_indices[i];
        tempAverageContacts += numberOfContactsFromSet(&beads_in_use_tree, pModel, temp);
    }

    tempAverageContacts = (float)(tempAverageContacts/(double)workingLimit);
    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_annealed_" + nameTo, this, pData, numberOfCoolingTempSteps, current_volume, tempAverageContacts);

    printParameters(&acceptanceRateDuringRun, &tempDuringRun, &divergenceDuringRun, &workingLimitDuringRun);

    // calculate Rg
    std::cout << " CALCULATED RG : " << calculateRgSym(bead_indices, workingLimit, pModel) << std::endl;

    if (currentNumberOfComponents > 1 ) {
        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return "failed";
    }

    /*
     * calculate distribution of contacts
     */
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
            //contactSum += numberOfContacts(lowest_bead_indices[i], &lowest_bead_indices, lowestWorkingLimit, contactCutOff, pModel, pDistance);
            if (numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]) == c){
                totalContactsAt++;
            }
        }
        std::printf("  CONTACTS : %4d => %.3f \n", c, totalContactsAt/(double)totalCounts);
    }


    // At end of each temp, update a probability model for volume?  Use this to select
    // perform positional refinement until delta E stabilizes?
    // final round to move any points close to body
    float volume_test = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);
    // pModel->updateBeadIndices(workingLimit, deadLimit, bead_indices);
    // pModel->setBeadAverageAndStdev(oldN, oldStdev);
    pData->printKLDivergence(binCount);

    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");
    //pModel->writeModelToFile(totalBeadsInSphere, bead_indices, "sphere");

//    string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_final_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);
//    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_" + nameTo, this, pData, numberOfCoolingTempSteps, volume_test, tempAverageContacts);

    //pModel->writeModelToFile(workingLimit, bead_indices, "refined_");
    //pModel->writeModelToFile(deadLimit, bead_indices, "hull_");

    //  pModel->writeModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_" + filenameprefix);
    //  pModel->writeSymModelToFile(lowestWorkingLimit, lowest_bead_indices, "best_lowest_sym_");

    return nameOfModel;
}
/**
 * find an object centered at origin, specify both translation vector and tilt angle
 * create initial model of complete object
 * does not require pre-computed distances, this is a direct method meaning all pairwise distances are calcualted on the fly
 */
bool Anneal::createInitialModelSymmetryEx(Model *pModel, Data *pData) {

    //highTempRounds*=1.62;
    srand(time(0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> distribution(0.0,1.0);

    auto lowTempStop =  (double)0.1;
    double inv_kb_temp = 1.0f/lowTempStop;
    float acceptRate = 0.5f, inv500 = 1.0f/500.0f;
    float inv500slash499 = 499.0f/500.0f;

    contactCutOff = interconnectivityCutOff;
    //violation_limit = 1.8*pModel->getBeadRadius();
    //violation_limit = pModel->getBeadRadius();
    //
    //
    //violation_limit = sqrt(3)*pModel->getBeadRadius();
    violation_limit = (float)(2.0/3.0*(sqrt(3)*pModel->getBeadRadius()));

    maxbin = pModel->getMaxBin(pData);
    // adjust this part
    maxbin += 7;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two

    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50
    std::vector<unsigned int> testBinCount(maxbin);    // smallish vector, typically < 50
    std::vector<unsigned int> workingBinCount(maxbin); // smallish vector, typically < 50
    std::vector<unsigned int> binCountBackUp(maxbin);  // smallish vector, typically < 50

    std::cout << "      TOTAL EXP N_S BINS : " << totalBins << std::endl;
    std::cout << "      MAX MODEL N_S BINS : " << maxbin << std::endl;
    std::cout << "                BINWIDTH : " << pData->getBinWidth() << std::endl;
    std::cout << "             BEAD RADIUS : " << pModel->getBeadRadius() << std::endl;

    unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();

    // as bead indices are discarded, set upper limit of vector
    std::vector<unsigned int> subUnit_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<unsigned int> test_indices(totalBeadsInSphere);        // large vector ~1000's
    std::vector<unsigned int> lowest_subUnit_indices(totalBeadsInSphere); // large vector ~1000's
    std::vector<unsigned int> active_indices(totalBeadsInSphere);         // large vector ~1000's
    std::vector<unsigned int> backUpState(totalBeadsInSphere);

    const unsigned int num = totalBeadsInSphere;
    unsigned int * const ptr = (num != 0) ? subUnit_indices.data() : nullptr;
    for(unsigned int i = 0; i < num; i++) {
        ptr[i] = i;
    }

//    std::set<unsigned int> universe(subUnit_indices.begin(), subUnit_indices.end());
//    std::string universename = "uni";
//    pModel->writeSetToFile(universe, universename);
//    universe.clear();

    lowerV = (unsigned int)(pData->getVolume()  - pData->getVolume() *0.21);
    upperV = (unsigned int)(pData->getVolume()  + pData->getVolume() *0.05);

    const double targetSubUnitVolume = (lowerV+upperV)*0.5/(float)pModel->getNumberOfSubUnits();
    const auto radius_larger = (float)(pModel->getBeadRadius()*std::sqrt(7.0)/2.0);
    const auto lowerN = (unsigned int)(std::round(lowerV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()))/3;
    const auto upperN = (unsigned int)(std::round(upperV/(4.0/3.0*M_PI*radius_larger*radius_larger*radius_larger)/(float)pModel->getNumberOfSubUnits()));

    std::uniform_int_distribution<unsigned int> number_of_beads_to_use (lowerN, upperN);   // number of beads in ASU
    unsigned int subUnitWorkingLimit = number_of_beads_to_use(gen);

    std::cout << "  Bead Search Limited to : " << lowerN << " <= N <= " << upperN << std::endl;
    std::cout << "        INITIAL MODEL WL : " << subUnitWorkingLimit << std::endl;
    std::cout << "                SYMMETRY : " << pModel->getSymmetryString() << std::endl;
    std::cout << "          TOTAL SUBUNITS : " << pModel->getNumberOfSubUnits() << std::endl;

    // randomize and take the workingLength as first set
    // shuffle beads in asymmetric unit
    std::cout << "        CREATING INITIAL RANDOM MODEL " << std::endl;
    std::shuffle(subUnit_indices.begin(), subUnit_indices.end(), gen);
    std::sort(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());

    std::set<unsigned int> beads_in_use_tree(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit);
    std::set<unsigned int> hull;  // sort hull points into bead_indices
    recalculateDeadLimit(subUnitWorkingLimit, subUnit_indices, hull, pModel, totalBeadsInSphere);
    std::string switched = "HULL";
    // randomly pick from each selected index to create monomer
    // each index maps back to beads in universe and also the symmetry grouping
    // calculate volume subunit
    coordT points[3*(upperN+10)];
    char flags[25];
    std::sprintf(flags, "qhull s FA");

    float subunit_test_volume, subunit_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel);
    // calculate starting energy
    unsigned int totalViolations, temp_violations;
    std::uniform_int_distribution<unsigned int> randomIndexInUniverse(0,totalBeadsInSphere-1); // guaranteed unbiased
    double tiltAngle=0, rotationAngle=0;
    const double deltar = pModel->getBeadRadius();
    const unsigned int upperLimit = upperN*3;
    vector3 tempTranslationVector;
    translationVector = vector3(pModel->getBead(randomIndexInUniverse(gen))->getVec().x,0,0);
    std::vector<vector3> coordsInUse(upperLimit), translatedCoordsInUse(upperLimit), tempCoordsInUse(upperLimit), backedUpTranslatedCoordinates(upperLimit), backedUpCoordinates(upperLimit);

    for(unsigned int ind =0; ind<subUnitWorkingLimit; ind++){
        coordsInUse[ind] = pModel->getBead(ind)->getVec();
        translatedCoordsInUse[ind] = pModel->getBead(ind)->getVec() + translationVector;
    }
    std::copy(translatedCoordsInUse.begin(), translatedCoordsInUse.end(), backedUpTranslatedCoordinates.begin() );
    std::copy(coordsInUse.begin(), coordsInUse.end(), backedUpCoordinates.begin() );

    std::vector<double> contactsDistributionOfModel(13);
    std::vector<double> tempContactsDistributionOfModel(13);

    calculateModelPrDistributionSymEX(translatedCoordsInUse, &binCount, contactsDistributionOfModel, subUnitWorkingLimit, pModel, pData, totalViolations);
//    calculateModelPrDistributionSym(&subUnit_indices, &binCount, subUnitWorkingLimit, violations, pModel, pData );

    // fill binCount for first time
    double currentKL = pData->getScore(binCount);

    float hlambda = pData->getIsIntensity() ? 10.0f : 0.01f;
    float muConstant = pData->getIsPr() ? 0.00001f : 0.1f; // chi tends to start in 100's whereas DKL is 0.1;

//    const float targetVolume = (0.5f*(lowerV+upperV));///(float)pModel->getNumberOfSubUnits();
    const auto lowerSubUnitV = (unsigned int)((pData->getVolume()  - pData->getVolume() *0.21)/(float)pModel->getNumberOfSubUnits());
    const auto upperSubUnitV = (unsigned int)((pData->getVolume()  + pData->getVolume() *0.05)/(float)pModel->getNumberOfSubUnits());

    const double invTarVolSubUnitmuConstant = muConstant/targetSubUnitVolume;
    unsigned int failures=0;

    bool updateCVX = false;
    // float muConstant = 0.0001;//mu*0.00001f/((1.0f-mu)*upperV/(float)pModel->getNumberOfSubUnits());
    // pModel->writeModelToFile(groupWorkingLimit, subUnitIndices, "symStart");

    double testKL, test_energy; // sets alpha as a constant during high temp eq

    unsigned int lowestWorkingLimit = subUnitWorkingLimit;
    unsigned int tempNumberOfComponents, currentNumberOfComponents;

    EulerTour eulerTour(subUnit_indices.begin(), subUnitWorkingLimit, pModel);
    currentNumberOfComponents = eulerTour.getNumberOfComponents();

    //double current_energy = currentKL + hlambda*(currentNumberOfComponents-1)*(currentNumberOfComponents-1) ;
    double current_energy = currentKL + lambda*connectivityPotential(currentNumberOfComponents);
    std::cout << "       INITIAL => ENERGY :  " << current_energy << std::endl;
    std::cout << "       INITIAL =>   D_KL : " <<  currentKL << std::endl;
    std::cout << "       INITIAL => VOLUME : " << subunit_test_volume << std::endl;
    std::cout << "              VIOLATIONS : " << totalViolations << std::endl;
    std::cout << "       INITIAL        WL : " << subUnitWorkingLimit << std::endl;

    char addRemoveText[50];
    float sum_x_squared=0;
    unsigned int volumeCount=0, testIndex;
    float volumeSum=0, workingLimitSum=0;


    // what is tolerable number of violations per subsunit?
    current_energy += beta * totalViolations;

    double lowest_energy = current_energy;

    std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin());  // unaltered P(r)
    std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup cop

    const unsigned int noNeigborIndex = pModel->getNeighborLimit();
    unsigned int high, hullsize=hull.size();
    std::uniform_int_distribution<unsigned int> randomHull(0,hullsize-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> randomIndex(0,subUnitWorkingLimit-1); // guaranteed unbiased
    std::uniform_int_distribution<unsigned int> sign(0,1); // guaranteed unbiased

    double temp_subunit_volume_energy, current_subunit_volume_energy = std::abs(subunit_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;
    std::clock_t start = std::clock();
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    std::vector<unsigned int>neighbors(totalNeighbors);
    double deltaAngle = (2.0/180.0*M_PI); // 2 degrees?
    double cosx=1, cosz=1, sinx=0, sinz=0, tcos, tsin;

    for (high=0; high < highTempRounds; high++){ // iterations during the high temp search

        if (distribution(gen) < 0.13161361 ){
            // create points from workSet to determine HULL
            std::sprintf(addRemoveText, "POSITIONAL");
            // setup parameters for hull

            unsigned int * const pSelections = subUnit_indices.data(); // initialized as empty in Model class
            for (unsigned int i = 0; i < subUnitWorkingLimit; i++) {
                beadToPoint(&points[i*3], pModel->getBead(pSelections[i]));
                active_indices[i] = pSelections[i];
            }

            // needs to be optimized
            qh_new_qhull(3, subUnitWorkingLimit, points, 0, flags, nullptr, nullptr);
            vertexT * vertices = qh vertex_list;
            auto totalV = (unsigned int)(qh num_vertices);

            // only move CVX hull points
            std::vector<unsigned int> indices_to_check(totalV); // large vector ~1000's
            for (unsigned int v=0; v < totalV; v++) { //
                indices_to_check[v] = active_indices[qh_pointid(vertices->point)];
                vertices = vertices->next;
            }

            qh_freeqhull(true);

            std::uniform_int_distribution<unsigned int> randomVertices(0, totalV-1);
            unsigned int swap1 = indices_to_check[randomVertices(gen)];
            auto itIndex = std::find(subUnit_indices.begin(), subUnit_indices.begin() + subUnitWorkingLimit, swap1);

            unsigned int locale = std::distance(subUnit_indices.begin(), itIndex);

            // find bead to swap in active set
            eulerTour.removeNode(swap1);
            // find position to swap to
            auto set_it=hull.begin();
            std::advance(set_it, randomHull(gen));
            unsigned int neighbor = *set_it;
            while (neighbor == swap1 || beads_in_use_tree.find(neighbor) != beads_in_use_tree.end()) { // find a new position
                set_it=hull.begin();
                std::advance(set_it, randomHull(gen));
                neighbor = *set_it;
            }

            vector3 swap1Vec = translatedCoordsInUse[locale]; //
            // rotate and translate
            translatedCoordsInUse[locale] = rotate((pModel->getBead(neighbor)->getVec()), cosx, sinx, cosz, sinz) + translationVector;

            // make the swap, sort and update P(r)
            auto pSwap2 = std::find(subUnit_indices.begin()+subUnitWorkingLimit, subUnit_indices.end(), neighbor);
            std::iter_swap(itIndex, pSwap2);

            calculateModelPrDistributionSymEX(translatedCoordsInUse, &binCount, tempContactsDistributionOfModel, subUnitWorkingLimit, pModel, pData, temp_violations);

            testKL = pData->getScore(binCount); // 100x faster than calculateKLEnergySymmetry
            tempNumberOfComponents = eulerTour.addNode(neighbor, pModel);

            subunit_test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel); // calculate volume of the subUnit
            temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

            test_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations;

            if (( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy-temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                coordsInUse[locale] = pModel->getBead(neighbor)->getVec(); // sync untransformed vector
                beads_in_use_tree.erase(swap1);
                beads_in_use_tree.insert(neighbor);
                updateCVX = true;
            } else {
                translatedCoordsInUse[locale] = swap1Vec;
                std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                std::copy(backUpState.begin(), backUpState.end(), subUnit_indices.begin());
                eulerTour.removeNode(neighbor);
                eulerTour.addNode(swap1, pModel);
            }

        } else { // add or remove

            if (distribution(gen) < 0.09131 && subUnitWorkingLimit < upperN){

                if (distribution(gen) < 0.5){ // translate
                    std::sprintf(addRemoveText, " TRANSLATE");

                    double plusminus = (sign(gen) == 1) ? 1.0f : -1.0f;
                    tempTranslationVector = vector3(deltar*plusminus, 0, 0) + translationVector;

                    // translate vector and recalculate DKL
                    for(unsigned int i=0; i<subUnitWorkingLimit; i++){
                        tempCoordsInUse[i] = rotate(coordsInUse[i], cosx, sinx, cosz, sinz) + tempTranslationVector;
                        //tempCoordsInUse[i] =  coordsInUse[i] + tempTranslationVector;
                    }
                    calculateModelPrDistributionSymEX(tempCoordsInUse, &binCount, tempContactsDistributionOfModel, subUnitWorkingLimit, pModel, pData, temp_violations);

                    testKL = pData->getScore(binCount); // 100x faster than calculateKLEnergySymmetry
                    test_energy = testKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*temp_violations;

                    if (test_energy < current_energy ||  (std::exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))){
                        tempNumberOfComponents = currentNumberOfComponents;
                        subunit_test_volume = subunit_volume;
                        temp_subunit_volume_energy = current_subunit_volume_energy;

                        std::copy(tempCoordsInUse.begin(), tempCoordsInUse.end(), translatedCoordsInUse.begin());
                        translationVector = tempTranslationVector;
                        updateCVX=true;
                    } else {
                        std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    }

                } else if (distribution(gen) < 0.5) { //rotate
                    std::sprintf(addRemoveText, "  ROTATE  ");
                    double plusminus = (sign(gen) == 1) ? 1.0f : -1.0f;
                    double newAngle = rotationAngle + deltaAngle*plusminus;

//                    if (newAngle < 0){
//                        double trotationAngle = newAngle < 0 ? (2*M_PI - newAngle) : newAngle;
//                        std::cout << "minus " << std::cos(newAngle) << " " << std::cos(trotationAngle) << std::endl;
//                        exit(0);
//                    }

                    // rotate coordinates in Use
                    tcos = std::cos(newAngle);
                    tsin = std::sin(newAngle);

                    for(unsigned int i=0; i<subUnitWorkingLimit; i++){
                        const auto & pVec = coordsInUse[i];
                        tempCoordsInUse[i] =  vector3((float)(pVec.x*tcos - pVec.y*tsin), (float)(pVec.x*tsin + pVec.y*tcos), pVec.z);
                    }

                    calculateModelPrDistributionSymEX(tempCoordsInUse, &binCount, tempContactsDistributionOfModel, subUnitWorkingLimit, pModel, pData, temp_violations);
                    testKL = pData->getScore(binCount); // 100x faster than calculateKLEnergySymmetry
                    test_energy = testKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*temp_violations;

                    if (test_energy < current_energy ||  (std::exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))){
                        tempNumberOfComponents = currentNumberOfComponents;
                        subunit_test_volume = subunit_volume;
                        temp_subunit_volume_energy = current_subunit_volume_energy;

                        std::copy(tempCoordsInUse.begin(), tempCoordsInUse.end(), translatedCoordsInUse.begin());
                        rotationAngle = newAngle < 0 ? (2*M_PI - newAngle) : newAngle;
                        cosz = tcos;
                        sinz = tsin;
                        updateCVX=true;
                    } else {
                        std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    }
                } else { //tilt (rotate along X)
                    std::sprintf(addRemoveText, "  TILT    ");
                    double plusminus = (sign(gen) == 1) ? 1.0f : -1.0f;
                    double newAngle = tiltAngle + deltaAngle*plusminus;
                    // rotate coordinates in Use
                    tcos = std::cos(newAngle);
                    tsin = std::sin(newAngle);

                    for(unsigned int i=0; i<subUnitWorkingLimit; i++){
                        const auto & pVec = coordsInUse[i];
                        tempCoordsInUse[i] =  vector3(pVec.x, pVec.y*tcos - pVec.z*tsin, pVec.y*tsin + pVec.z*tcos);
                    }

                    calculateModelPrDistributionSymEX(tempCoordsInUse, &binCount, tempContactsDistributionOfModel, subUnitWorkingLimit, pModel, pData, temp_violations);
                    testKL = pData->getScore(binCount); // 100x faster than calculateKLEnergySymmetry
                    test_energy = testKL + lambda*connectivityPotential(currentNumberOfComponents) + beta*temp_violations;

                    if (test_energy < current_energy ||  (std::exp((current_energy - test_energy) * inv_kb_temp) > distribution(gen))){

                        tempNumberOfComponents = currentNumberOfComponents;
                        subunit_test_volume = subunit_volume;
                        temp_subunit_volume_energy = current_subunit_volume_energy;

                        std::copy(tempCoordsInUse.begin(), tempCoordsInUse.end(), translatedCoordsInUse.begin());
                        tiltAngle = newAngle < 0 ? (2*M_PI - newAngle) : newAngle;
                        cosx = tcos;
                        sinx = tsin;
                        updateCVX=true;
                    } else {
                        std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin());
                    }
                }

            } else {
                // hard limit on number of beads
                if ( ((distribution(gen) < 0.5) && subUnitWorkingLimit > lowerN) || subUnitWorkingLimit > upperN ) { // REMOVE beads from sorted list into useable range < deadLimit
                    // randomly swap positions with end of workingLength, could remove CVX Hull Point
                    std::sprintf(addRemoveText, "  REMOVE  ");

                    unsigned int locale = randomIndex(gen);
                    testIndex = subUnit_indices[locale];

                    subUnitWorkingLimit -=1;
                    std::iter_swap(coordsInUse.begin() + locale, coordsInUse.begin() + subUnitWorkingLimit);
                    std::iter_swap(translatedCoordsInUse.begin() + locale, translatedCoordsInUse.begin() + subUnitWorkingLimit);
                    std::iter_swap(subUnit_indices.begin()+locale, subUnit_indices.begin() + subUnitWorkingLimit);

                    calculateModelPrDistributionSymEX(translatedCoordsInUse, &binCount, tempContactsDistributionOfModel, subUnitWorkingLimit, pModel, pData, temp_violations);
                    testKL = pData->getScore(binCount);

                    subunit_test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel); // calculate volume of the subUnit
                    temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                    tempNumberOfComponents = eulerTour.removeNode(testIndex);

                    test_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations;

                    if (( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy - temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                        beads_in_use_tree.erase(testIndex);
                        updateCVX=true;
                    } else { // undo changes and move to next bead (rejecting)
                        std::copy(backedUpTranslatedCoordinates.begin(), backedUpTranslatedCoordinates.end(), translatedCoordsInUse.begin() );
                        std::copy(backedUpCoordinates.begin(), backedUpCoordinates.end(), coordsInUse.begin() );
                        eulerTour.addNode(testIndex, pModel);
                        subUnitWorkingLimit += 1;
                        std::copy(binCountBackUp.begin(), binCountBackUp.end(), binCount.begin()); //copy to bin count
                        std::copy(backUpState.begin(), backUpState.end(), subUnit_indices.begin());
                    }

                } else { // ADD beads
                    std::sprintf(addRemoveText, "    ADD   ");

                    std::set<unsigned int>::const_iterator set_it=hull.begin();
                    std::advance(set_it, randomHull(gen));
                    testIndex = *set_it;

                    while ( beads_in_use_tree.find(testIndex) != beads_in_use_tree.end()) { // find a new position
                        set_it=hull.begin();
                        std::advance(set_it, randomHull(gen));
                        testIndex = *set_it;
                    }

//                testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
//                while ( testIndex == noNeigborIndex) { // find a new position
//                    testIndex = getUseableNeighborFromSet(&beads_in_use_tree, pModel, subUnit_indices[randomIndex(gen)]);
//                }

                    auto itIndex = std::find(subUnit_indices.begin(), subUnit_indices.end(), testIndex);
                    // new vector gets added on to end
                    coordsInUse[subUnitWorkingLimit] = pModel->getBead(testIndex)->getVec();
                    translatedCoordsInUse[subUnitWorkingLimit] = rotate((pModel->getBead(testIndex)->getVec()), cosx, sinx, cosz, sinz) + translationVector;

                    std::iter_swap(subUnit_indices.begin() + subUnitWorkingLimit, itIndex); // this swaps to working position and changes backup
                    // increment workingLimit to include new position
                    subUnitWorkingLimit += 1;

                    calculateModelPrDistributionSymEX(translatedCoordsInUse, &binCount, tempContactsDistributionOfModel, subUnitWorkingLimit, pModel, pData, temp_violations);
                    testKL = pData->getScore(binCount);

                    subunit_test_volume = calculateCVXHULLVolume(flags, &subUnit_indices, subUnitWorkingLimit, pModel); // calculate volume of the subUnit
                    temp_subunit_volume_energy = std::abs(subunit_test_volume-targetSubUnitVolume)*invTarVolSubUnitmuConstant;

                    // waste of time since connectivity will remain constant?  No, if connectivity is 2, I could add a position that changes to 1
                    tempNumberOfComponents = eulerTour.addNode(testIndex, pModel);

                    test_energy = testKL + lambda*connectivityPotential(tempNumberOfComponents) + beta*temp_violations;

                    if (( (test_energy + temp_subunit_volume_energy) < (current_energy + current_subunit_volume_energy)) || (std::exp((current_energy - test_energy + (current_subunit_volume_energy-temp_subunit_volume_energy)) * inv_kb_temp) > distribution(gen)) ) {
                        beads_in_use_tree.insert(testIndex);
                        updateCVX=true;
                    } else { // undo changes (rejecting)
                        auto beginBinCount = binCount.begin();
                        restoreAddingFromBackUp(&subUnit_indices, &backUpState, &subUnitWorkingLimit, &binCountBackUp, &beginBinCount);
                        eulerTour.removeNode(testIndex);
                    }
                }
            }
        } // END OF ADD?REMOVE


        if (updateCVX){

            current_energy = test_energy;
            currentKL = testKL;
            currentNumberOfComponents = tempNumberOfComponents;
            totalViolations = temp_violations;
            subunit_volume = subunit_test_volume;
            current_subunit_volume_energy = temp_subunit_volume_energy;

            acceptRate = inv500slash499*acceptRate+inv500;
            updateCVX = false;
            failures=0;

//            if ( hullsize < 1.3*subUnitWorkingLimit || currentNumberOfComponents == 1){
                populateLayeredDeadlimitUsingSet(beads_in_use_tree, hull, pModel);
                switched = "NEAREST NEIGHBOR";
//            } else {
//            switched = "HULL";
//            recalculateDeadLimit(subUnitWorkingLimit, subUnit_indices, hull, pModel, totalBeadsInSphere);
//            }

            hullsize = hull.size();
            randomHull = std::uniform_int_distribution<unsigned int>(0,hullsize-1);

            std::copy(subUnit_indices.begin(), subUnit_indices.end(), backUpState.begin());   // make backup copy
            std::copy(translatedCoordsInUse.begin(), translatedCoordsInUse.end(), backedUpTranslatedCoordinates.begin() );
            std::copy(coordsInUse.begin(), coordsInUse.end(), backedUpCoordinates.begin() );
            std::copy(binCount.begin(), binCount.end(), binCountBackUp.begin()); // make backup copy
            std::copy(tempContactsDistributionOfModel.begin(), tempContactsDistributionOfModel.end(), contactsDistributionOfModel.begin());
            randomIndex = std::uniform_int_distribution<unsigned int>(0,subUnitWorkingLimit-1); // guaranteed unbiased

        } else {
            acceptRate = inv500slash499*acceptRate;
            failures++;
        }

        updateASAConstantTemp(high, highTempRounds, acceptRate, lowTempStop, inv_kb_temp);

        // if no positional refinement using populate, the deadlimit space is fine
        // if using positional refinement with populate, get artefacts

        std::printf("*******************             %s                 ******************* \n", addRemoveText);
        std::printf("      MAXSTEPS : %7i (%7i) ACCEPTRATE : %.3f TEMP : %.2E\n", highTempRounds, high, acceptRate, lowTempStop);
        std::printf("         GRAPH : %7i  MODE => %s \n", currentNumberOfComponents, switched.c_str());
        std::printf("         LIMIT : %7i UPPER <= %i  LOWER >= %i (%i)\n", subUnitWorkingLimit, upperN, lowerN, hullsize);
        std::printf("        VOLUME : %.0f (%.0f)  MU => %.1E  MU*VOL => %.2E\n", subunit_volume, targetSubUnitVolume, muConstant, current_subunit_volume_energy);
        std::printf("   ORIENTATION : TILT %.4f YAW %.4f TRANS %.2f \n", tiltAngle, rotationAngle, translationVector.length());
        std::printf("    VIOLATIONS : %7d  BETA => %.1E (%i)  \n", totalViolations, beta, volumeCount);
        std::printf("          D_KL : %.4E ENRGY: %.4E (%.4E) \n", currentKL, current_energy, lowest_energy);
        std::cout << "*******************                                        *******************" << std::endl;

        // write to file to animate search
        if (currentNumberOfComponents == 1 && current_energy < lowest_energy ){
            std::copy(subUnit_indices.begin(), subUnit_indices.end(), lowest_subUnit_indices.begin());
            lowestWorkingLimit = subUnitWorkingLimit;

            lowest_energy = current_energy;
            workingLimitSum += subUnitWorkingLimit;
            sum_x_squared += subUnitWorkingLimit*subUnitWorkingLimit;
            volumeSum += subunit_volume;
            volumeCount++;
        }

//        checkSetAndVector(subUnitWorkingLimit, &subUnit_indices, &beads_in_use_tree);

//        if (checkForRepeats(subUnit_indices)){
//            std::cout << " STOPPED POSITIONAL " << " WL: " << subUnitWorkingLimit << " D_KL " << currentKL << " <=> " << std::endl;
//            exit(0);
//        }

//        violations = getViolations(&subUnit_indices, subUnitWorkingLimit, pModel, pData);
//        if (violations != totalViolations){
//            std::cout << " totalViolations " << totalViolations << " actual " << violations << std::endl;
//            testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, violations, pModel, pData );
//            if (currentKL != testKL){
//                printf("*******************             %s                 ******************* \n", addRemoveText);
//                std::cout << high  << "           Failed " << pData->getScore(binCount) << " != " << testKL << std::endl;
//                std::cout << high  << " Failed currentKL " << currentKL << " != " << testKL << std::endl;
//                printf("CKL %.8E  TKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//            }
//            exit(0);
//        }
//
//        calculateModelPrDistributionSym(&subUnit_indices, &testBinCount, subUnitWorkingLimit, violations, pModel, pData );
//        testKL = pData->getScore(testBinCount); // 100x faster than calculateKLEnergySymmetry
//        if (currentKL != testKL){
//            printf("*******************             %s                 ******************* \n", addRemoveText);
//            std::cout << high  << "           Failed " << pData->getScore(binCount) << " != " << testKL << std::endl;
//            std::cout << high  << " Failed currentKL " << currentKL << " != " << testKL << std::endl;
//            printf("CKL %.8E  TKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//            exit(0);
//        }


//        testKL = calculateKLEnergySymmetry(&subUnit_indices, &testBinCount, subUnitWorkingLimit, violations, pModel, pData );
//        if (currentKL != testKL){
//            // if this fails, update of P(r) is wrong or subUnit_indices is corrupted
//            printf("*******************             %s                 ******************* \n", addRemoveText);
//            std::cout << high  << "           Failed " << pData->getScore(binCount) << " != " << testKL << std::endl;
//            std::cout << high  << " Failed currentKL " << currentKL << " != " << testKL << std::endl;
//            printf("CKL %.8E  TKL %.8E  DIFF => %.8E\n", currentKL, testKL, (testKL - currentKL));
//
//            int bins = testBinCount.size();
//            std::cout << " BIN : " << std::endl;
//            for(int i=0; i<bins; i++){
//                if (testBinCount[i] != binCount[i]){
//                    std::cout << i << " : " << testBinCount[i] << " == " << binCount[i] << std::endl;
//                }
//            }
//
//            exit(0);
//        }

    } // end of HIGH TEMP EQUILIBRATION


    switched = "HULLFinal";
    recalculateDeadLimit(lowestWorkingLimit, lowest_subUnit_indices, hull, pModel, totalBeadsInSphere);
    pModel->writeSetToFile(hull, switched);



    highTempStartForCooling = (float)lowTempStop;
    // pModel->updateBeadIndices(lowestWorkingLimit, lowestDeadLimit, lowest_subUnit_indices);
    // calculate average volume and standard deviation
    float volumeAverage = workingLimitSum/(float)volumeCount;
    float volumeStdev = std::sqrt(sum_x_squared/(float)volumeCount - volumeAverage*volumeAverage);

    // remove points close to hull

    pModel->setStartingSet(lowest_subUnit_indices);
    pModel->setStartingWorkingLimit(lowestWorkingLimit);
    pModel->setBeadAverageAndStdev(volumeAverage, volumeStdev);

    std::string tempName = "initial_CVX_subunit_" + filenameprefix;
    pModel->writeModelToFile(lowestWorkingLimit, lowest_subUnit_indices, tempName , high);
    pModel->writeSymModelToFile(currentKL, lowestWorkingLimit, lowest_subUnit_indices, binCount, "initial_CVX_sym_" + filenameprefix, this, pData, high, subunit_volume, 4.1);



    std::cout << "AVERAGE # BEADS EST HIGH TEMP SELECTION: " << (int)volumeAverage << " SIGMA: " << (int)volumeStdev << std::endl;
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    //std::printf("   AVERAGE => %0.f (%0.f) \n", pModel->getVolumeAverage(), pModel->getVolumeStdev());
    std::printf("    LOWEST => %d \n", lowestWorkingLimit);
    std::cout << "*******************                                        *******************" << std::endl;

    EulerTour finalEulerTour(lowest_subUnit_indices.begin(), lowestWorkingLimit, pModel);
    currentNumberOfComponents = finalEulerTour.getNumberOfComponents();

    if (currentNumberOfComponents == 1 ) {
        return true;
    } else {
        tempName = "failed_subunit_" + filenameprefix;
        pModel->writeModelToFile(subUnitWorkingLimit, subUnit_indices, tempName , high);
        pModel->writeSymModelToFile(currentKL, subUnitWorkingLimit, subUnit_indices, binCount, "failed_sym" + filenameprefix, this, pData, high, subunit_volume, 4.1);

        std::cout << "SEARCH TOO SHORT, EULER TOUR > 1 " << std::endl;
        std::cout << "INCREASE highTempRounds, g" << std::endl;
        return false;
    }
}


unsigned int Anneal::minimizeViolations(std::vector<unsigned int> *bead_indices, unsigned int indicesWorkingLimit, Model *pModel){
    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;
    Bead * tempBead;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(bead_indices->begin(), bead_indices->begin()+indicesWorkingLimit, gen);

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        tempBead = pModel->getBead((*bead_indices)[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }


    std::set<unsigned int> violated;
    float distance_to;
    unsigned int violations = 0;


    for(unsigned int i=0; i<indicesWorkingLimit; i++) {

        const vector3 *baseVec = &coordinates[i];
        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
            distance_to = ((*baseVec) - coordinates[s]).length();
            if (distance_to < violation_limit){
                violated.insert((*bead_indices)[i]);
                violations++;
            }
        }
    }

    unsigned int swapTo = indicesWorkingLimit;

    for(auto vio : violated){
        auto it = std::find(bead_indices->begin(), bead_indices->begin()+swapTo, vio);
        if (it != (bead_indices->begin()+swapTo)){
            swapTo-=1;
            std::iter_swap(it, bead_indices->begin()+swapTo);
        }
    }

    return swapTo;
}

void Anneal::printSymModel(std::vector<unsigned int> *subUnit_indices, const unsigned int indicesWorkingLimit, Model *pModel, Data *pData){
    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;
    Bead * tempBead;

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        tempBead = pModel->getBead((*subUnit_indices)[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    std::vector<std::string> chains(25);
    chains[0] = "A";
    chains[1] = "B";
    chains[2] = "C";
    chains[3] = "D";
    chains[4] = "E";
    chains[5] = "F";
    chains[6] = "G";
    chains[7] = "H";
    chains[8] = "I";
    chains[9] = "J";
    chains[10] = "K";
    chains[11] = "L";
    chains[12] = "M";
    chains[13] = "N";
    chains[14] = "O";
    chains[15] = "P";
    chains[16] = "R";
    chains[17] = "S";
    chains[18] = "T";
    chains[19] = "U";
    chains[20] = "V";
    chains[21] = "W";
    chains[22] = "X";
    chains[23] = "Y";
    chains[24] = "Q";

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

    std::string nameOf = "symModel.pdb";
    FILE * pFile = fopen(nameOf.c_str(), "w");
    std::string chainID = chains[0];

    unsigned int residueCnt=1;
    for(unsigned int i=0; i<totalCoordinates; i++) {

        if (i%indicesWorkingLimit == 0 && i > 0){ // reinitialize counter
            residueCnt = 1;

            unsigned int locale = std::floor((double)i/(double)indicesWorkingLimit);
            locale = (locale >= 25) ? 25-locale : locale;
            if (locale > 0){
                chainID = chains[locale];
            }
        }

        pModel->printAtomLine(pFile, residueCnt, chainID, residueCnt, coordinates[i].x, coordinates[i].y, coordinates[i].z);
    }

    fclose(pFile);

}
/**
 * count the number of violations of parent subunit
 * @param bead_indices
 * @param beadIndiciesWorkingLimit
 * @param pModel
 * @param pData
 * @return
 */
unsigned int Anneal::getViolations(std::vector<unsigned int> *subUnit_indices, const unsigned int indicesWorkingLimit, Model *pModel){


    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;
    Bead * tempBead;

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        tempBead = pModel->getBead((*subUnit_indices)[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

//    std::string nameOf = "violations.pdb";
//    FILE * pFile = fopen(nameOf.c_str(), "a");

    float distance_to;
    unsigned int violations = 0;
    // count violations
    for(unsigned int i=0; i<indicesWorkingLimit; i++) {

        const vector3 *baseVec = &coordinates[i];
        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
            distance_to = ((*baseVec) - coordinates[s]).length();
            if (distance_to < violation_limit){
//                pModel->printAtomLine(pFile, s, "A", s, coordinates[s].x, coordinates[s].y, coordinates[s].z);
                violations++;
            }
        }
    }

//    fclose(pFile);
    // check for connectivity
    return violations;
}


/**
 * Calculate both the number of violations (clashes) from symmetry mates and the contact distirbution ofthe subunit
 * @param subUnit_indices (sorted due to SET container)
 * @param contactsDistributionOfModel
 * @param indicesWorkingLimit
 * @param pModel
 * @return
 */
unsigned int Anneal::getViolationsFromSet(std::set<unsigned int> *subUnit_indices, std::vector<double> & contactsDistributionOfModel, const unsigned int indicesWorkingLimit, Model *pModel){


    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);


    // create first subunit from selected indices and populate coordinates
    vector3 * const ptr = coordinates.data();
    unsigned int index=0;
    for(auto & ind : *subUnit_indices){ // sorted in order because of set
        ptr[index] = pModel->getBead(ind)->getVec();
        index++;
    }

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }


//    std::string nameOf = "violations.pdb";
//    FILE * pFile = fopen(nameOf.c_str(), "a");

    float distance_to;
    unsigned int violations = 0;
    std::fill(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), 0.0d);
    // count violations
    index=0;
    for(auto & ind : *subUnit_indices){
        unsigned int tempContacts = numberOfContactsFromSet(subUnit_indices, pModel, ind);
        const vector3 *baseVec = &ptr[index];

        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
            distance_to = ((*baseVec) - ptr[s]).length();

 //           if (distance_to < contactCutOff){ // this includes additional symmetry mate contacts
 //               tempContacts++;
                if (distance_to < violation_limit){
//                    pModel->printAtomLine(pFile, s, "A", s, ptr[s].x, ptr[s].y, ptr[s].z);
                    violations++;
                }
 //           }
        }
        index++;
 //       tempContacts = (tempContacts > 12) ? 12 : tempContacts;
        ++contactsDistributionOfModel[tempContacts];
    }

//    fclose(pFile);
    // entire Symmetry made molecule
//    for(unsigned int index=0;index<totalCoordinates; index++){
//        const vector3 *baseVec = &ptr[index];
//        unsigned int counter=0;
//        for(unsigned int index2=index+1;index2<totalCoordinates; index2++){
//            const vector3 *baseVec2 = &ptr[index2];
//            distance_to = (*baseVec - *baseVec2).length();
//            if (distance_to < contactCutOff){
//                counter++;
//                if (counter >= 12){
//                    break;
//                }
//            }
//        }
//        ++contactsDistributionOfModel[counter];
//    }


//    for(auto & ind : *subUnit_indices){
//        unsigned int tempContacts = numberOfContactsFromSet(subUnit_indices, pModel, ind);
//        const vector3 *baseVec = &ptr[index];
//
//        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
//            distance_to = ((*baseVec) - ptr[s]).length();
//
//            if (distance_to < contactCutOff){
//                tempContacts++;
//                if (distance_to < violation_limit){
//                    violations++;
//                }
//            }
//        }
//        index++;
//
//        tempContacts = (tempContacts > 12) ? 12 : tempContacts;
//        ++contactsDistributionOfModel[tempContacts];
//    }

    // check for connectivity
    return violations;
}


/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 * calculation is from the coordinates
 * Steps:
 * 1. create coordinates of first subunit (workingLimit)
 * 2. create symmetry mates and add to a mster list of coordinates
 * 3. calculate distance distribution for each pair of coordinates
 */
void Anneal::calculateModelPrDistributionSymCE(std::vector<unsigned int> *subUnit_indices,
                                               std::set<unsigned int> & subUnit_indices_tree,
                                               std::map<unsigned int, std::vector<vector3> > &map,
                                               std::vector<double> & contactsDistributionOfModel,
                                               std::vector<unsigned int> *binCount,
                                               const unsigned int indicesWorkingLimit,
                                               unsigned int &violations,
                                               Model *pModel, Data *pData) {

    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0);

    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3 *> coordinates(totalCoordinates);

    vector3 ** const pCoords = coordinates.data();
    unsigned int * const pSubInd = (*subUnit_indices).data();
    // create first subunit from selected indices and populate coordinates
//    for (unsigned int i=0; i < indicesWorkingLimit; i++){
//        tempBead = pModel->getBead(pSubInd[i]);
//        ptr[i] = vector3(tempBead->getX(), tempBead->getY(), tempBead->getZ());
//    }

    unsigned int index=0;
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        auto mit = map.find(pSubInd[i]);
        vector3 * const pVecBase = (*mit).second.data();

        for(unsigned int s=0; s<pModel->getNumberOfSubUnits(); s++){ // iterate through map, add all symmetry related points
            pCoords[index] = &pVecBase[s];
            index++;
        }
    }

//    unsigned int count = indicesWorkingLimit;
//    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
//        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
//    }

    // calculate Pr, order is irrelavant (all pairwise)
    // parallelize this code

    float distance_to;
    for (unsigned int i=0; i < totalCoordinates; i++) {

        const vector3 * tempVec1 = pCoords[i];
        /*
         * each thread to have its own copy of tempVec1
         */

        for (unsigned int next_i = i + 1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = (*tempVec1 - *pCoords[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // unrolling could lead to a RACE condition of updating same location at same time
        }
    }


    /*
     * for each lattice point within parent, determine how many violations it makes with other points in symmetry mates
     */
//    violation=0;
//    for (unsigned int i=0; i < indicesWorkingLimit; i++) {
//
//        const vector3 * ptempVec = &coordinates[i];
//
//        for (unsigned int next_i = indicesWorkingLimit; next_i < totalCoordinates; next_i++) {
//            // calculate distance and convert to bin
//            distance_to = ((*ptempVec) - coordinates[next_i]).length();
//            if (distance_to < violation_limit){
//                violation++;
//            }
//        }
//    }

    std::fill(contactsDistributionOfModel.begin(), contactsDistributionOfModel.end(), 0.0d);
    violations=0;
    unsigned int tempContacts;
    for(unsigned int ind=0; ind<indicesWorkingLimit; ind++){
        tempContacts = numberOfContactsFromSet(&subUnit_indices_tree, pModel, pSubInd[ind]);
        const vector3 *baseVec = pCoords[ind];

        for (unsigned int s = indicesWorkingLimit; s < totalCoordinates; s++) {
            distance_to = (*baseVec - *pCoords[s]).length();
            if (distance_to < contactCutOff){ // since violation_limit is always < contactCutOff, only check if true
                tempContacts++;
                if (distance_to < violation_limit){
                    violations++;
                }
            }
        }
        tempContacts = (tempContacts > 12) ? 12 : tempContacts;
        ++contactsDistributionOfModel[tempContacts];
    }

}


/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 * calculation is from the coordinates
 * Steps:
 * 1. create coordinates of first subunit (workingLimit)
 * 2. create symmetry mates and add to a mster list of coordinates
 * 3. calculate distance distribution for each pair of coordinates
 */
void Anneal::calculateModelPrDistributionSym(std::vector<unsigned int> *subUnit_indices, std::vector<unsigned int> *binCount, const unsigned int indicesWorkingLimit, unsigned int &violation,  Model *pModel, Data *pData) {
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0.0);

    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;
    vector3 * const ptr = coordinates.data();
    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        ptr[i] = pModel->getBead((*subUnit_indices)[i])->getVec();
    }


    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

    // calculate Pr, order is irrelavant
    float distance_to;
    violation=0;
    for (unsigned int i=0; i < indicesWorkingLimit; i++) {

        tempVec1 = &coordinates[i];

        for (unsigned int next_i = i + 1; next_i < indicesWorkingLimit; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }

        for (unsigned int next_i = indicesWorkingLimit; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
            if (distance_to < violation_limit){
                violation++;
            }
        }
    }

    for (unsigned int i=indicesWorkingLimit; i < totalCoordinates; i++) {

        tempVec1 = &coordinates[i];

        for (unsigned int next_i = i + 1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }
}

/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 * calculation is from the coordinates
 * Steps:
 * 1. create coordinates of first subunit (workingLimit)
 * 2. create symmetry mates and add to a mster list of coordinates
 * 3. calculate distance distribution for each pair of coordinates
 */
void Anneal::calculateModelPrDistributionSymEX(std::vector<vector3> & translatedVector, std::vector<unsigned int> *binCount, std::vector<double> & contactsDistributionOfModel, const unsigned int indicesWorkingLimit,  Model *pModel, Data *pData, unsigned int &violations) {
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0.0);

    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;

    vector3 * const ptr = coordinates.data();
    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        ptr[i] = translatedVector[i];
    }

    unsigned int count = indicesWorkingLimit, tempContacts;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

    // calculate Pr, order is irrelavant
    float distance_to;
    violations=0;
    for (unsigned int i=0; i < indicesWorkingLimit; i++) {

        tempVec1 = &coordinates[i];
        tempContacts=0;
        for (unsigned int next_i = i + 1; next_i < indicesWorkingLimit; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
            if (distance_to < contactCutOff){
                tempContacts++;
            }
        }

        //tempContacts = (tempContacts > 12) ? 12 : tempContacts;
        ++contactsDistributionOfModel[tempContacts];

        for (unsigned int next_i = indicesWorkingLimit; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
//            if (distance_to < violation_limit){
//                violations++;
//            }
        }
    }


    for (unsigned int i=indicesWorkingLimit; i < totalCoordinates; i++) {

        tempVec1 = &coordinates[i];

        for (unsigned int next_i = i + 1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }

}

/**
 * The new position to Add MUST BE WITHIN THE WORKING LIMIT
 * beadsInUse must be sorted up to workingLimit
 *
 * @param addMeSubUnitIndex
 * @param beadsInUse
 * @param workingLimit
 * @param prBins
 * @param pModel
 * @param pData
 * @return
 */
inline unsigned int Anneal::addToPrSym(unsigned int addMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData){
    // each bead has a sym related partner in pModel
    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin

    unsigned int findIt=0;
    const unsigned int * const ptr = beadsInUse.data();
    for(unsigned int i=(workingLimit-1); i>=0; i--){
        if (ptr[i] == addMeSubUnitIndex){
            findIt = i;
            break;
        }
    }

    const int stopAt = findIt;
    // add interdomain distances
    //
    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * const pC = coordinates.data();

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < workingLimit; i++){
        pC[i] = pModel->getBead(ptr[i])->getVec();
    }

    unsigned int count = workingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, workingLimit, count, coordinates);
    }

    // adjust Pr
    float distance_to;
    for (unsigned int s=0; s < subUnits; s++){  //
        unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &pC[basis];

        for (unsigned int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - pC[next_i]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }

        for (unsigned int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - pC[next_i]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }

    /*
     * correct for double counting
     */
    for (unsigned int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &pC[basis];

        for (unsigned int ss=s+1; ss < subUnits; ss++){
            distance_to = ((*prefVec) - pC[stopAt+ss*workingLimit]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }
    /*
     * Should include contacts distribution
     *
     *
     * count any symmetry related clashes
     */
    unsigned int totalviolations=0;
    for (unsigned int i=0; i < workingLimit; i++){
        const vector3 * pBase = &pC[i];
        for (unsigned int next=workingLimit; next < totalCoordinates; next++){
            distance_to = ((*pBase) - pC[next]).length();
            if (distance_to < violation_limit){
                ++totalviolations;
            }
        }
    }


    return totalviolations;
}






/**
 * beadsInUse does not have to be sorted
 *
 * prBins is the model P(r) distribution
 *
 * returns the number of violations adjusted with respect to removeMeSubUnitIndex
 *
 *
 */
unsigned int Anneal::removeFromPrSym(unsigned int const &removeMeSubUnitIndex, std::vector<unsigned int> & beadsInUse, unsigned int const &workingLimit, std::vector<unsigned int> & prBins, Model *pModel, Data *pData){

    unsigned int findIt=0;
    // find index of bead to remove
    const unsigned int * const ptr = beadsInUse.data();
    for(unsigned int i=0; i < workingLimit; i++){
        if (ptr[i] == removeMeSubUnitIndex){
            findIt = i;
            break;
        }
    }

    const unsigned int stopAt = findIt;
    //
    // remove interdomain distances
    //
    float distance_to;

    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < workingLimit; i++){
        coordinates[i] = pModel->getBead(ptr[i])->getVec();
    }

    unsigned int count = workingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, workingLimit, count, coordinates);
    }


    // adjust Pr
    for (unsigned int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (unsigned int next_i = 0; next_i < basis; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }

        for (unsigned int next_i = basis+1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*prefVec) - coordinates[next_i]).length();
            --prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }


    /**
     *
     * we double remove(count) the symmetry related vectors, need to add back
     */
    for (unsigned int s=0; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        unsigned int basis = stopAt+s*workingLimit;
        const vector3 * prefVec = &coordinates[basis];

        for (unsigned int ss=s+1; ss < subUnits; ss++){
            distance_to = ((*prefVec) - coordinates[stopAt+ss*workingLimit]).length();
            ++prBins[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }

    return 0;
}


/**
 *
 *
 */
unsigned int Anneal::removeLatticePositionToModelSym(
                                                   std::vector<unsigned int> & bead_indices,
                                                   std::vector<unsigned int> & modelPrBins,
                                                   unsigned int * pWorkingLimit,
                                                   const unsigned int * pLatticePointToRemove, Model * pModel, Data *pData){
    auto pBeginIt = bead_indices.begin();
    auto itIndex = std::find(pBeginIt, pBeginIt + *pWorkingLimit, *pLatticePointToRemove);
    // remove original from P(r)
    // copy(beginBinCount, endBinCount, binCountBackUp.begin()); //copy to bin count
    unsigned int violations = removeFromPrSym(*pLatticePointToRemove, bead_indices, *pWorkingLimit, modelPrBins, pModel, pData);
    // reduce the workingLimit
    // if wl = 10
    // 0 1 2 3 4 5 6 7 8 9 10
    // remove 4
    // 0 1 2 3 9 5 6 7 8 4 10
    //
    // 0 1 2 3 5 6 7 8 9 4 10
    //
    *pWorkingLimit -= 1;
    // swap selected point to the workingLimit
    std::iter_swap(itIndex, pBeginIt + *pWorkingLimit);
    // still need to sort, swap changes the order
    // to save some time, sort should only be from swapped point to workingLimit
    //std::sort(pBeginIt, pBeginIt + *pWorkingLimit);
    return violations;
}







float Anneal::calculateRgSym(std::vector<unsigned int> &beadsInUse, unsigned int const &workingLimit, Model *pModel) {
    //
    // add interdomain distances
    //
    // get coordinate of selectedBead and calculate distance to all other selected beads and convert to bin
    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*workingLimit;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;


    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < workingLimit; i++){
        Bead * tempBead = pModel->getBead(beadsInUse[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    unsigned int count = workingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, workingLimit, count, coordinates);
    }

    vector3 com = vector3(0,0,0);
    for (unsigned int i=0; i < totalCoordinates; i++){
        com += coordinates[i];
    }

    float inv= 1.0f/(float)totalCoordinates;
    com.x *= inv;
    com.y *= inv;
    com.z *= inv;

    float rgsum=0, diff;
    for (unsigned int i=0; i < totalCoordinates; i++){
        diff = (coordinates[i] - com).length();
        rgsum += diff*diff;
    }

    return sqrt(rgsum*inv);
}


float Anneal::surfaceAreaSymObject(const unsigned int workingLimit, const unsigned int totalSubUnits, float sa){
    // Surface area calculations
    const unsigned int totalToDo = workingLimit*totalSubUnits;
    std::vector<float> weights(totalToDo);
    std::fill(weights.begin(), weights.end(), probe_radius);

    return surfaceToVolume(totalToDo, weights, totalCoordinatesInObject);
}

/**
 * recalculate P(r) distribution then compare against dataset for KL divergence
 * calculation is from the coordinates
 * Steps:
 * 1. create coordinates of first subunit (workingLimit)
 * 2. create symmetry mates and add to a mster list of coordinates
 * 3. calculate distance distribution for each pair of coordinates
 */
double Anneal::calculateKLEnergySymmetry(std::vector<unsigned int> *subUnit_indices, std::vector<unsigned int> *binCount, const unsigned int indicesWorkingLimit, unsigned int &violation, Model *pModel, Data *pData) {
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    std::fill(binCount->begin(), binCount->end(), 0.0);

    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;
    Bead * tempBead;

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        tempBead = pModel->getBead((*subUnit_indices)[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

    // calculate Pr, order is irrelavant
    float distance_to;
    violation=0;
    for (unsigned int i=0; i < indicesWorkingLimit; i++) {

        tempVec1 = &coordinates[i];

        for (unsigned int next_i = i + 1; next_i < indicesWorkingLimit; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }


        for (unsigned int next_i = indicesWorkingLimit; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax

            if (distance_to < violation_limit){
                violation++;
            }
        }
    }

    for (unsigned int i=indicesWorkingLimit; i < totalCoordinates; i++) {

        tempVec1 = &coordinates[i];

        for (unsigned int next_i = i + 1; next_i < totalCoordinates; next_i++) {
            // calculate distance and convert to bin
            distance_to = ((*tempVec1) - coordinates[next_i]).length();
            ++(*binCount)[pData->convertToBin(distance_to)]; // some distances will exceed dmax
        }
    }

    return pData->getScore(*binCount);
}


float Anneal::calculateCVXVolumeSymmetry(std::vector<unsigned int> *subUnit_indices, const unsigned int indicesWorkingLimit, Model *pModel) {
    // calculate distribution
    // go through entire distance vector, count only those that are kept
    // reset binCount
    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);
    vector3 * tempVec1;
    Bead * tempBead;

    // create first subunit from selected indices and populate coordinates
    for (unsigned int i=0; i < indicesWorkingLimit; i++){
        tempBead = pModel->getBead((*subUnit_indices)[i]);
        tempVec1 = &coordinates[i];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
    }

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

    // calculate CVX Volume, order is irrelavant
    char flags[25];
    std::sprintf(flags, "qhull s FA");

    unsigned int numpoints = 3*totalCoordinates;
    coordT points[numpoints];

    const vector3 * ptr = coordinates.data();
    unsigned int next = 0;
    for (unsigned int i=0; i < totalCoordinates; i++){
        next = 3*i;
          points[next] = ptr[i].x;
        points[next+1] = ptr[i].y;
        points[next+2] = ptr[i].z;
    }

    qh_new_qhull (3, totalCoordinates, points, 0, flags, nullptr, nullptr);
    auto volume_test = (float)(qh totvol);

    //qh totarea;
    qh_freeqhull(true);
    //float calcVol = pModel->getBeadVolume()*upTo;
    return volume_test;///calcVol;
}


unsigned int Anneal::numberOfContactsFromSetSym(std::set<unsigned int> *beads_in_use, Model *pModel, unsigned int selectedIndex){

    unsigned int interContacts = numberOfContactsFromSet(beads_in_use, pModel, selectedIndex);

    auto it = pModel->getPointerToNeighborhood(selectedIndex);
    unsigned int neighborContacts = 0;

    // For the selectedIndex, determine number of contacts it makes
    auto endOfSet = beads_in_use->end();
    unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();
    const unsigned int neighhborLimit = pModel->getNeighborLimit();

    for (unsigned int i=0; i< totalNeighbors; i++){

        unsigned int neighbor = *(it+i);

        if ((neighbor < neighhborLimit ) && beads_in_use->find(neighbor) != endOfSet){
            neighborContacts += 1;
        } else if (neighbor == neighhborLimit) {
            break;
        }
    }


    // for each subunit, calculate additional contacts
    // create first subunit from selected indices and populate coordinates
    auto indicesWorkingLimit = (unsigned int) beads_in_use->size();
    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;
    Bead * tempBead;

    unsigned int index =0;
    // assemble coordinates of base subunit
    for (auto bit = beads_in_use->begin(); bit != endOfSet; ++bit){
        tempBead = pModel->getBead(*bit);
        tempVec1 = &coordinates[index];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
        ++index;
    }

    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }

    // calculate distance to all other beads
    const vector3 * pBase = &pModel->getBead(selectedIndex)->getVec();

    for (unsigned int s=indicesWorkingLimit; s<totalCoordinates; s++){
        float diff = ((*pBase) - coordinates[s]).length();
        if ( diff < contactCutOff ){
            ++neighborContacts;
        }
    }

    return neighborContacts - interContacts;
}

void Anneal::updateContactsDistribution(std::set<unsigned int> *beads_in_use_tree, Model *pModel){

    /*
     * calculate distribution of contacts
     */
    double alphaDecay = 0.63;
    unsigned int totalCounts=0;
    std::vector<double> temp(13);
    std::fill(temp.begin(), temp.end(), 0.0d);

    for (unsigned int c=0; c<13; c++){
        for(auto it = beads_in_use_tree->begin(); it != beads_in_use_tree->end(); ++it){
            if (numberOfContactsFromSet(beads_in_use_tree, pModel, *it) == c){
                totalCounts++;
            }
        }
        temp[c] = totalCounts;
    }

    double inv = 1.0d/(double)totalCounts;
    double oneminus = 1.0d - alphaDecay;

    for (unsigned int c=0; c<13; c++){
        contactsDistribution[c] = oneminus*contactsDistribution[c] + alphaDecay*temp[c]*inv;
    }
}

void Anneal::populateContactsDistributionSym(std::vector<double> & distribution, std::set<unsigned int> *beads_in_use, Model *pModel){

    for(unsigned int i=0; i<13; i++){
        distribution[i]=0.0d;
    }

    auto indicesWorkingLimit = (unsigned int)beads_in_use->size();
    const unsigned int subUnits = pModel->getNumberOfSubUnits();
    const unsigned int totalCoordinates = subUnits*indicesWorkingLimit;
    std::vector<vector3> coordinates(totalCoordinates);

    vector3 * tempVec1;
    Bead * tempBead;

    unsigned int index =0;
    for (auto it = beads_in_use->begin(); it != beads_in_use->end(); ++it){
        tempBead = pModel->getBead(*it);
        tempVec1 = &coordinates[index];
        (*tempVec1).x = tempBead->getX();
        (*tempVec1).y = tempBead->getY();
        (*tempVec1).z = tempBead->getZ();
        ++index;
    }


    unsigned int count = indicesWorkingLimit;
    for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
        pModel->transformCoordinatesBySymmetryPreCalc(s, indicesWorkingLimit, count, coordinates);
    }


    /*
     * for each bead, determine number of contacts
     */
    for (unsigned int i=0; i < indicesWorkingLimit; i++){  // create sym related subunits and add to coordinates vector
        unsigned int neighborContacts = 0;
        unsigned int stopAt = i;
        const vector3 * pBase = &coordinates[i];

        for(unsigned int k = 0; k<stopAt ; ++k){
            if ( ((*pBase) - coordinates[k]).length() < contactCutOff){
                ++neighborContacts;
            }
        }

        ++stopAt;

        for(unsigned int k=stopAt; k < totalCoordinates; ++k){
            if ( ((*pBase) - coordinates[k]).length() < contactCutOff){
                ++neighborContacts;
            }
        }

        if (neighborContacts > 12){
            ++distribution[ 12 ];// += 1;
        } else {
            ++distribution[neighborContacts];
        }

    }
}

/**
 * Model is remapped to existing lattice defined by bin_width in Data
 *
 * @param pModel
 * @param pData
 * @param name
 * @param PDBFilename
 * @return
 */

bool Anneal::initializeModelToRefineSym(Model *pModel, Data *pData, std::string name, std::string PDBFilename) {

    srand(time(0));
    contactCutOff = pModel->getNieghborCutOffLimit();
    violation_limit = pModel->getBeadRadius();
    // float * pDistance = pModel->getPointerToDistance();
    // convert distances within the large search space to bins based on input P(R)-DATA file
    // this->fillPrBinsAndAssignTotalBin( pModel, pData);
    maxbin = pModel->getMaxBin(pData);
    // adjust this part
    maxbin += 2;
    totalBins = pData->getShannonBins(); // distances greater than number of shannon are assigned to last bin for D_KL calculation
    maxbin = (maxbin > totalBins) ? maxbin : totalBins; // choose the greater of the two
    // create working observed probability distribution that encompasses search sphere
    pData->setScoringFunction(maxbin);

    // initialize Universe and fill indices
    const unsigned int totalBeadsInSphere = pModel->getTotalNumberOfBeadsInUniverse();
    std::vector<unsigned int> bead_indices(totalBeadsInSphere); // large vector ~1000's

    logger("TOTAL PTS IN UNIVERSE", std::to_string(totalBeadsInSphere));
    unsigned int * ptr = (totalBeadsInSphere != 0) ? bead_indices.data() : nullptr;
    for(unsigned int i = 0; i < totalBeadsInSphere; i++) {
        ptr[i] = i;
    }

    // returns only seed, needs to be mapped to entire Universe
    std::vector<double> prPDB(maxbin);
    pModel->createSeedFromPDBSym(PDBFilename, pData, maxbin, &prPDB);  // binCount and target prPDB is same size
    unsigned int workingLimit = pModel->getTotalInSeed();

    // sort reassigned bead model into working universe
    unsigned int locale=0;
    for(auto  it = pModel->getSeedBegin(); it != pModel->getSeedEnd(); ++it) {
        auto fit = std::find(bead_indices.begin(), bead_indices.end(), *it); // if itTrueIndex == endTrue, it means point is not within the set
        std::iter_swap(bead_indices.begin()+locale, fit);
        locale++;
    }

    if (locale != workingLimit){
        std::cout << " Locale does not equal working limit " << std::endl;
        exit(0);
    }

    std::sort(bead_indices.begin(), bead_indices.begin() + locale);

    // number of shannon bins for the model is calculated over the Universe (not the data)
    std::vector<unsigned int> binCount(maxbin);        // smallish vector, typically < 50

    logger("EXP N_S BINS", std::to_string(totalBins));
    logger("MODEL N_S BINS", std::to_string(maxbin));
    logger("BINWIDTH", formatNumber(pData->getBinWidth(), 2));
    logger("BEAD RADIUS", formatNumber(pModel->getBeadRadius(), 2));

    unsigned int violations;
    calculateModelPrDistributionSym(&bead_indices, &binCount, workingLimit, violations, pModel, pData );
    // fill binCount for first time
    float currentKL = pData->getScore(binCount);
    //
    std::set<unsigned int> beads_in_use_tree(bead_indices.begin(), bead_indices.begin() + workingLimit);

    logger("WORKINGLIMIT", std::to_string(workingLimit));
    logger("VIOLATIONS", std::to_string(violations));
    // parse violations?
    // setup parameters for hull
    // area of the model P(r) distribution
    // populate deadLimit
    // layer of beads within interconnectivity cutOff
    // randomize points within deadlimit
    // changing working limit to add/remove
    // for each randomized arrangement, does volume get smaller?, connectivity improve?  D_KL?
    EulerTour eulerTour(bead_indices.begin(), workingLimit, pModel);
    unsigned int currentNumberOfComponents  = eulerTour.getNumberOfComponents();

    if (currentNumberOfComponents > 1){
        logger("ERROR UNFIT STARTING MODEL", PDBFilename);
        return false;
    }
    logger("STARTING D_KL", formatNumber(currentKL, 7));
    char flags[] = "qhull FA";
    float current_volume = calculateCVXHULLVolume(flags, &bead_indices, workingLimit, pModel);

    // int high;
    // pModel->writeModelToFile(workingLimit, bead_indices, name, high);
    pModel->setReducedSeed(workingLimit, bead_indices);
    pModel->setStartingSet(bead_indices);
    pModel->setStartingWorkingLimit(workingLimit);
    pModel->setBeadAverageAndStdev(workingLimit, 0);

    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************        ESTIMATED LATTICE POINTS        *******************" << std::endl;
    std::cout << "*******************                                        *******************" << std::endl;
    std::cout << "*******************                 VOLUME                 *******************" << std::endl;
    logger("CONVEX HULL VOLUME", formatNumber(current_volume,1));
    std::cout << "*******************                                        *******************" << std::endl;

    logger("SHANNON BINS IN DATA", std::to_string(pData->getShannonBins()));
    logger("ZERO BIN", std::to_string( pData->getZeroBin()));
    logger("","EXITING INITIAL MODEL BUILD");

    float    tempAverageContacts=0.0;
    for (unsigned int i=0; i<workingLimit; i++){
        unsigned int temp = numberOfContactsFromSet(&beads_in_use_tree, pModel, bead_indices[i]);
        tempAverageContacts += temp;
    }

//    double runningContactsSum = calculateTotalContactSumPotentialSym(&beads_in_use_tree, &bead_indices, workingLimit, pModel);

    float average_number_of_contacts = tempAverageContacts/(float)workingLimit;
    logger("AVERAGE NUMBER CONTACTS", formatNumber(average_number_of_contacts,1));
    std::string nameOfModel = pModel->writeModelToFile2(currentKL, workingLimit, bead_indices, binCount, "subunit_remapped_x", this, pData, 0, current_volume, average_number_of_contacts);
    pModel->writeSymModelToFile(currentKL, workingLimit, bead_indices, binCount, "sym_remapped_x", this, pData, 0, current_volume, average_number_of_contacts);

    isRefine = true;
    return true;
}

void Anneal::calculateModelParametersSymmetry(std::set<unsigned int> *subUnit_indices, Model *pModel) {

    //const unsigned int subUnits = pModel->getNumberOfSubUnits();

    const vector3 vecXY(1,1,0);
    const vector3 vecZ(0,0,1);

    float maxXY=0, maxZ=0;

    for(auto & ind : *subUnit_indices){ // sorted in order because of set
        const vector3 * tempVec1 = &(pModel->getBead(ind)->getVec());
        vector3 temp = (*tempVec1)*vecXY;
        if (temp.length() > maxXY){
            maxXY=temp.length();
        }

        temp = (*tempVec1)*vecZ;
        if (temp.length() > maxZ){
            maxZ=temp.length();
        }
//        for (unsigned int s=1; s < subUnits; s++){  // create sym related subunits and add to coordinates vector
//            temp = pModel->transformVectorBySymmetry(s, *tempVec1)*vecXY;
//            if (temp.length() > maxXY){
//                maxXY=temp.length();
//            }
//        }
    }

    logger("SEARCH SPACE PARAMETERS", "FOR A CYLINDER");
    logger("HEIGHT", formatNumber(2*(maxZ + pModel->getBeadRadius()),2));
    logger("RADIUS", formatNumber(maxXY + pModel->getBeadRadius(),2));

}

vector3 Anneal::rotate(vector3 const & vec, double cosx, double sinx, double cosz, double sinz){
    return vector3(cosz*vec.x-sinz*vec.y, cosx*sinz*vec.x+cosx*cosz*vec.y - sinz*vec.z, sinx*sinz*vec.x+sinx*cosz*vec.y+cosx*vec.z);
}