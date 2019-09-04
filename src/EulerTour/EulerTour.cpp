//
// Created by xos81802 on 11/07/2018.
//

#include "EulerTour.h"
#include "../Model.h"

// Constructor
EulerTour::EulerTour(){};

EulerTour::EulerTour(const EulerTour &e2){

    for(auto & node : e2.nodes){
        nodes.emplace(node.first, Node(node.second));
    }

    for(auto & e2Tour : e2.tours){
        auto pair = tours.emplace(e2Tour.first, std::list< Node *>());
        std::list<Node *> * pListIt = &(*pair.first).second;
        for(auto & lit : e2Tour.second){
            unsigned int keytocopy = lit->getKey();
            pListIt->push_back(&(*nodes.find(keytocopy)).second);
        }
    }

    totalComponents = tours.size();
}

EulerTour::EulerTour(std::vector<unsigned int> & indices, const unsigned int workingLimit, Model *pModel){
    this->createInitialTour(workingLimit, pModel, indices);
}

EulerTour::EulerTour(std::vector<unsigned int>::const_iterator beginIt, const unsigned int workingLimit, Model *pModel){
    this->createInitialTour(workingLimit, pModel, beginIt);
}

EulerTour::EulerTour(std::set<unsigned int> & indices, Model *pModel){
    // for each node, make adjacency list
    const unsigned int neighborLimit = pModel->getNeighborLimit();

    for(std::set<unsigned int>::iterator sit = indices.begin(); sit != indices.end(); ++sit){
        // make new Node from lattice/bead
        unsigned int nodeToInsert = *sit;

        nodes.emplace(nodeToInsert, Node(nodeToInsert));
        Node * pNode = &(nodes.find(nodeToInsert)->second); // retrieve the newly made node
        //Node * pNode = &(*pNit.first).second;
        // create the neighborhood for nodes in nodes list
        // as I add each node, create neighborhood
        auto vit = pModel->getPointerToNeighborhood(pNode->getKey());

        for (unsigned int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is in node list, add to adjacency list of node (use this list to build Euler Tour)
            unsigned int neighbor = *(vit+j);
            // if not found, itIndex will report last
            if (neighbor < neighborLimit && nodes.find(neighbor) != nodes.end()){
                pNode->addNeighbor(&(nodes.find(neighbor)->second)); // if node is found to have a neighbor, add to node's adjacency list
            } else if (neighbor == neighborLimit) {
                break;
            }
        }

        addToTour(nodeToInsert);
    } // end of adding beads

    //validateNodes("from createInitialTour");
    totalComponents = (unsigned int)tours.size();
}



void EulerTour::createInitialTour(unsigned int workingLimit, Model *pModel, std::vector<unsigned int> & indices) {
    // for each node, make adjacency list
    unsigned int neighbor, nodeToInsert;
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    const unsigned int * ptr = indices.data();
    const unsigned int * neighborhood = pModel->getDirectPointerToNeighborhood();
    const unsigned int totalNeighbors = pModel->getSizeOfNeighborhood();

    for(unsigned int i=0; i < workingLimit; i++){
        // make new Node from lattice/bead
        nodeToInsert = ptr[i];

        auto pNit = nodes.emplace(nodeToInsert, Node(nodeToInsert));
        //Node * pNode = &(nodes.find(newNode)->second); // retrieve the newly made node
        Node * pNode = &(*pNit.first).second;
        //nodes.insert ( std::pair<unsigned int, Node>(nodeToInsert, Node(nodeToInsert) ) );
        //Node * pNode = &(nodes.find(nodeToInsert)->second);
        // create the neighborhood for nodes in nodes list
        // as I add each node, create neighborhood
        unsigned int locale = totalNeighbors*nodeToInsert;

        for (unsigned int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is in node list, add to adjacency list of node (use this list to build Euler Tour)
            //neighbor = *(it+j);
            neighbor = neighborhood[locale + j];
            // if not found, itIndex will report last
            auto nit = nodes.find(neighbor);
            if (neighbor < neighborLimit && nit != nodes.end()){
                pNode->addNeighbor(&(nit->second)); // if node is found to have a neighbor, add to node's adjacency list
            } else if (neighbor == neighborLimit) {
                break;
            }
        }

        addToTour(nodeToInsert);
    } // end of adding beads

//    validateList("FINAL FROM EulerTour::createInitialTour ");
    totalComponents = (unsigned int)tours.size();
}

/**
 * critical, new node must already be added to nodes!
 */
void EulerTour::createInitialTour(unsigned int workingLimit, Model *pModel, std::vector<unsigned int>::const_iterator beginIt) {

    // for each node, make adjacency list
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    nodes.clear();
    tours.clear();

    for(unsigned int i=0; i < workingLimit; i++){
        // make new Node from lattice/bead
        unsigned int nodeToInsert = *(beginIt + i);

        auto pNit = nodes.emplace(nodeToInsert, Node(nodeToInsert));
        Node * pNode = &(*pNit.first).second;

        //nodes.emplace(nodeToInsert, Node(nodeToInsert));
        //Node * pNode = &(nodes.find(nodeToInsert)->second);

        // create the neighborhood for nodes in nodes list
        // as I add each node, create neighborhood
        auto it = pModel->getPointerToNeighborhood(pNode->getKey());
        for (unsigned int j=0; j < pModel->getSizeOfNeighborhood(); j++){
            // if neighbor is in node list, add to adjacency list of node (use this list to build Euler Tour)
            unsigned int neighbor = *(it+j);
            // if not found, itIndex will report last
            auto nit = nodes.find(neighbor);
            if (neighbor < neighborLimit && nit != nodes.end()){
                pNode->addNeighbor(&(nit->second)); // if node is found to have a neighbor, add to node's adjacency list
            } else if (neighbor == neighborLimit) {
                break;
            }
        }
        //std::cout << "  " << i << " ADDING NODE -> " << nodeToInsert << std::endl;
        addToTour(nodeToInsert);
    } // end of adding beads

    totalComponents = (unsigned int)tours.size();
}


/**
 * create subTour rooted at Node
 * reroot to node found in current EulerTour
 * splice into EulerTour
 * Having a neighbor means tour is already present to be added to
 */
bool EulerTour::addToTour(unsigned int nodeToAdd){
    // get the new node from NODES map list
    Node *pNeighbor, * pNode = &(nodes.find(nodeToAdd)->second);
    // add to Tour
    if (pNode->getTotalNeighbors() == 0){ // if no existing neighbors, make a new tour
        /*
         * add node to tours as new tour
         * set node tour index mapping
         * create tour (empty)
         */
        auto pNit = tours.emplace(pNode->getKey(), std::list< Node *>());
        pNode->setRootNodeOfTour(pNode->getKey());
        (*pNit.first).second.push_back(&(*pNode));
//        tours.emplace (pNode->getKey(), std::list< Node *>());  // key is the root of the tour
//        pNode->setRootNodeOfTour(pNode->getKey());
//        tours.find(pNode->getKey())->second.push_back(&(*pNode));

    } else {
        // create subtour of new Node using Node's found neighbors
        std::list< Node * > subTour;  // subtour is created and has no membership to existing tour
        createSubTour(pNode, &subTour); // create subtour using pNode's neighbors
        unsigned int keyOfNeighborNode = pNode->getKey();
        // new Node could bridge 1 or more prior nodes
        std::set<unsigned  int > rootOfNeighborhoods;

        for (auto pit : *pNode->getIteratorToIndices()){
            pNeighbor = &nodes.find(pit)->second;
            // get tour index to merge into
            rootOfNeighborhoods.insert(pNeighbor->getRootNodeOfTour()); // how many neighborhoods does newNode belong to?
            // get pointer to list of neighbor, is it the same?
            keyOfNeighborNode = pNeighbor->getKey(); // use last neighbor as start
        }
        // what happens if node that is removed is root node for a tour?
        // set tour membership of new Node

        pNeighbor = &(nodes.find(keyOfNeighborNode)->second);

        rerootSubTour(pNeighbor, &subTour); //prepare for merging in pNeighbor tour (doesn't change tour membership of nodes, still points to old location)
        // subTour nodes do not have root at this point
        //pNode->setPointerToTour(pNeighbor->getPointerToTour()); // sets pointer to list< Node *> in vector
        //nodes.find(nodeToAdd)->second.setPointerToTour(&(tours.find(pNeighbor->getRootNodeOfTour())->second)); //update node
//        pNode->setPointerToTour(&(tours.find(pNeighbor->getRootNodeOfTour())->second)); //update node
        pNode->setRootNodeOfTour(pNeighbor->getRootNodeOfTour()); // new
        const unsigned int rootToBaseTour = pNode->getRootNodeOfTour();
        /*
         * Merge the new subTour with existing tour specified by Neighbor
         * re-root to node found incommon to a previous tour in tours
         * merge in tours[i]
         *  get base tour of neighbor to merge into
         */
        //std::list< Node * > * pExistingNeighborTour = pNeighbor->getPointerToTour();
        std::list< Node * > * pExistingNeighborTour = &tours.find(pNeighbor->getRootNodeOfTour())->second;
        //std::cout << pNeighbor->getRootNodeOfTour() << " " << tours.find(pNeighbor->getRootNodeOfTour())->first << std::endl;
        unsigned int existingRootNode = tours.find(pNeighbor->getRootNodeOfTour())->first;
        // subTour is rooted to pNode
        //               subTour => 58565
        //
        // pExistingNeighborTour => 124313426797621   (root node is 1)
        //
        // reroot        subTour => 65856
        //
        //std::cout << "EulerTour::addToTour finding node by key : " << keyOfNeighborNode << std::endl;
        auto inTour = std::find_if(pExistingNeighborTour->begin(), pExistingNeighborTour->end(), find_Node_by_key(pNeighbor->getKey()));
        //
        // locate 6 at ...267...
        //
        subTour.back() = nullptr;
        subTour.pop_back();
        //
        // pop_back      subTour => 6585
        // merge tour
        // merge subTour of newNode with existing tour
        pExistingNeighborTour->splice(inTour, subTour); // Add additional tours to base
        // subTour should be empty now
        //
        // pExistingNeighborTour => 1243134265856797621   (root node is 1)
        //
        // nodes in this tour will have mixed roots with pNeighbor
        // if new node bridges two or more tours, pExistingNeighborTour will be mixed
        if (rootOfNeighborhoods.size() > 1){ // if number of neighborhoods is greater than 1, it means the new node bridges at least two
            // find tour greater than minTourIndex, reroot tour and merge
            for(auto pit : *pNode->getIteratorToIndices()){
//            for (int j=0; j < (totalNeighbors - 1); j++){//possible underflow
                // if two neighbors are in same neighborhood, and I add the first one
                // then i have a problem if rootNodeOfTour is not updated before next neighbor is checked
//                pNeighbor = pNode->getPointerToNeighborByIndex(j);
                pNeighbor = &nodes.find(pit)->second;
                unsigned int rootOfCurrentTour = pNeighbor->getRootNodeOfTour();
//                std::cout << j << " neighbor " << pNeighbor->getKey() << std::endl;
                // if two nodes belong to same tour, only need to add tour once
                if (rootOfCurrentTour != rootToBaseTour){ // prevents neighbor that was reassigned to lower index from being reassigned again
                    std::list< Node * > * pList = &(tours.find(rootOfCurrentTour)->second);
                    // pList should not be zero size
                    rerootSubTour(pNeighbor, pList);
                    inTour = std::find_if(pExistingNeighborTour->begin(), pExistingNeighborTour->end(),find_Node_by_key(pNeighbor->getKey()));
                    // better to iterate over node list and check for tour membership
                    // for all nodes in this tour, reset tourIndex
                    for(auto pIt = pList->begin(); pIt != pList->end(); ++pIt){
                        // SLOW STEP (CAN BE OPTIMIZED)
                        //if ( (*pIt)->getPointerToTour() != pExistingNeighborTour ){
                        //  (*pIt)->setPointerToTour(pExistingNeighborTour);
                        //nodes.find((*pIt)->getKey())->second.setPointerToTour(&(tours.find(pNeighbor->getRootNodeOfTour())->second)); //update node
//                        nodes.find((*pIt)->getKey())->second.setPointerToTour(pExistingNeighborTour); //update node
                        nodes.find((*pIt)->getKey())->second.setRootNodeOfTour(existingRootNode);
                        //(*pIt)->setRootNodeOfTour(pExistingNeighborTour->front()->getKey()); // new
                        //}
                    }
                    // reassign before popping, in case pList size is 1
                    pList->back() = nullptr;
                    pList->pop_back();
                    // merge pList into existingNeighborTour
                    pExistingNeighborTour->splice(inTour, *pList);
                    // remove, erasing from vector changes address of the other elements,
                    // pList must be zero in size
                    // std::map<unsigned int, std::list< Node *> > tours;
                    auto lastTime = tours.find(rootOfCurrentTour);
                    for(auto & entries : (*lastTime).second){
                        entries = nullptr;
                    }
                    tours.erase(rootOfCurrentTour); // calls destructor on list held by tours map
                }
            }
        }

//        printList("NEWLY ADDED", &tours[nodes.find(nodeToAdd)->second.getRootNodeOfTour()]);
        // after all neighbors have been added, we then update the nodes tour
        // only want to iterate over the unique nodes in the tour rather than all members of the tour!
        //
//        for(std::list< Node * >::iterator pIt = pExistingNeighborTour->begin(); pIt != pExistingNeighborTour->end(); ++pIt){
//            // SLOW STEP (CAN BE OPTIMIZED)
//            if ( (*pIt)->getPointerToTour() != pExistingNeighborTour ){
//                (*pIt)->setPointerToTour(pExistingNeighborTour);
//                std::cout << " UPDATING " << std::endl;
//                exit(0);
//            }
//        }

    } // end of creating subtour
//    if (validateList("FINAL FROM EulerTour::addToTour ")){
//        exit(0);
//    }
    return true;
}



/**
 *               9
 *               |
 *   1 - 2 - 6 - 7
 *   |   |
 *   3 - 4
 *   124313426797621
 *  smallest tour size is 1, then 3, there is no two
 *
 *  Only reroot the tour, there is no updating of the nodes root tour
 */
void EulerTour::rerootSubTour(Node * newRoot, std::list< Node * > * subTourToLoad){

    // 124313426797621 => reroot to 6
    std::list< Node * >::iterator inTour = std::find_if(subTourToLoad->begin(), subTourToLoad->end(), find_Node_by_key(newRoot->getKey()));
    // inTour points to 6 at ..4267...
    // remove first node
    // delete (subTourToLoad->front());
    // if front of tour equals newRoot, do nothing
    auto sizeOfTour = (unsigned int)subTourToLoad->size();
    if (subTourToLoad->front()->getKey() != newRoot->getKey() && (sizeOfTour > 2)){
       // delete subTourToLoad->front();
        subTourToLoad->front() = nullptr;
        subTourToLoad->pop_front();
        // subTourToLoad => 24313426797621
        std::list< Node * > tempList;
        tempList.splice(tempList.begin(), *subTourToLoad, (*subTourToLoad).begin(), inTour);
        // tempList => 2431342
        // subTourToLoad => 6797621
        tempList.push_back(&(*newRoot));
        //      tempList => 24313426
        // transfer components from tempList to subTour
        subTourToLoad->splice((*subTourToLoad).end(), tempList);
        //printList("rerooted ", subTourToLoad);
        // subTourToLoad => 6797621-24313426
        // clear tempList
        tempList.clear();
    }

//    else if (sizeOfTour == 1) {
//        std::cout << " TOUR SIZE is ONE !!!!!  " << std::endl;
//        std::cout << " TOUR SIZE is ONE !!!!!  " << subTourToLoad->front()->getKey() << std::endl;
//    }
}


/**
 * Creates proper Euler Tour where first and last Node are the same
 */
void EulerTour::createSubTour(Node * pNode, std::list< Node * > * subTourToLoad){

    subTourToLoad->push_back(pNode);
    unsigned int totalNeighbors = pNode->getTotalNeighbors();
    unsigned int totalIterations = 2*totalNeighbors + 1;
    unsigned int index=0;

    auto pIt = pNode->getIteratorToIndices()->begin();

    for(unsigned int i=1; i<totalIterations; i++){
        if (i%2 == 0){
            subTourToLoad->push_back(pNode);
        } else {
            // add from neighborhood and not Node

//            unsigned int neighborIndex = pNode->getPointerToNeighborByIndex(index);
            unsigned int neighborIndex = *pIt;
            std::advance(pIt, 1);

            subTourToLoad->push_back(&nodes.find(neighborIndex)->second);

//            subTourToLoad->push_back(&nodes.find(pNode->getPointerToNeighborByIndex(index)->getKey())->second);
            //subTourToLoad->push_back(pNode->getPointerToNeighborByIndex(index));
            index++;
        }
    }
}

/**
 * add a node to existing set of nodes
 * first, create neighborhood from existing nodes
 * find which tour to tours to add to and merge any tours that are bridged by new node
 * returns number of nodes in tour
 *
 * Possible undefined condition, using find on nodes can lead to nodes.end()
 *
 * the control of active nodes is maintained by beads_in_use and bead_indices
 */
unsigned int EulerTour::addNode(unsigned int newNode, Model * pModel) {

    unsigned int neighbor;
    //validateNodes("FROM ADD NODE  : " +  std::to_string(newNode));

    //nodes.emplace(newNode, Node(newNode));
    //Node * pNode = &(nodes.find(newNode)->second); // retrieve the newly made node
    auto pNit = nodes.emplace(newNode, Node(newNode));
    Node * pNode = &(*pNit.first).second;
    //
    // can lead to nodes.end() which doesn't contain a node
    // somehow getting a nonsense node to return
    // build neighborhood (check from existing nodes or beads in use)
    const unsigned int neighborLimit = pModel->getNeighborLimit();
    auto it = pModel->getPointerToNeighborhood(newNode);
    for (unsigned int j=0; j < pModel->getSizeOfNeighborhood(); j++){
        // if neighbor is inside workinglimit, don't add
        neighbor = *(it+j); // index of bead used as key
        // if not found, itIndex will report last
        if (neighbor < neighborLimit && nodes.find(neighbor) != nodes.end()){ // look in nodes list
            pNode->addNeighbor(&(nodes.find(neighbor)->second)); // add neighbor to both parent and child
            // again nothing prevents me from finding nodes.end()
            // has no idea which list it is in
        } else if (neighbor == neighborLimit){
            break;
        }
    }
    //printTourSizes();
    //std::cout << " NEW NODE " << newNode << " => " << pNode->getTotalNeighbors() <<  "  key of new node " << (nodes.find(newNode)->second).getKey() << std::endl;
    // add to Tour
    addToTour(newNode);
    totalComponents = tours.size();
    //totalComponents = simpleTours.size();
//    validateNodes("AFTER ADDING NODE : " + std::to_string(newNode));
    return totalComponents;
}

/**
 *
 * @param indexOfNode
 * @return total number of components in the tour
 */
unsigned int EulerTour::removeNode(unsigned int indexOfNode){

    Node * pNodeToRemove = &(nodes.find(indexOfNode)->second);
    unsigned int indexOfNodeToRemove = pNodeToRemove->getKey();
    // check neighbors
    std::list< Node * > * pTour;// = &tours[pNodeToRemove->getRootNodeOfTour()];  // tour with node to remove
    // break the euler tour and remove the pNodeToRemove
    if (pNodeToRemove->getTotalNeighbors() > 0){

        pTour = &tours.find(pNodeToRemove->getRootNodeOfTour())->second;    // get tour that contains node to remove
        // reroot tour to NodeToRemove and remove old tour
        rerootAndReassignSubTour(pNodeToRemove, pTour);  // remove old tour and set as newly rooted tour
        pTour = &tours.find(pNodeToRemove->getRootNodeOfTour())->second;    // get tour that contains node to remove

        unsigned int keyOfNodeToRemove = pNodeToRemove->getKey();

        unsigned int rounds = 1;
        //std::cout << indexOfNode << " " << pNodeToRemove->getTotalNeighbors() << std::endl;
        while (pNodeToRemove->getTotalNeighbors() > 0) {
            //rerootSubTour(pNodeToRemove, pTour);
//            pNodeToRemove->printNeighbors();
//            std::cout << " => " << pNodeToRemove->getKey() << std::endl;
            // who is the neighbor
            auto pRightNeighbor = std::next(pTour->begin(), 1);
            auto pStartHere = std::next(pTour->begin(), 2);
//                std::cout <<  " right =>" << (*pRightNeighbor)->getKey() <<  " next => " << (*pStartHere)->getKey() << std::endl;
            Node * pNeighbor = *pRightNeighbor; // get the pointer to the node that is in the Node list and not List of the Tour
            unsigned int indexOfNeighborNode = pNeighbor->getKey();
//            std::cout <<"ADDRESS BEFORE : "<< &(*pNeighbor) << " " << pNeighbor->getKey()<< std::endl;
            // find next neighbor pair
            std::list< Node * >::iterator pNext = std::next(pTour->begin(), 2);

            if ((*pStartHere)->getKey() != keyOfNodeToRemove) { // should be XYX where X is being removed, but if not, need to find first occurrence of X after Y

                std::list< Node * >::iterator pLeftNeighbor = std::find(pStartHere, pTour->end(), pNeighbor); // look for occurrences of pNeighbor
                while(pLeftNeighbor != pTour->end()){ // find neighbor and check to see if next to it is NodeToRemove

                    pNext = std::next(pLeftNeighbor, 1); // could be .end() ?
                    if ((*pNext)->getKey() == pNodeToRemove->getKey()) {
                        break;
                    }
                    // find next
                    pStartHere = std::next(pLeftNeighbor, 1);
                    pLeftNeighbor = std::find(pStartHere, pTour->end(), pNeighbor);
                }
            }

            // potential pLeftNeighbor will reach end?
            pTour->front() = nullptr;
            pTour->pop_front();

            std::list< Node * > subTour;
            subTour.splice(subTour.begin(), *pTour, pTour->begin(), pNext); // upTo but not including pNext, pNext should point to node to remove
            //
            // if pNext is last node, then pTour will have only one element
            //
            // now I have two tours
            // 1. subTour
            // 2. pTour
            // check if I can sub back into remaining pTour
            std::set<unsigned int> nodesToCheck; // how could subTour be empty?
            for(auto & sit : subTour){
                nodesToCheck.insert(sit->getKey());
            }

            // merge if possible, check each node in nodesToCheck
            for(auto sit = nodesToCheck.begin(); sit != nodesToCheck.end(); ++sit){
                std::list< Node * >::iterator locationOfCommonNode = std::find_if(pTour->begin(), pTour->end(), find_Node_by_key(*sit));
                if (locationOfCommonNode != pTour->end()){
                    // merge subTour with pTour
                    rerootSubTour(*locationOfCommonNode, &subTour);
                    // root pTour to common node
                    subTour.back() = nullptr;
                    subTour.pop_back();
                    pTour->splice(locationOfCommonNode, subTour);
                    break;
                }
            }
            /*
             * subTour should be empty if node in common exists within pTour
             * make a new Tour from subTour if it couldn't be merged(linked)
             */
            if (subTour.size() > 0){
                unsigned int newRootNodeForTour = subTour.front()->getKey();
                auto pTourToAdd = tours.find(newRootNodeForTour);
                if (pTourToAdd == tours.end()){ // if tour doesn't exist, make it using emplace and grab the returned pair
                    pTourToAdd = tours.emplace(newRootNodeForTour, std::list<Node *>()).first;
                } else {
                    std::cout << "pExisting Tour " << std::endl;
                    exit(0);
                }
                // get the pointer to the new Node List
                //std::list< Node *> * pNewRootNodeTour = &(tours.find(newRootNodeForTour)->second); // std::map<int, std::list< Node *> > tours;
                std::list< Node *> * pNewRootNodeTour = &(*pTourToAdd).second;
                // newly create tour is empty so add subTour to it and reset RootNodes for each member

                pNewRootNodeTour->splice(pNewRootNodeTour->begin(), subTour);

                for(auto & it : *pNewRootNodeTour){
                    //std::cout << newRootNodeForTour << " BEFORE => " << (*(nodes.find(it->getKey()))).second.getRootNodeOfTour() << std::endl;
                    it->setRootNodeOfTour(newRootNodeForTour);
                    //(*(nodes.find(it->getKey()))).second.setRootNodeOfTour(newRootNodeForTour);
                    //std::cout << newRootNodeForTour << "  AFTER => " << (*(nodes.find(it->getKey()))).second.getRootNodeOfTour() << std::endl;
                }
            }
            // pTour should be rooted at NodeToRemove at this point
            // should stay that way until pNodeToRemove is gone
            //(nodes.find(indexOfNodeToRemove)->second).removeNeighbor(&nodes.find(indexOfNeighborNode)->second);
            pNodeToRemove->removeNeighbor(pNeighbor); // should be in balance after removing
//            if ((*nodes.find(indexOfNodeToRemove)).second.isNeighborPresent(indexOfNeighborNode)){
//                std::cout << "IT didn't remove " << std::endl;
//                exit(0);
//            }
            // if no more pNodeToRemove in pTour, can't reroot to node to remove
            rounds++;
        }

        // remove pTour (should contain only a single node and be rooted to nodeToRemove)
        auto findIt = tours.find(pNodeToRemove->getRootNodeOfTour());
        if (findIt != tours.end()){
            for(auto & lit : (*findIt).second){
                lit = nullptr;
            }
            (*findIt).second.clear();
            tours.erase(pNodeToRemove->getRootNodeOfTour()); // remove from tour
        }

    } else { // no neighbors
        //
        // tours is a map => std::map<int, std::list< Node *> > tours; // key is the root of the tour
        // list must be cleared first
        auto pTour = tours.find(pNodeToRemove->getRootNodeOfTour());
        if (pTour != tours.end()){

            if ((*pTour).second.size() > 1){
                std::cout << "FAULTY tour size " << (*pTour).second.size() << std::endl;
                exit(0);
            }

            auto pList = &tours.find(pNodeToRemove->getRootNodeOfTour())->second; // get iterator to list list

            for(auto & lit : *pList){
                lit = nullptr;
            }
            pList->erase(pList->begin(), pList->end()); // erase invalidates existing iterators
            pList->clear();
            tours.erase(pNodeToRemove->getRootNodeOfTour());
        }
    }

    // remove node
    // remove the node from nodes
    if (pNodeToRemove->getTotalNeighbors() == 0){
//        std::cout << " REMOVING NODE : " << pNodeToRemove->getKey() << " ROOTNODE :" << pNodeToRemove->getRootNodeOfTour() << std::endl;
//        if (nodes.find(pNodeToRemove->getKey()) == nodes.end()){
//            std::cout << "ERROR NOT FOUND in NODES LIST!  CANT REMOVE " << pNodeToRemove->getKey() << std::endl;
//            exit(0);
//        }
        // nodes is a map => std::map<int, Node> nodes;
        nodes.erase(pNodeToRemove->getKey()); // calls destructor on Node
    } else {
        std::cout << " NODE NEIGHBOR LIST NOT EMPTY PROBLEM WITH CODE?" << std::endl;
        exit(0);
    }

    //printTourSizes();
//    validateNodes("FROM REMOVE NODE AT END : " +  std::to_string(indexOfNode));
//    if (validateList("FINAL FROM EulerTour::removeNode ")){
//        exit(0);
//    }
    totalComponents = tours.size();
    return totalComponents;
}

/**
 * passing a pointer to a list of Nodes
 *
 */
void EulerTour::rerootAndReassignSubTour(Node * newRoot, std::list< Node * > * subTourToLoad){

    //std::cout << newRoot->getKey() << " subTour size : " << subTourToLoad->size() << std::endl;
    if (subTourToLoad->front()->getKey() != newRoot->getKey()){
        // printList("================> REROOT_AND_REASSIGN_SUBTOUR REROOT Before", subTourToLoad);
        //
        // subTourToLoad => 124313426797621   (root node is 1)
        //       newRoot => 6
        auto inTour = std::find_if(subTourToLoad->begin(), subTourToLoad->end(), find_Node_by_key(newRoot->getKey()));
        Node * pN1 = subTourToLoad->front();
        subTourToLoad->front() = nullptr;
        subTourToLoad->pop_front();
        //
        //        inTour => 6 at ...4267...
        // subTourToLoad => 24313426797621
        //
        std::list< Node * > tempList;
        tempList.splice(tempList.begin(), *subTourToLoad, (*subTourToLoad).begin(), inTour); // exclude inTour position, transfer is upto
        //      tempList => 2431342
        // subTourToLoad => 6797621
        tempList.push_back(newRoot);
        //      tempList => 24313426
        // copy to subTour
        subTourToLoad->splice((*subTourToLoad).end(), tempList);
//        printList("          => rerootAndReassignSubTour", subTourToLoad);
        // subTourToLoad => 679762124313426
        resetRootNodesInSubTourOfTour(subTourToLoad);
        // clear tempList
        tempList.clear();
    }
}


/**
 * subTour is a pointer to tour
 * removes old key from the tours map.
 * what happens if they oldNodeKey and newRootNodeForTour are the same?
 */
void EulerTour::resetRootNodesInSubTourOfTour(std::list<Node *> * subTour){

    unsigned int oldNodeKey = subTour->front()->getRootNodeOfTour();
    unsigned int newRootNodeForTour = subTour->front()->getKey();

    if (oldNodeKey == newRootNodeForTour){
        for(auto it = subTour->begin(); it != subTour->end(); ++it) {
            auto pTempNode = &nodes.find((*it)->getKey())->second;
            pTempNode->setRootNodeOfTour(newRootNodeForTour);
        }
    } else {

        auto pTempTour = tours.find(newRootNodeForTour);
        if (pTempTour == tours.end()) {
            pTempTour = tours.emplace(newRootNodeForTour, std::list<Node *>()).first;
            //tours.insert(std::pair<unsigned int, std::list<Node *> >(newRootNodeForTour, std::list<Node *>()));  // key is the root of the tour
        } else {
            std::cout  << "resetRootNodesInSubTourOfTour SUBTOUR " << std::endl;
            exit(0);
        }

        // get the pointer to the new Node List
        //std::list< Node *> * pNewRootNodeTour = &(tours.find(newRootNodeForTour)->second); // std::map<int, std::list< Node *> > tours;
        std::list< Node *> * pNewRootNodeTour = &(*pTempTour).second; // std::map<int, std::list< Node *> > tours;
        /*
         * Add nodes to newly created tour while setting Root node
         */
        for(auto it = subTour->begin(); it != subTour->end(); ++it) {
            auto pTempNode = &nodes.find((*it)->getKey())->second;
            pTempNode->setRootNodeOfTour(newRootNodeForTour);
            pNewRootNodeTour->push_back(&(*pTempNode));
        }

        // if oldNodeKey and newRootNodeForTour are not the same, then insert will not insert, must delete first and create new entry
        auto pOldKeyTour = tours.find(oldNodeKey);
        if (oldNodeKey != newRootNodeForTour && pOldKeyTour != tours.end()){
            pOldKeyTour->second.erase(pOldKeyTour->second.begin(), pOldKeyTour->second.end());
            pOldKeyTour->second.clear(); // clear the list
            tours.erase(oldNodeKey); // remove tour from MAP based on key.
        }
    }

}


// reuse Tour and reset
unsigned int EulerTour::newTour(std::vector<unsigned int>::iterator beginIt, unsigned int workingLimit, Model *pModel) {

    // clear pointers to nodes
    for (auto it=tours.begin(); it!=tours.end(); ++it){
        // delete points in tour
//        for(std::list<Node *>::iterator lit = it->second.begin(); lit != it->second.end(); ++lit){
//            *lit = nullptr;
//            //it->second.erase(lit);
//        }

        for(auto & lit : (*it).second){
            lit = nullptr;
        }
        while(!(*it).second.empty()){
            delete (*it).second.front();
            (*it).second.pop_front();
        }
        it->second.clear();
    }

    // remove all nodes
    auto it = nodes.begin();
    while(it != nodes.end()){
        it = nodes.erase(it);
    }

    tours.clear();
    nodes.clear();

    //pSelectedLattice = &beginIt;
    this->createInitialTour(workingLimit, pModel, beginIt);
    std::cout << "Tour size after initialization of new tour: " << tours.size() << std::endl;
    return totalComponents;
}

bool EulerTour::validateTour(std::list<Node *> * tourtocheck){
    // each pair of numbers needs a reverse sequence
    // go through each element of list and validate it is a node
    //std::map<int, int> pairs;
    std::list<Node *>::iterator stopIt = std::next(tourtocheck->end(),-1);
    for(std::list<Node *>::iterator nodeIt=tourtocheck->begin(); nodeIt != stopIt; nodeIt++){

        std::list<Node *>::iterator next = std::next(nodeIt, 1);

        int firstKey = (*nodeIt)->getKey();
        int secondKey = (*next)->getKey();

        bool isTrue = true;

        for(std::list<Node *>::iterator nodeIt2=tourtocheck->begin(); nodeIt2 != stopIt; nodeIt2++){

            std::list<Node *>::iterator next2 = std::next(nodeIt2, 1);
            int firstKey2  = (*nodeIt2)->getKey();
            int secondKey2 = (*next2)->getKey();
            if (firstKey == secondKey2 && secondKey == firstKey2){
                isTrue = false;
                break;
            }
        }

        if (isTrue){
            std::cout << "InVALID TOur " << firstKey << " " << secondKey << std::endl;
            return false;
        }
    }

    return true;
}

bool EulerTour::validateList(std::string text){

    //std::cout << " VALIDATING LISTS IN TOURS :::::: " << text << std::endl;
    int count = 1;
    for(auto it = tours.begin(); it != tours.end(); it++){

        unsigned int tourIndex = it->first; // tour index should be an existing node
        std::list<Node *> * pNodelist = &it->second;

        // check if the root of tour is in the nodes list
        if (nodes.find(tourIndex) == nodes.end()){
            std::cout << " @ " << text << " Node Does Not exist => " << tourIndex << " tour size : " << pNodelist->size()<< std::endl;
            return true;
        }

        // check if nodelist size violates limit, can only be 1 or greater than 2 such as a-b-a
        if (pNodelist->size() == 2){
            std::cout << " @ " << text << " LIST SIZE == 2 ************************************************* tourindex => " << tourIndex << std::endl;
            return true;
        }

        // check if nodelist contains a legitimate tour (every edge contains reciprocal pair)
        if (pNodelist->size() > 2 && !validateTour(pNodelist)){
            std::cout << " @ " << text << " INVALID LIST  " << tourIndex << std::endl;
            printList("INVALID LIST or TOUR ", pNodelist);
            return true;
        }

        // checks that first and last nodes of tour are the same
        if (*(pNodelist->begin()) != pNodelist->back()){
            printList(" DOES NOT EQUAL => FIRST AND LAST ELEMENT OF TOUR ", pNodelist);
            std::cout << " DOES NOT EQUAL => FIRST AND LAST ELEMENT OF TOUR " << std::endl;
            return false;
        }


        // go through each element of list and validate it is a node
        for(std::list<Node *>::iterator nodeIt=pNodelist->begin(); nodeIt != pNodelist->end(); nodeIt++){
            int key = (*nodeIt)->getKey();
            if (nodes.find(key) == nodes.end()){ // if not found we have an extra node
                std::cout << "  =>>> LIST " << count << " LOCALE : " << text << " ======>        tour MAP IND : " << tourIndex << " of " << tours.size() << std::endl;
                std::cout << "  =>>> LIST " << count << " LOCALE : " << text << " ======>       NODELIST SIZE : " << pNodelist->size() << std::endl;
                std::cout << "  =>>> LIST " << count << " LOCALE : " << text << " ======>                KEY  : " << key << std::endl;
                std::cout << "  =>>> NONSENSE KEY  : " << key << std::endl;
                printList("MESSED UP TOUR ", pNodelist);
                return true;
            }
        }
        count++;
    }

    return false;
}


bool EulerTour::validateNodes(std::string str){

    for (auto & it : nodes){
        // check for repeats
        Node * tempNode = &it.second;

        if (!(it.second).validate()){
            std::cout << str << std::endl;
            std::cout << " INVALID NODE " << it.first << " " << tempNode->getKey() << std::endl;
            exit(0);
        }

        // for each node, check that his neighbors are present in node list
        unsigned int totalNeighbors = tempNode->getTotalNeighbors();
        std::set<unsigned int> * pSet = tempNode->getIteratorToIndices();
        int count = 0;

        for(std::set<unsigned int>::iterator nit = pSet->begin(); nit != pSet->end(); ++nit){ // go through neighbors of selected node and check for reciprocity
            if (nodes.find(*nit) == nodes.end()){ // listed neighbor in Node does not exist in workingset of Euler Tour
                std::cout << str << std::endl;
                std::cout << count << " ( " << totalNeighbors << " ) " << " NEIGHBOR NODE NOT FOUND IN NODES LIST " << tempNode->getKey() << " neigh => " << *nit<< std::endl;
                std::cout << count << " BASE NODE NEIGHBORS : " << std::endl;
                tempNode->printNeighbors();
                exit(0);
            }

            // check if neighbor has current node as neighbor
            if ( !(nodes.find(*nit)->second.isNeighborPresent(tempNode->getKey())) ){
                std::cout << str << std::endl;
                std::cout << count << " NODE - NEIGHBOR relationship invalid " << std::endl;
                std::cout << " BASE NODE : " << tempNode->getKey() << " has neighbor <=> " << nodes.find(*nit)->second.getKey() << " lack of reciprocity" << std::endl;
                // print neighborhood of node
                std::cout << count << " BASE NODE NEIGHBORS : " << std::endl;
                tempNode->printNeighbors();
                std::cout << count << "  NEIGHBOR NEIGHBORS : " << std::endl;
                nodes.find(*nit)->second.printNeighbors();
                exit(0);
            }
            count++;
        }
    }
    return true;
}


/**
 * Print a tour, diagnostic tool
 */
void EulerTour::printList(std::string text, std::list< Node * > * list){
    int tourLength = list->size();

    int count=1;
    std::cout << text << std::endl;
    std::cout << "      LIST START Length : " << tourLength << std::endl;
    std::cout << "              ROOT NODE : " << list->front()->getKey() << std::endl;
    for(std::list< Node * >::iterator iterator = list->begin(); iterator != list->end(); iterator++) {
        std::cout << "                   LIST : (" << count << ") " << (*iterator)->getKey() << " root -> " << (*iterator)->getRootNodeOfTour() << std::endl;
        count++;
    }
    std::cout << " _________________________ " << std::endl;
}