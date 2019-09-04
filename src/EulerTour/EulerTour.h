//
// Created by xos81802 on 11/07/2018.
//

#ifndef SUGATAMA_EULERTOUR_H
#define SUGATAMA_EULERTOUR_H

#include <string>
#include <vector>
#include <map>
#include <list>
#include <cstdlib>
#include <memory>
#include <unordered_set>
#include "Node.h"
#include <ctime>
#include <set>

class Model;

/**
 *   Class is used to monitor the graph connectivity
 *   EulerTour eulerTour(beginIt, subUnitWorkingLimit, pModel);
 *   currentNumberOfComponents = eulerTour.getNumberOfComponents();
 */

class EulerTour {

    struct find_Node_by_key : std::unary_function< Node *, bool>{
        find_Node_by_key(unsigned int keyToFind) : key(keyToFind){}
        bool operator () (Node * p) { return p->getKey() == key; }
    private:
        unsigned int key;
    };

    unsigned int totalComponents;
    std::map<unsigned int, Node> nodes; // use shared pointer? make instance on heap
    std::map<unsigned int, std::list< Node *> > tours; // key is the root of the tour

    void createInitialTour(unsigned int workingLimit, Model *pModel, std::vector<unsigned int>::const_iterator beginIt);
    bool addToTour(unsigned int nodeToAdd);

    void createSubTour(Node * pNode, std::list< Node * > * subTourToLoad);
    void rerootSubTour(Node * newRoot, std::list< Node * > * subTourToLoad);
    void rerootAndReassignSubTour(Node * newRoot, std::list< Node * > * subTourToLoad);

    void printList(std::string text, std::list< Node * > * list);
    //void resetRootNodesInSubTour(std::list<Node *> * subTour);
    void resetRootNodesInSubTourOfTour(std::list<Node *> * subTour);
    bool validateTour(std::list<Node *> * tourtocheck);

public:
    // pass in iterator to beginning of selected lattice points that are sorted upto WorkingLimit
    // want to use this Class to determine connectivity of a graph
    // Connectivity is the numberOfComponents
    //
    EulerTour();
    EulerTour(std::vector<unsigned int>::const_iterator beginIt, const unsigned int workingLimit, Model *pModel);
    EulerTour(std::vector<unsigned int> & indices, const unsigned int workingLimit, Model *pModel);
    EulerTour(std::set<unsigned int> & beginIt, Model *pModel);

    EulerTour(const EulerTour &e2);
    ~EulerTour(){
        for(auto & tour : tours){
            for(auto & ele : tour.second){
                ele = nullptr;
            }
            tour.second.clear();
        }
        tours.clear();
        nodes.clear();
    }

    void swap(EulerTour & other) noexcept {
        using std::swap;
        std::swap(nodes, other.nodes);
        other.totalComponents = totalComponents;
        std::swap(tours, other.tours);
    }


    EulerTour & operator=(const EulerTour & tourtomove) {

        EulerTour tmp(tourtomove);
        tmp.swap(*this);

        return *this;
    }


    unsigned int addNode(unsigned int latticePoint, Model *pModel);
    unsigned int removeNode(unsigned int indexOfNode);

    unsigned int getNumberOfComponents(){return totalComponents;}

    unsigned int newTour(std::vector<unsigned int>::iterator beginIt, unsigned int workingLimit, Model *pModel);

    void createInitialTour(unsigned int workingLimit, Model *pModel, std::vector<unsigned int> &indices);

    bool validateList(std::string text);
    bool validateNodes(std::string str);

    void printTourInfo(){
        for(auto & tour : tours){
            std::cout << "TOUR : " << tour.first << " SIZE : " << tour.second.size() << std::endl;
        }
    }

    const std::list<Node *> * const getPointerToLargestTour(){
        unsigned int max = 0;
        std::list< Node *> * pList = nullptr;
        for(auto & tour : tours){
            if (tour.second.size() > max){
                max = tour.second.size();
                pList = &tour.second;
            }
        }
        return pList;
    }

    std::map<unsigned int, std::list< Node *> > * getTours(){ return &tours;} // key is the root of the tour
};

#endif //SUGATAMA_EULERTOUR_H
