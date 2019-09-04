//
// Created by xos81802 on 11/07/2018.
//

#ifndef SUGATAMA_NODE_H
#define SUGATAMA_NODE_H


#include <string>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <set>
#include <cstdlib>
#include <iostream>
/**
 * Key is unique to the node
 */


class Node {

private:
    unsigned int key;
    std::set< unsigned int > adjacencyListIndex;

    unsigned int totalNeighbors = 0;
    unsigned int rootNode;

    bool accessed = false;
    void printAdjacencyList(std::string text);

public:
    //Node();
    explicit Node(unsigned int key);
    Node(const Node &n2){
        this->key = n2.key;
        this->rootNode = n2.rootNode;
        for(auto & it : n2.adjacencyListIndex){
            this->adjacencyListIndex.insert(it);
        }
        this->totalNeighbors = adjacencyListIndex.size();
        this->accessed = n2.accessed;
    }

    void swap(Node & other) noexcept {
        using std::swap;
        other.key = key;
        other.accessed = accessed;
        std::swap(adjacencyListIndex, other.adjacencyListIndex);
        other.totalNeighbors = totalNeighbors;
        other.rootNode = rootNode;
    }

    Node & operator=(const Node & nodeToCopy) {

        Node tmp(nodeToCopy);
        tmp.swap(*this);

        return *this;
    }

    Node(Node && other) noexcept  {
        other.swap(*this);
    }

    ~Node(){
        adjacencyListIndex.clear();
    }

    unsigned int getKey() const {return key;}
    void addNeighbor(Node * pNode);
    void removeNeighbor(Node * pNode);
    void printNeighbors();

    void setAccessed(bool value){this->accessed = value;}

    //void setPointerToTour(std::list < Node * > * pointer);//{ pointerToTour = pointer;} // this points to a map which should not invalidated
//    std::list < Node * > * getPointerToTour(){ return pointerToTour;}

    void setRootNodeOfTour(unsigned int index){rootNode = index;}
    unsigned int getRootNodeOfTour(){return rootNode;}

    bool getAccessed(){return this->accessed;}

    unsigned int getTotalNeighbors(){return adjacencyListIndex.size();}
//    Node * getPointerToNeighborByIndex(int index){return adjacencyList[index];}
    bool isNeighborPresent(unsigned int index);

    std::set<unsigned int > * getIteratorToIndices(){return &adjacencyListIndex;}

    bool validate();
};


#endif //SUGATAMA_NODE_H
