//
// Created by xos81802 on 21/08/2018.
//

#include <iostream>
#include "Example.h"

Example::Example(std::string what) : type(what){

    if (type == "anchor"){
        this->printAnchor();
    }

    if(type=="helical"){
        this->printHelical();
    }

    if(type=="help"){
        this->printManual();
    }
}

void Example::printAnchor(){

    std::cout << "Created anchor_example.txt  " << std::endl;
    std::cout << "Please modify file specific to your modeling problem  " << std::endl;

    std::string nameOf = "anchor_example.txt";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    fprintf(pFile, "# This example describes a single unknown component designated B \n");
    fprintf(pFile, "# Component B has a volume 8000 Angstroms^3 and must be specified \n");
    fprintf(pFile, "# COMPONET B has two known contacts with input SEED PDB \n");
    fprintf(pFile, "# Input SEED PDB has a chain A and B contacts A at resid locations 223 and 249 \n");
    fprintf(pFile, "# All anchor points are specified by RESID and CHAIN \n");
    fprintf(pFile, "COMPONENT_ID B VOLUME %d\n",8000);
    fprintf(pFile, "RESID %d CHAIN A COMPONENT_ID B\n",223);
    fprintf(pFile, "RESID %d CHAIN A COMPONENT_ID B\n",249);

    fclose(pFile);
}


void Example::printHelical(){

    std::cout << "Created helical_params.txt  " << std::endl;
    std::cout << "Please modify file specific to your modeling problem  " << std::endl;

    std::string nameOf = "helical_params.txt";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    fprintf(pFile, "# This example describes a fiber with helical symmetry \n");
    fprintf(pFile, "#  PITCH : the distance between Valleys or Hills (Angstroms) \n");
    fprintf(pFile, "# RADIUS : the estimated width of the fiber divided by 2 (Angstroms) \n");
    fprintf(pFile, "#   RISE : the distance between successive subunits (Angstroms) => pitch/number_of_subunits_per_turn \n");
    fprintf(pFile, "#  THETA : angle that describe rotation of next subunit => 360/number_of_subunits_per_turn \n");
    fprintf(pFile, "# VOLUME : estimated volume of subunit (Angstroms^3) \n");
    fprintf(pFile, "#    PGS : Point Group Symmetry, C1 is default, if 2 intertwinned fibers => C2, etc  \n");
    fprintf(pFile, "PITCH %d\n",80);
    fprintf(pFile, "RADIUS %d \n",45);
    fprintf(pFile, "RISE %.1f \n",10.8);
    fprintf(pFile, "THETA %.1f \n",30);
    fprintf(pFile, "VOLUME %d \n",24718);
    fprintf(pFile, "PGS C1 \n");


    fclose(pFile);
}


void Example::printManual(){

    std::cout << "Created helical_params.txt  " << std::endl;
    std::cout << "Please modify file specific to your modeling problem  " << std::endl;

    std::string nameOf = "manual.txt";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    fprintf(pFile, "# This example describes a fiber with helical symmetry \n");
    fprintf(pFile, "#  PITCH : the distance between Valleys or Hills (Angstroms) \n");
    fprintf(pFile, "# RADIUS : the estimated width of the fiber divided by 2 (Angstroms) \n");
    fprintf(pFile, "#   RISE : the distance between successive subunits (Angstroms) => pitch/number_of_subunits_per_turn \n");
    fprintf(pFile, "#  THETA : angle that describe rotation of next subunit => 360/number_of_subunits_per_turn \n");
    fprintf(pFile, "# VOLUME : estimated volume of subunit (Angstroms^3) \n");
    fprintf(pFile, "#    PGS : Point Group Symmetry, C1 is default, if 2 intertwinned fibers => C2, etc  \n");
    fprintf(pFile, "PITCH %d\n",80);
    fprintf(pFile, "RADIUS %d \n",45);
    fprintf(pFile, "RISE %.1f \n",10.8);
    fprintf(pFile, "THETA %.1f \n",30);
    fprintf(pFile, "VOLUME %d \n",24718);
    fprintf(pFile, "PGS C1 \n");


    fclose(pFile);
}