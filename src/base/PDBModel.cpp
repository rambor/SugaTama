//
// Created by xos81802 on 13/07/2018.
//

#include "PDBModel.h"

PDBModel::PDBModel(std::string file, bool rna, bool discardWaters, float lower) {

    filename = file;
    ifRNA = rna;
    ifWaters = false;
    watersPerResidue = 30;

    std::string line, tempResi;
    // open file and read in contents
    std::ifstream scxFile (filename.c_str());

    if(scxFile.fail()){
        //File does not exist code here
        std::string alert="";
        char buffer[80];
        std::snprintf(buffer, 80, " ******* ERROR => File does not exist :  %s\n", filename);
        alert.append(buffer);
        throw std::invalid_argument(alert);
    }

    unsigned int fileLength = 0;
    volume=0.0;
    cutoff = lower*lower;

    atomType.reserve(5000);
    chainID.reserve(5000);
    resID.reserve(5000);

    waterLines.reserve(500*watersPerResidue);

    x.reserve(5000);
    y.reserve(5000);
    z.reserve(5000);

    boost::regex pdbStart("ATOM");
    boost::regex hetatm("HETATM");
    boost::regex wat("HOH");
    boost::regex pdbX("-*[0-9]+.[0-9]+");
    boost::regex numberFormat("[0-9]+.[0-9]+");
    boost::regex ifHydrogen("H[A-Z]+");
    boost::regex edgeRadiusFormat("EDGE RADIUS");

    std::cout << "  =>   READING PDB FILE :  " << filename << std::endl;
    if (scxFile.is_open()) {
        while(!scxFile.eof()) {
            getline(scxFile, line); //this function grabs a line and moves to next line in file
            // string::substring(position,length)
            // Check if line starts with ATOM or HETATM and exclude hydrogens
            if (line.length() > 0 && boost::regex_search(line, edgeRadiusFormat)){
                std::vector<std::string> contents;
                boost::split(contents, line, boost::is_any_of("\t  "), boost::token_compress_on);
                auto it = std::find(contents.begin(), contents.end(), "RADIUS");
                auto index = (unsigned int)std::distance(contents.begin(), it) + 1;
                for(unsigned int i=index; i<contents.size(); i++){
                    if (boost::regex_search(contents[i], numberFormat)){
                        edge_radius = std::strtof(contents[i].c_str(), nullptr);
                        break;
                    }
                }

                if (edge_radius > 1.3){
                    found_edge_radius = true;
                } else {
                    std::string alert="";
                    char buffer[80];
                    std::snprintf(buffer, 80, " ******* ERROR => CHECK EDGE RADIUS :  %f\n", edge_radius);
                    alert.append(buffer);
                    throw std::invalid_argument(alert);
                }
            }

            if ((line.length() > 0 && (boost::regex_search(line.substr(0, 6), pdbStart) || boost::regex_search(line.substr(0, 6), hetatm)) && !boost::regex_search(line.substr(17,3), wat)) && boost::regex_search(line.substr(31,8),pdbX) && !boost::regex_search(line.substr(12,4), ifHydrogen)) {

                x.push_back(strtof((char*)line.substr(30,8).c_str(), nullptr));
                y.push_back(strtof((char*)line.substr(38,8).c_str(), nullptr));
                z.push_back(strtof((char*)line.substr(46,8).c_str(), nullptr));

                //Atom type taken from rows 13-16 needs to be converted to the actual atom, e.g., CB1 is C
                atomType.push_back(line.substr(12,4));//

                // residue name, three letter for protein, two for nucleic acids
                tempResi = line.substr(17,3);

                // reassign residue abbreviations for RNA
                if (ifRNA){  // A => ALA, G => GLY, C => CYS
                    if ((tempResi == "  A") || (tempResi == "ADE") || (tempResi == " rA") || (tempResi == "A  ")){
                        tempResi = " rA";
                    } else if ((tempResi == "  G") || (tempResi == "GUA") || (tempResi == " rG") || (tempResi == "G  ")) {
                        tempResi = " rG";
                    } else if ((tempResi == "  U") || (tempResi == "URI") || (tempResi == " rU") || (tempResi == "U  ")) {
                        tempResi = " rU";
                    } else if ((tempResi == "  C") || (tempResi == "CYT") || (tempResi == " rC") || (tempResi == "C  ")) {
                        tempResi = " rC";
                    }
                }

                resi.push_back(tempResi);             // residue name Protein (ALA, GLY, ...), RNA (rA, rG, rU, rC), DNA (dA, dG, dU, rC)
                resID.push_back( atoi( line.substr(22,4).c_str()) ); // residue sequence number
                chainID.push_back(line.substr(21,1));

                atomVolume.push_back((float)SUGATAMA_FUNCTIONS_H::residueToVolume( atomType.back(), tempResi ));

                atomicRadii.push_back((float)std::cbrt(atomVolume[fileLength]*0.75/M_PI));

                volume += atomVolume.back(); //FUNCTIONS_RPR::residueToVolume( line.substr(17,3), atomType.back() );
                fileLength++;

            } else if (!discardWaters && line.length() > 20 && (line.substr(17,3) == "HOH")) {
                waterLines.push_back(line);
            }
            // keep HETATM flag - say you have a heme?
            // WATERS r in lines containing HETATM
            // if include waters is set, must
        }

        totalAtoms = fileLength;
        std::cout << "\tTotal atoms: " << totalAtoms << " ALGEBRAIC VOLUME " << volume << std::endl;
    }

    scxFile.close();
    scxFile.clear();

    trimmedAtomType = atomType;
    trimmedResi = resi;

    // center coordinates and convert to r, theta, phi
    rThetaPhiAtomType = new float[totalAtoms*5];

    atomicNumbers = new unsigned int[totalAtoms];
    atomicExpVolume = new float[totalAtoms];

    centeredX = new float[totalAtoms];
    centeredY = new float[totalAtoms];
    centeredZ = new float[totalAtoms];

    // send in x,y,z coordintes return centered coordinates
    SUGATAMA_FUNCTIONS_H::dmaxFromPDB(x,y,z, &dmax, centeredX, centeredY, centeredZ, totalAtoms);


    unsigned int locale;
    unsigned int atomicNumber;

    float sumx, sumy, sumz=0;

    for (unsigned int n=0; n < totalAtoms; n++) {
        sumx+= x[n];
        sumy+= y[n];
        sumz+= z[n];
        /*
         Convert to spherical coordinates
         xyz_to_rtp(x[],y[],z[])
         pass in &address to the array elements and dereference in the function to get value
         */
        // float * xyz_to_rtp(const float &tempx, const float &tempy, const float &tempz);

        float * pconvertXYZ = SUGATAMA_FUNCTIONS_H::xyz_to_rtp((const float &) centeredX[n], (const float &) centeredY[n], (const float &) centeredZ[n]);

        atomicNumber = convertToAtomicNumber(atomType[n], resi[n]);

        locale = n*5;
        rThetaPhiAtomType[locale] = *pconvertXYZ;               // [0] r
        rThetaPhiAtomType[locale+1] = cosf(*(pconvertXYZ+1));   // [1] cos(theta)
        rThetaPhiAtomType[locale+2] = *(pconvertXYZ+2);         // [2] phi
        rThetaPhiAtomType[locale+3] = atomicNumber;             // [3] atomic number
        rThetaPhiAtomType[locale+4] = 1.0;                      // [4]

        atomList.insert( atomicNumber );                        // Set of unique atomic numbers for calculating f(q)
        atomicNumbers[n] = atomicNumber;

        atomicExpVolume[n] = cbrt( atomVolume[n]*atomVolume[n] )*invFourPI; // volume^(2/3)*1/(4*PI)

        // call trim on resi names and atom types
        boost::algorithm::trim (trimmedAtomType[n]);
        boost::algorithm::trim (trimmedResi[n]);
    }

    centeringX=sumx/(float)totalAtoms;
    centeringY=sumy/(float)totalAtoms;
    centeringZ=sumz/(float)totalAtoms;
    centeringVector = vector3(centeringX, centeringY, centeringZ);

    // atomList is by atomic number, but water will be a molecule, give it a special number of 99
    atomList.insert(99); // Add unique for water.
    totalResidues = (unsigned int)resID.size();
    // maximum theoretical amount
    keptWaters = new Coords[watersPerResidue*totalResidues];
    waterCount = 0;

    if (!discardWaters && (waterLines.size() > 0)){
        waterCount = (unsigned int)waterLines.size();
        float invwaterCount = 1.0f/waterCount;
        float waterX=0, waterY=0, waterZ=0, aveX, aveY, aveZ;

        for(int w=0; w < waterCount; w++){
            keptWaters[w].x = (float)strtod((char*)waterLines[w].substr(30,8).c_str(), nullptr);
            keptWaters[w].y = (float)strtod((char*)waterLines[w].substr(38,8).c_str(), nullptr);
            keptWaters[w].z = (float)strtod((char*)waterLines[w].substr(46,8).c_str(), nullptr);
            keptWaters[w].type = "O";
            keptWaters[w].occ = 1.0;
            waterX += keptWaters[w].x;
            waterY += keptWaters[w].y;
            waterZ += keptWaters[w].z;
        }

        aveX= waterX*invwaterCount;
        aveY= waterY*invwaterCount;
        aveZ= waterZ*invwaterCount;

        for(unsigned int w=0; w < waterCount; w++){
            keptWaters[w].x += -aveX;
            keptWaters[w].y += -aveY;
            keptWaters[w].z += -aveZ;
        }

        Coords * waterCoordinate;
        rWater = new float[waterCount];
        phiWater = new float[waterCount];
        waterCosTheta = new float[waterCount];

        for(unsigned int rt=0; rt<waterCount; rt++){
            waterCoordinate = &keptWaters[rt];
            //pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp((const float &) centeredX[n], (const float &) centeredY[n], (const float &) centeredZ[n]);
            //pconvertXYZ = FUNCTIONS_RPR::xyz_to_rtp(&(*waterCoordinate).x, &(*waterCoordinate).y, &(*waterCoordinate).z);
            float * pconvertXYZ = SUGATAMA_FUNCTIONS_H::xyz_to_rtp((const float &) (*waterCoordinate).x, (const float &) (*waterCoordinate).y, (const float &) (*waterCoordinate).z);
            rWater[rt] = *pconvertXYZ;
            waterCosTheta[rt] = cosf(*(pconvertXYZ+1));
            phiWater[rt] = *(pconvertXYZ+2);
        }

        ifWaters = true;
        std::cout << "\tUsing => " << waterCount << " waters from PDB" << std::endl;
    }
} // end of constructor

void PDBModel::writeKeptWaters() const {

    int count = 1;
    for(int i=0; i< waterCount; i++){
        if (keptWaters[i].type.size() == 0){
            break;
        }

        printf("%-3s%7i%4s%5s%2s%4i     %7.3f %7.3f %7.3f  %4.2f %4.2f\n", "ATOM", count, "O", "HOH", "A", count, keptWaters[i].x, keptWaters[i].y, keptWaters[i].z, keptWaters[i].occ, 100.0 );
        count++;
    }
}


void PDBModel::writeCenteredCoordinatesToFile(std::string name) const {

    std::string residue_index;

    name = name + ".pdb";
    //const char * outputFileName = name.c_str() ;
    //const char * originalPDBFilename = this->filename.c_str();
    FILE * pFile;
    pFile = fopen(name.c_str(), "w");

    fprintf(pFile,"REMARK  CENTERED COORDINATES : %s\n", this->filename.c_str());
    for (unsigned int n=0; n < totalAtoms; n++) {
        residue_index = boost::lexical_cast<std::string>(resID[n]);
        //fprintf(pFile, "%-3s%7i%4s%5s%2s%4s     %7.3f %7.3f %7.3f  1.00 100.00\n", "ATOM", n+1, trimmedAtomType[n].c_str(), resi[n].c_str(), chainID[n].c_str(), residue_index.c_str(), centeredX[n], centeredY[n], centeredZ[n] );
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, trimmedAtomType[n].c_str(), resi[n].c_str(), chainID[n].c_str(), residue_index.c_str(), centeredX[n], centeredY[n], centeredZ[n] );
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
}

/**
 * copy constructor
 * @param e2
 */
PDBModel::PDBModel(const PDBModel &e2) {
    filename = e2.filename;
    ifRNA = e2.ifRNA;
    ifWaters = e2.ifWaters;
    totalAtoms = e2.totalAtoms;
    totalResidues = e2.totalResidues;
    waterCount = e2.waterCount;
    watersPerResidue = e2.watersPerResidue;
    // Following vectors will all be same length
    atomType = std::vector<std::string>(e2.atomType.begin(), e2.atomType.end());
    resi = std::vector<std::string>(e2.resi.begin(), e2.resi.end());
    trimmedAtomType = std::vector<std::string>(e2.trimmedAtomType.begin(), e2.trimmedAtomType.end());
    trimmedResi = std::vector<std::string>(e2.trimmedResi.begin(), e2.trimmedResi.end());
    chainID = std::vector<std::string>(e2.chainID.begin(),e2.chainID.end());
    waterLines = std::vector<std::string>(e2.waterLines.begin(), e2.waterLines.end());
    // waters
    resID = std::vector<unsigned int>(e2.resID.begin(), e2.resID.end());
    std::vector < float > x,y,z,atomVolume,atomicRadii;
    x = std::vector < float > (e2.x.begin(), e2.x.end());
    y = std::vector < float > (e2.y.begin(), e2.y.end());
    z = std::vector < float > (e2.z.begin(), e2.z.end());
    atomVolume = std::vector < float > (e2.atomVolume.begin(), e2.atomVolume.end());
    atomicRadii = std::vector < float > (e2.atomicRadii.begin(), e2.atomicRadii.end());


    Coords * keptWaters = nullptr;

    atomList = std::set<unsigned int>(atomList.begin(), atomList.end());
    volume = e2.volume;
    dmax = e2.dmax;
    cutoff = e2.cutoff;

    // setup for dynamically allocated arrays
    rThetaPhiAtomType = new float[totalAtoms*5];
    for (unsigned int i = 0; i<(totalAtoms*5); i++){
        rThetaPhiAtomType[i] = e2.rThetaPhiAtomType[i];
    }

    atomicNumbers = new unsigned int[totalAtoms];
    atomicExpVolume = new float[totalAtoms];
    centeredX = new float[totalAtoms];
    centeredY = new float[totalAtoms];
    centeredZ = new float[totalAtoms];
    for (unsigned int i = 0; i<(totalAtoms); i++){
        atomicNumbers[i] = e2.atomicNumbers[i];
        atomicExpVolume[i] = e2.atomicExpVolume[i];
        centeredX[i] = e2.centeredX[i];
        centeredY[i] = e2.centeredY[i];
        centeredZ[i] = e2.centeredZ[i];
    }

    keptWaters = new Coords[e2.waterCount];
    rWater = new float[waterCount];
    phiWater = new float[waterCount];
    waterCosTheta = new float[waterCount];
    for (unsigned int i = 0; i<waterCount; i++){
        keptWaters[i] = Coords(e2.keptWaters[i]);
        rWater[i] = e2.rWater[i];
        phiWater[i] = e2.phiWater[i];
        waterCosTheta[i] = e2.waterCosTheta[i];
    }
}



