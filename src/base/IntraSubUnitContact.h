//
// Created by xos81802 on 20/09/2018.
//

#ifndef SUGATAMA_INTRASUBUNITCONTACT_H
#define SUGATAMA_INTRASUBUNITCONTACT_H

class IntraSubUnitContact {
    unsigned int subUnit1Index, subUnit2Index;

public:
    IntraSubUnitContact(unsigned int one, unsigned int two) {this->subUnit1Index = one; this->subUnit2Index=two;}
    int getSubUnit1(){return subUnit1Index;}
    int getSubUnit2(){return subUnit2Index;}
};

#endif //SUGATAMA_INTRASUBUNITCONTACT_H
