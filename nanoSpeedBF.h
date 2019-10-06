//
// Created by ttbond on 19-6-27.
//

#ifndef SRCFIND_NANOSPEEDBF_H
#define SRCFIND_NANOSPEEDBF_H
#include <vector>
#include <set>
#include "bedFile.h"
#include "basicInfo.h"
using namespace std;

class nanoSpeedBF: bedFile{
public:
    nanoSpeedBF();
    void loadNanoSpeedBF(const vector<basicInfo> &source);
    void staRepRegions();
    set<basicInfo> staRel;
    int *staNum;
};


#endif //SRCFIND_NANOSPEEDBF_H
