//
// Created by ttbond on 18-12-1.
//

#ifndef SRCFIND_NONBPOTDETECTOR_H
#define SRCFIND_NONBPOTDETECTOR_H
#include<vector>
#include<algorithm>
#include"bedFile.h"
#include"nonB.h"

class nonbPotDetector {
public:
    nonbPotDetector(char *bedFileName,char *agct);
    void detect();
    void printF(const char *fileName);
    void printF_point(const char *fileName);
    ~nonbPotDetector();
    bedFile *bedF;
    std::vector<nonB> detectRel;
    int maxSep;
};


#endif //SRCFIND_NONBPOTDETECTOR_H
