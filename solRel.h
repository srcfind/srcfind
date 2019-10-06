//
// Created by ttbond on 18-9-26.
//

#ifndef SRCFIND_SOLREL_H
#define SRCFIND_SOLREL_H
#include"basicInfo.h"
#include"detectRegion.h"


class solRel {
public:
    solRel();
    solRel(const solRel &source);
    solRel(detectRegion &source);
    ~solRel();
    void printMe(FILE *fp);

    double getMax(double *st,double *ed);
    double getSum(double *st,double *ed);
    double getMean(double *st,double *ed);
    long long getSpeLen(double *st,double *ed,double bar);
    basicInfo *info;
    double sumRevComScore,sumDirRepScore,sumMirRepScore;
    double maxRevComScore,maxDirRepScore,maxMirRepScore;
    double meanRevComScore,meanDirRepScore,meanMirRepScore;
    long long lengthMaxRev,lengthMaxDir,lengthMaxMir;
    double maxSelectP;
};


#endif //SRCFIND_SOLREL_H
