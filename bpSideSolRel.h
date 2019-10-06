//
// Created by ttbond on 18-9-28.
//

#ifndef SRCFIND_BPSIDESOLREL_H
#define SRCFIND_BPSIDESOLREL_H
#include"solRel.h"
#include "bpDecReg.h"

class bpSideSolRel:public solRel {
public:
    bpSideSolRel(detectRegion &source);
    bpSideSolRel(bpSideSolRel &source);
    long long getDis2bp(double *st,double *ed,double maxVal);
    long long getMaxDis(double *st,double *ed,double maxVal);
    long long dis2bpRev,dis2bpDir,dis2bpMir;
    long long maxDisRev,maxDisDir,maxDisMir;
    long long bpPos;
};


#endif //SRCFIND_BPSIDESOLREL_H
