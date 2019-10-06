//
// Created by ttbond on 18-9-29.
//

#ifndef SRCFIND_BPSOLREL_H
#define SRCFIND_BPSOLREL_H
#include"bpSideSolRel.h"

class bpSolRel {
public:
    bpSolRel(bpSideSolRel &_left,bpSideSolRel &_right);
    bpSideSolRel *left,*right;
    long long getDis2bp(char side,char type);
    long long getMaxDis(char side,char type);
    double getMax(char side,char type);
};


#endif //SRCFIND_BPSOLREL_H
