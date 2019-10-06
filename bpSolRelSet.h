//
// Created by ttbond on 18-9-29.
//

#ifndef SRCFIND_BPSOLRELSET_H
#define SRCFIND_BPSOLRELSET_H
#include<vector>
#include"bpSolRel.h"

class bpSolRelSet {
public:
    std::vector<bpSolRel> bpSols;

    bpSolRelSet();
    bpSolRelSet(std::vector<bpSolRel> &source);
    void addBpSol(bpSolRel &tmp);
    void printMax(char side,char type,FILE *fp=NULL);
    void printMaxDis(char side,char type,FILE *fp=NULL);
    void printDis2bp(char side,char type,FILE *fp=NULL);

};


#endif //SRCFIND_BPSOLRELSET_H
