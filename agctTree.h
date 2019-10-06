//
// Created by ttbond on 18-10-17.
//

#ifndef SRCFIND_AGCTTREE_H
#define SRCFIND_AGCTTREE_H
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<algorithm>
#include<string.h>
#include"agctNode.h"
#include"gcSta.h"
class agctTree {
public:

    agctTree();
    ~agctTree();
    void release(agctNode *tar,int level);
    void release(char *&tmp);
    void addSeq(int *lenSt,int lenEd,char *seq);
    void addSeq(agctNode *nowNode,char *from,int len,int level);
    inline int nu2int(char *from);
    inline char int2nu(int from);
    void printMostSeq(int rank,int _leastLen,FILE *fp);
    void getTotalNum(std::vector<int> &rel,agctNode *nowNode,int level);
    void selectByFre(agctNode *nowNode,int level);
    agctNode *root;
    int maxLevel,leastLen,minSelectFre,minSelectLen;
    char *charCache;
    gcSta *myGcSta;
    FILE *fpCache;
};


#endif //SRCFIND_AGCTTREE_H
