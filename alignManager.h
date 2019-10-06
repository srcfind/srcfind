//
// Created by ttbond on 19-8-21.
//

#ifndef SRCFIND_ALIGNMANAGER_H
#define SRCFIND_ALIGNMANAGER_H
#include <iostream>
#include <map>
#include <vector>
#include <cstring>
#include <Rcpp.h>
#include <RInside.h>
#include <stdlib.h>
#include "omp.h"
#include "oneFastq.h"
#include "bamFile.h"
#include "faFile.h"
#include "detectRegion.h"
#include "alignInfo.h"
using namespace std;

class alignManager {
public:
    alignManager(oneFastq &read,char *refFile, const char *bamFileName,char *refCache,const char *outF,bool plot=false);
    alignManager(oneFastq &read,const char *bamFileName,const char *outF);
    alignManager(oneFastq &read);
    inline int code(const char tmp,bool com=false);
    void generateKmer(map<long long, vector<long long> > &rel,const char *seq,bool com=false);
    vector< pair<int,int> > getXY(bool refCom,bool readCom);
    void getCluster(vector< pair<int,int> > &pos);
    void getCluster(vector< pair<int,int> > &pos,int k,int *&index,int *&num);
    int getClusterId(int tar,int *fa);
    long long getClusterRefLen(int *Index,int *numIndex,vector<pair<int, int>> &pos,int refLen);
    pair<long long,long long> getClusterSamLen(int *Index,int *numIndex,vector<pair<int, int>> &pos,int refLen);
    void printPos(int *Index,int *numIndex,vector<pair<int, int>> &pos,int type);
    vector<pair<int,int> > getPolePos(vector<pair<int,int> > &posSet);
    void printPos(vector <pair <int,int>> &pos,int type);
    void storeMainRel(int *Index,int *numIndex,vector<pair<int, int>> &rel,vector<pair<int,int> > &mainRel);
    void initSegTree(int range);
    void dfsInitSegTree(int l,int r,int id);
    void dfsInsSeg(int l,int r,int id);
    int getTotalSegLen(int l,int r,int id);
    void copyCluster(int segId,int *index,int *num,int * &des_index,int * &des_num,int relLen);
    vector< pair<int,int> > getMaxClusterOnSam(int *Index,int *numIndex,vector<pair<int, int>> &pos);
    pair<double,double> getLinearRegressionPara(vector<pair<int,int> > &pos);
    void plotDotplot(vector<pair<int,int> > pointSet1,int colTyp1,vector<pair<int,int> >pointSet2,int colTyp2,string &qname,int refLen,int readLen,int id);
    void plotRefBoundary(long long refBoundaryPos,char *refCache,string &qname);
    void getEndBoundaryInfo(int chr,long long st,long long ed,oneFastq &read);
    void getDotplot(char *refFile, const char *bamFileName,char *refCache,const char *outD,const char *wd);
    void getDotplot(basicInfo refPos,const char *refCache,const char *wd);
    basicInfo getRefPosFromAlign(alignInfo source,alignInfo tar,bool rev);
    ~alignManager();
    oneFastq myFastq;
    const char *seq1,*seq2;
    char *cache;
    int *posId,*idNum,*posId2,*idNum2,*secPosId,*secIdNum;
    int *segPosId[3],*segIdNum[3];
    map<long long,vector<long long> > map1,map2;
    int *segTreeL,*segTreeR,*segNodesHave,*segNodesFull;
    static RInside R;
    clock_t bornT;
    long long refBoundaryPos;
};


#endif //SRCFIND_ALIGNMANAGER_H
