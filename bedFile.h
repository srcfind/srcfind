//
// Created by ttbond on 18-10-11.
//

#ifndef SRCFIND_BEDFILE_H
#define SRCFIND_BEDFILE_H

#include<vector>
#include<algorithm>
#include<iostream>
#include"basicInfo.h"
#include"detectRegion.h"
#include"ttbond_fa.h"
#include"solRel.h"
using namespace std;
enum fileT {BEDF,NONBF};

class bedFile {
public:

    bedFile(const char *fileName,char *_agct);
    bedFile(std::vector<basicInfo> _source,char *_agct);
    bedFile();
    bool load(const char * fileName,fileT fileType);
    bool loadNonBF(FILE *fp);
    //~bedFile();
    void printSingleRcdScore(char *fileName);
    void solAndPrint(char *fileName,char *writeType);
    void detectPotentialDSB();
    bool atRegionEnd();
    void init();
    basicInfo getNowRegion();
    char *getChrAgct(basicInfo &tar);
    void sortMyRegion();
    std::vector<basicInfo> myRegion;
    std::vector<solRel> myRel;
    char *agct;
    int nowPos;
    int nowRefChr;
};


#endif //SRCFIND_BEDFILE_H
