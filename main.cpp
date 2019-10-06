//
// Created by ttbond on 18-10-9.
//
#include<stdio.h>
#include<string.h>
#include<vector>
#include<stdlib.h>

#include"geneRegionFile.h"

#include"basicInfo.h"
#include"bedFile.h"
#include"faFile.h"
#include"backgroundGenerator.h"
#include"agctTree.h"
#include"nonbPotDetector.h"

char agct[300000000];
char agctForBedFile[300000000];
std::vector<basicInfo> tmpVec;
int main(){

    nonbPotDetector det("test.bed",agct);
    det.detect();
    det.printF_point("testF.dat");
    //backgroundGenerator bg(agct,agctForBedFile);


    /*
    geneRegionFile grf("./geneExons/geneExonNoCancer.bed",agct);
    grf.getSeqSta();
    */


    //bedFile bf("./dsb/dsb.bed",agct);
    //bf.solAndPrint("./dsb/allScore.scr","w");

    //bf.printSingleRcdScore("./dsb/single");
}
