//
// Created by ttbond on 18-9-26.
//
#ifndef SRCFIND_SOLREL_CPP
#define SRCFIND_SOLREL_CPP
#include "solRel.h"

solRel::solRel() {

}

solRel::solRel(detectRegion &source) {
    maxSelectP=0.8;
    info=new basicInfo(source.chr,source.st,source.ed,source.mySv);
    maxRevComScore=getMax(source.reverseComScore,source.reverseComScore+source.length);
    maxDirRepScore=getMax(source.directRepScore,source.directRepScore+source.length);
    maxMirRepScore=getMax(source.mirrorRepScore,source.mirrorRepScore+source.length);
    lengthMaxRev=getSpeLen(source.reverseComScore,source.reverseComScore+source.length,maxRevComScore*maxSelectP);
    lengthMaxDir=getSpeLen(source.directRepScore,source.directRepScore+source.length,maxDirRepScore*maxSelectP);
    lengthMaxMir=getSpeLen(source.mirrorRepScore,source.mirrorRepScore+source.length,maxMirRepScore*maxSelectP);
    sumRevComScore=getSum(source.reverseComScore,source.reverseComScore+source.length);
    sumDirRepScore=getSum(source.directRepScore,source.directRepScore+source.length);
    sumMirRepScore=getSum(source.mirrorRepScore,source.mirrorRepScore+source.length);
    meanRevComScore=getMean(source.reverseComScore,source.reverseComScore+source.length);
    meanDirRepScore=getMean(source.directRepScore,source.directRepScore+source.length);
    meanMirRepScore=getMean(source.mirrorRepScore,source.mirrorRepScore+source.length);
    /*if(meanMirRepScore<0.001){
        info->printMe();
        printf("%s\n",source.agct);
        getchar();
    }*/
}

solRel::solRel(const solRel &source) {
    maxSelectP=0.8;
    info=new basicInfo(*(source.info));
    maxRevComScore=source.maxRevComScore;
    maxDirRepScore=source.maxDirRepScore;
    maxMirRepScore=source.maxMirRepScore;
    lengthMaxRev=source.lengthMaxRev;
    lengthMaxDir=source.lengthMaxDir;
    lengthMaxMir=source.lengthMaxMir;
    sumRevComScore=source.sumRevComScore;
    sumDirRepScore=source.sumDirRepScore;
    sumMirRepScore=source.sumMirRepScore;
    meanRevComScore=source.meanRevComScore;
    meanDirRepScore=source.meanDirRepScore;
    meanMirRepScore=source.meanMirRepScore;
    /*if(meanMirRepScore==0.00){
        info->printMe();
        getchar();
    }*/
}

solRel::~solRel(){
    delete info;
}

double solRel::getMax(double *st, double *ed) {
    double rel=-1;
    for(double *i=st;i<ed;i++){
        if((*i)>rel){
            rel=(*i);
         }
    }
    return rel;
}

double solRel::getSum(double *st, double *ed) {
    double sum=0;
    for(double *i=st;i<ed;i++){
        sum+=(*i);
    }
    return sum;
}

long long solRel::getSpeLen(double *st, double *ed, double bar) {
    long long len=0;
    for(double *i=st;i<ed;i++){
        if((*i)>=bar){
            len++;
        }
    }
    return len;
}

double solRel::getMean(double *st, double *ed) {
    double sum=getSum(st,ed);
    return sum/(ed-st);
}

void solRel::printMe(FILE *fp){
    fprintf(fp,"%.2lf\t%.2lf\t%.2lf\t",meanDirRepScore,meanMirRepScore,meanRevComScore);
}

#endif