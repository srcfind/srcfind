//
// Created by ttbond on 18-12-1.
//
#ifndef SRCFIND_NONBPOTDETECTOR_CPP
#define SRCFIND_NONBPOTDETECTOR_CPP

#include "nonbPotDetector.h"
nonbPotDetector::nonbPotDetector(char *bedFileName,char *agct){
    bedF=new bedFile(bedFileName,agct);
    maxSep=7;
}
nonbPotDetector::~nonbPotDetector(){
    if(bedF!=NULL){
        delete bedF;
    }
}
void nonbPotDetector::detect() {
    bedF->init();
    std::vector<std::vector<basicInfo> >tmpRel;
    while(!bedF->atRegionEnd()){
        basicInfo nowR=bedF->getNowRegion();
        detectRegion region(bedF->getChrAgct(nowR),nowR);
        tmpRel=region.getPotDirPos();
    }
    int tsiz=tmpRel.size();
    int jishu=0;
    for(int i=0;i<tsiz;i++){
        std::vector<basicInfo> &tvec=tmpRel[i];
        sort(tvec.begin(),tvec.end());
        int siz=tvec.size();
        for(int j=0;j<siz;j++){
            for(int k=j+1;k<siz;k++){
                if(tvec[j].ed>=tvec[k].st){
                    continue;
                }
                if(tvec[k].st-tvec[j].ed>=maxSep){
                    break;
                }
                jishu=jishu+1;
                std::vector<basicInfo> toConsNonb;
                toConsNonb.push_back(tvec[j]);
                toConsNonb.push_back(tvec[k]);
                nonB tmpNonb(toConsNonb,bedF->getChrAgct(toConsNonb[0]),DIRREP);
                if(tmpNonb.getMisNum()!=1||tmpNonb.isBrinkDiverge()){
                    continue;
                }
                detectRel.push_back(tmpNonb);
            }
        }
    }
    sort(detectRel.begin(),detectRel.end());
    std::vector<nonB>:: iterator it;
    it=unique(detectRel.begin(),detectRel.end());
    while(detectRel.end()!=it){
        detectRel.pop_back();
    }
    printf("after unique:%d\n",detectRel.size());
    puts("finish");
}
void nonbPotDetector::printF(const char *fileName){
    int siz=detectRel.size();
    FILE *relF=fopen(fileName,"w");
    for(int i=0;i<siz;i++){
        detectRel[i].printMe(relF);
        fprintf(relF,"\n");
    }
    fclose(relF);
}
void nonbPotDetector::printF_point(const char *fileName){
    std::vector<basicInfo> tmp;
    int siz=detectRel.size();
    for(int i=0;i<siz;i++){
        tmp.insert(tmp.end(),detectRel[i].misMachPoint.begin(),detectRel[i].misMachPoint.end());
    }
    sort(tmp.begin(),tmp.end());
    std::vector<basicInfo>::iterator it;
    it=unique(tmp.begin(),tmp.end());
    while(tmp.end()!=it){
        tmp.pop_back();
    }
    FILE *fp=fopen(fileName,"w");
    siz=tmp.size();
    for(int i=0;i<siz;i++){
        tmp[i].printMe(fp);
        fprintf(fp,"\n");
    }
    fclose(fp);
}
#endif