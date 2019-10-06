//
// Created by ttbond on 18-9-28.
//

#include "bpSideSolRel.h"


bpSideSolRel::bpSideSolRel(detectRegion &source):
                    solRel(source){
    bpDecReg &_source=(bpDecReg &)source;
    bpPos=_source.bpPos;
    dis2bpRev=getDis2bp(source.reverseComScore,source.reverseComScore+source.length,maxRevComScore);
    dis2bpDir=getDis2bp(source.directRepScore,source.directRepScore+source.length,maxDirRepScore);
    dis2bpMir=getDis2bp(source.mirrorRepScore,source.mirrorRepScore+source.length,maxMirRepScore);
    maxDisRev=getMaxDis(source.reverseComScore,source.reverseComScore+source.length,maxRevComScore);
    maxDisDir=getMaxDis(source.directRepScore,source.directRepScore+source.length,maxDirRepScore);
    maxDisMir=getMaxDis(source.mirrorRepScore,source.mirrorRepScore+source.length,maxMirRepScore);
}


bpSideSolRel::bpSideSolRel(bpSideSolRel &source):
                    solRel((solRel&)source){
    bpPos=source.bpPos;
    dis2bpRev=source.dis2bpRev;
    dis2bpDir=source.dis2bpDir;
    dis2bpMir=source.dis2bpMir;
    maxDisRev=source.maxDisRev;
    maxDisDir=source.maxDisDir;
    maxDisMir=source.maxDisMir;
}


long long bpSideSolRel::getDis2bp(double *st, double *ed,double maxVal) {
    long long rel=1e15;
    for(double *i=st;i<ed;i++){
        if((*i)==maxVal){
            long long nowPos=i-st+(info->st);
            long long partialPos=nowPos>bpPos?nowPos-bpPos:bpPos-nowPos;
            rel=rel>nowPos?nowPos:rel;
        }
    }
    return rel;
}

long long bpSideSolRel::getMaxDis(double *st,double *ed,double maxVal){
    long long rel=1e15;
    double tmp1=-1,tmp2=-1;
    double *pos1,*pos2=NULL,*pos3=NULL;
    for(double *i=st;i<ed;i++){
        if((*i)==maxVal){
            pos1=i+(i-st);
            for(double *j=i-1;j>=st;j--){
                if((*j)>(*(j+1))){
                    for(double *k=j;k>=st;k--){
                        if(tmp1<(*k)){
                            tmp1=*k;
                            pos2=k;
                        }
                    }
                    break;
                }
            }
            for(double *j=i+1;j<ed;j++){
                if((*j)>(*(j-1))){
                    for(double *k=j;k<ed;k++){
                        if(tmp2<(*k)){
                            tmp2=*k;
                            pos3=k;
                        }
                    }
                    break;
                }
            }
            break;
        }
    }
    if(pos2!=NULL){
        if(pos3!=NULL){
            if((*pos2)>(*pos3)){
                rel=pos1-pos2;
            }
            else{
                rel=pos3-pos2;
            }
        }
        else{
            rel=pos1-pos2;
        }
    }
    else{
        if(pos3!=NULL){
            rel=pos3-pos2;
        }
        else{
            rel=-1;
        }
    }
    return rel;
}