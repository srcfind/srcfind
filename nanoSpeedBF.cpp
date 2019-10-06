//
// Created by ttbond on 19-6-27.
//
#include <set>
#include "nanoSpeedBF.h"


nanoSpeedBF::nanoSpeedBF():bedFile() {
    staNum=NULL;
}

void nanoSpeedBF::loadNanoSpeedBF(const vector<basicInfo> &source) {

    myRegion=source;
}

void nanoSpeedBF::staRepRegions() {
    staRel.clear();
    delete []staNum;
    int siz=myRegion.size();
    staNum=new int[siz];
    memset(staNum,0,sizeof(int)*siz);
    for(int i=0;i<siz;i++){
        if(myRegion[i].isIllegal()){
            continue;
        }
        for(int j=0;j<i;j++){
            if(myRegion[j].isIllegal()){
                continue;
            }
            if(myRegion[i].overlapRate(myRegion[j])>0){
                staRel.insert(myRegion[i]);
                staRel.insert(myRegion[i]);
                staNum[i]++,staNum[j]++;
            }
        }
    }
    int num=0;
    for(int i=0;i<siz;i++){
        if(staNum[i]>0){
            num++;
            myRegion[i].printMe();
        }
    }
    printf("total %d regions found\n",num);
}


