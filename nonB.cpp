//
// Created by ttbond on 18-11-15.
//
#ifndef SRCFIND_NONB_CPP
#define SRCFIND_NONB_CPP

#include "nonB.h"

nonB::nonB(std::vector<basicInfo> _myRegion,char *agct,nonBTyp _myNonBTyp){
    myRegion=_myRegion;
    myNonBTyp=_myNonBTyp;
    std::sort(myRegion.begin(),myRegion.end());
    regionNum=myRegion.size();
    for(int i=0;i<regionNum;i++){
        myRegion[i].getAgct(agct);
    }
    for(int i=0;i<regionNum-1;i++){
        for(int j=i+1;j<regionNum;j++){
            char *stri=myRegion[i].cache;
            char *strj=myRegion[j].cache;
            int regionLeni=strlen(stri);
            int regionLenj=strlen(strj);
            if(regionLeni!=regionLenj){
                printf("Error in nonB regions from length\n");
                exit(0);
            }
            int regionLen=regionLeni;
            for(int k=0;k<regionLen;k++){
                if(stri[k]!=strj[k]){
                    relativePos.push_back(std::make_pair(k,regionLen-1-k));
                    misMachPoint.push_back(basicInfo(myRegion[i].chr,myRegion[i].st+k,myRegion[i].st+k));
                    misMachPoint.push_back(basicInfo(myRegion[j].chr,myRegion[j].st+k,myRegion[j].st+k));
                    std::vector<char> tmp;
                    tmp.push_back(stri[k]);
                    tmp.push_back(strj[k]);
                    misMachAgct.push_back(tmp);
                }
            }
        }
    }

    for(int i=0;i<regionNum;i++){
        myRegion[i].deleteCache();
    }
}

bool nonB::isBrinkDiverge(){
    int siz=relativePos.size();
    for(int i=0;i<siz;i++){
        if(relativePos[i].first==0||relativePos[i].second==0){
            return true;
        }
    }
    return false;
}

int nonB::getMisNum() {
    return relativePos.size();
}

bool nonB::operator ==(const nonB &source){
    if(myRegion.size()!=source.myRegion.size()){
        return false;
    }
    int siz=myRegion.size();
    for(int i=0;i<siz;i++){
        if(!(myRegion[i]==(basicInfo &)source.myRegion[i])){
            return false;
        }
    }
    return true;
}

bool nonB::operator <(const nonB &source){
    if(myRegion.size()!=source.myRegion.size()){
        return myRegion.size()<source.myRegion.size();
    }
    else{
        int siz=myRegion.size();
        for(int i=0;i<siz;i++){
            if(!(myRegion[i]==(basicInfo &)source.myRegion[i])){
                return myRegion[i]<(basicInfo &)source.myRegion[i];
            }
        }
    }
    return true;
}

void nonB::printMe(FILE *relFp){
    int siz=myRegion.size();
    for(int i=0;i<siz;i++){
        myRegion[i].printMe(relFp);
    }
    for(int i=0;i<siz;i++){
        misMachPoint[i].printMe(relFp);
    }
}
#endif
