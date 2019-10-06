//
// Created by ttbond on 18-10-15.
//
#ifndef SRCFIND_BACKGROUNDGENERATOR_CPP
#define SRCFIND_BACKGROUNDGENERATOR_CPP
#include "backgroundGenerator.h"

long long * backgroundGenerator::lengthArr=new long long[25]{1200,1300,1400,1500,1600,1700,1800,1900,2100,2200,2300,2400,2500,2600,2700,2800,2900};
backgroundGenerator::backgroundGenerator(char *agctArr,char *agctArr2,double _selectP,long long _length,int _chr){
    length=_length;
    chr=_chr;
    selecP=_selectP;
    std::vector<basicInfo> tmpVec;
    if(length==-1){
        for(int i=0;i<17;i++){
            long long nowLength=lengthArr[i];
            if(chr==-1){
                for(int j=1;j<24;j++){
                    faFile testF(agctArr);
                    testF.getValidatedRegionByChr(j,"GRCh38.d1.vd1.fa");
                    tmpVec.clear();
                    testF.getRandomBackgroundRegion(nowLength,selecP,tmpVec);
                    bedFile bf(tmpVec,agctArr2);
                    char cache[50]="./bg/backGround_";
                    if(j==1){
                        bf.solAndPrint(strcat(cache,std::to_string((int)nowLength).c_str()),"w");
                    }
                    else{
                        bf.solAndPrint(strcat(cache,std::to_string((int)nowLength).c_str()),"a");
                    }
                }
            }
        }
    }

}



#endif