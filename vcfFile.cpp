//
// Created by ttbond on 18-9-16.
//

#include "vcfFile.h"
vcfFile::vcfFile(char *fileName){
    printf("loading vcfFile from %s...\n",fileName);
    basicVCF=new oneVCF(fileName);
    FILE *fp=fopen(fileName,"r");
    int maxPerLineLength=150000;
    char cache[maxPerLineLength+5];
    while(fgets(cache,maxPerLineLength,fp)!=NULL){
        if((*cache)=='c'){
            if(strlen(cache)==maxPerLineLength){
                continue;
            }
            vcfs.push_back(oneVCF(cache));
        }
    }
    vcfNum=vcfs.size();
    printf("successfully load %lld svs\n",vcfNum);
}
vcfFile::vcfFile(std::vector<oneVCF> vcfsVec){
    oneVCF testIfHaveBasicVcf;
    if(testIfHaveBasicVcf.knowFormat==-1){
        printf("there are no oneVCF before");
        exit(-1);
    }
    vcfs=vcfsVec;
    vcfNum=vcfs.size();
}
std::vector<oneVCF> vcfFile::selectByChr(int chrNum,bool notReverse){
    printf("selecting by chr from %lld svs...\n",vcfNum);
    std::vector<oneVCF>::iterator it;
    std::vector<oneVCF> rel;
    for(it=vcfs.begin();it!=vcfs.end();it++){
        if(((*it).chr==chrNum)==notReverse){
            rel.push_back((*it));
        }
    }
    printf("successfully select %lld svs\n",rel.size());
    return rel;
}

//notReverse is false means to except sv
std::vector<oneVCF> vcfFile::selectBySv(svType sv,bool notReverse){
    printf("selecting by svType from %lld svs...\n",vcfNum);
    std::vector<oneVCF>::iterator it;
    std::vector<oneVCF> rel;
    for(it=vcfs.begin();it!=vcfs.end();it++){
        if(((*it).sv==sv)==notReverse){
            rel.push_back((*it));
        }
    }
    printf("successfully select %lld svs\n",rel.size());
    return rel;
}


void vcfFile::getDetetionRegions2(std::vector<int>&chr,std::vector<long long>&st,std::vector<long long>&st2,std::vector<long long>&ed,std::vector<long long>&ed2,long long maxIdx){
    std::vector<oneVCF>::iterator it;
    for(it=vcfs.begin();it!=vcfs.end();it++){
        oneVCF &tmp=(*it);
        int width=std::min(std::max((long long)1000,tmp.length*2),(long long)5000);
        chr.push_back(tmp.chr);
        st.push_back(std::max((long long)0,tmp.st-width));
        ed.push_back(std::min(maxIdx,tmp.ed+width));
    }
}

void vcfFile::getDetetionRegions(std::vector<int>&chr,std::vector<long long>&st,std::vector<long long>&ed,std::vector<svType>&svs,long long maxIdx){
    std::vector<oneVCF>::iterator it;
    for(it=vcfs.begin();it!=vcfs.end();it++){
        oneVCF &tmp=(*it);
        int width=200;
        chr.push_back(tmp.chr);
        st.push_back(std::max((long long)0,tmp.st-width-1));
        ed.push_back(std::min(maxIdx,tmp.ed+width-1));
        svs.push_back(tmp.sv);
    }
}

void vcfFile::getDetRegBpSides(char side,long long size,std::vector<int>&chr,std::vector<long long>&st,std::vector<long long>&ed,std::vector<svType>&svs,long long maxIdx){
    if(side!='L'&&side!='R'){
        printf("para side in getDetRegBpSide wrong!\n");
        exit(-1);
    }
    std::vector<oneVCF>::iterator it;
    for(it=vcfs.begin();it!=vcfs.end();it++){
        oneVCF &tmp=(*it);
        chr.push_back(tmp.chr);
        if(side=='L') {
            st.push_back(std::max((long long) 0, tmp.st - size - 1));
            ed.push_back(std::min(tmp.st + size - 1, tmp.ed));
        }
        if(side=='R'){
            st.push_back(std::max(tmp.st,tmp.ed-size-1));
            ed.push_back(std::min(tmp.ed+size-1,maxIdx));
        }
        svs.push_back(tmp.sv);
    }
}