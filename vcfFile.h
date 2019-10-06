//
// Created by ttbond on 18-9-16.
//

#ifndef SRCFIND_VCFFILE_H
#define SRCFIND_VCFFILE_H
#include<vector>
#include<stdlib.h>
#include<algorithm>

#include"oneVCF.h"




class vcfFile {
public:
    vcfFile(char *fileName);
    vcfFile(std::vector<oneVCF> vcfsVec);
    std::vector<oneVCF> selectByChr(int chrNum,bool notReverse=true);
    std::vector<oneVCF> selectBySv(svType sv,bool notReverse=true);
    void getDetetionRegions(std::vector<int>&chr,std::vector<long long>&st,std::vector<long long>&ed,std::vector<svType>&svs,long long maxIdx);
    void getDetetionRegions2(std::vector<int>&chr,std::vector<long long>&st,std::vector<long long>&st2,std::vector<long long>&ed,std::vector<long long>&ed2,long long maxIdx);
    void getDetRegBpSides(char side,long long size,std::vector<int>&chr,std::vector<long long>&st,std::vector<long long>&ed,std::vector<svType>&svs,long long maxIdx);
    oneVCF *basicVCF;
    std::vector<oneVCF> vcfs;
    int vcfNum;
};


#endif //SRCFIND_VCFFILE_H
