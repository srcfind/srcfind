//
// Created by ttbond on 18-10-15.
//

#ifndef SRCFIND_FAFILE_H
#define SRCFIND_FAFILE_H
#include<vector>
#include<set>
#include<random>
#include<time.h>
#include"basicInfo.h"

class faFile {
public:
    faFile(char *_agct);
    inline bool isRegularChr(const char *str);
    //if the symbol is 'A/G/C/T/N' return true
    inline bool isNormalAgctSymbol(const char &tmp);
    //if the specific chrNum is found then it's sequencing will be loaded to array agct
    // and return the length of the chr .... if not exist return -1
    long long loadAgctByChr(int chrNum,char *fileName);
    void getValidatedRegionByChr(int chrNum,char *fileName);
    void getRandomBackgroundRegion(long long length,double percent,std::vector<basicInfo> &dest);
    void getRandomPos(int posNum,std::vector<basicInfo> &dest);
    char *agct;
    long long agctLen,totalValidatedNum;
    std::vector<basicInfo> validatedRegion;
};


#endif //SRCFIND_FAFILE_H
