//
// Created by ttbond on 18-10-15.
//

#ifndef SRCFIND_BACKGROUNDGENERATOR_H
#define SRCFIND_BACKGROUNDGENERATOR_H
#include<vector>
#include<string>
#include"bedFile.h"
#include"faFile.h"
#include"basicInfo.h"


class backgroundGenerator {
public:
    backgroundGenerator(char *agctArr,char *agctArr2,double _selectP=1.0/30000.0,long long _length=-1,int _chr=-1);
    double selecP;
    long long length;
    int chr;
    static long long *lengthArr;
};


#endif //SRCFIND_BACKGROUNDGENERATOR_H
