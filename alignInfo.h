//
// Created by ttbond on 19-9-23.
//

#ifndef SRCFIND_ALIGNINFO_H
#define SRCFIND_ALIGNINFO_H
#include <string>
#include <time.h>
#include "basicInfo.h"
using namespace std;

class alignInfo {
public:
    alignInfo(int _chr1,long long _st1,long long _ed1,long long st2,long long ed2,char* qname,int flag);
    const alignInfo& operator =(const alignInfo &source);
    alignInfo(const alignInfo &source);
    basicInfo refInfo,readInfo;
    int flag;
};


#endif //SRCFIND_ALIGNINFO_H
