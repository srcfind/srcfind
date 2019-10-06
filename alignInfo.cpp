//
// Created by ttbond on 19-9-23.
//

#include "alignInfo.h"

alignInfo::alignInfo(int _chr1, long long _st1, long long _ed1, long long _st2, long long _ed2, char *qname, int _flag) {
    refInfo=basicInfo(_chr1,_st1,_ed1);
    readInfo=basicInfo(-1,_st2,_ed2,qname);
    flag=_flag;
}

const alignInfo &alignInfo::operator=(const alignInfo &source) {
    refInfo=source.refInfo;
    readInfo=source.readInfo;
    flag=source.flag;
    return *this;
}

alignInfo::alignInfo(const alignInfo &source) {
    refInfo=source.refInfo;
    readInfo=source.readInfo;
    flag=source.flag;
}


