//
// Created by ttbond on 18-11-15.
//

#ifndef SRCFIND_BASICINFO_H
#define SRCFIND_BASICINFO_H


#include<stdio.h>
#include<string.h>
#include"ttbond_math.h"
using namespace std;
class basicInfo {
public:
    basicInfo();
    basicInfo(const basicInfo &source);
    basicInfo(int _chr,long long _st,long long _ed,svType _sv=OTHER);
    basicInfo(int _chr,long long _st,long long _ed,const char * qname,svType _sv=OTHER);
    ~basicInfo();
    bool operator <(const basicInfo &right)const;
    bool operator ==(basicInfo &right);
    basicInfo &operator =(const basicInfo &source);
    void printMe(FILE *fp=NULL);
    void printRange(FILE *fp=NULL);
    long long getLength();
    char *getAgct(char *ref);
    void deleteCache();
    double overlapRate(const basicInfo &tmp);
    bool isIllegal();
    bool myStGreaterThan(const basicInfo &tmp);
    bool myEdSmallerThan(const basicInfo &tmp);
    bool myEdLargerThanEdOf(const basicInfo &tmp);
    bool myEdSmallerThanStOf(const basicInfo &tmp);
    int chr;
    long long st,ed;
    svType sv;
    char *cache;
};


#endif //SRCFIND_BASICINFO_H
