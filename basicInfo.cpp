//
// Created by ttbond on 18-11-15.
//
#ifndef SRCFIND_BASICINFO_CPP
#define SRCFIND_BASICINFO_CPP
#include <iostream>
#include "basicInfo.h"
#include "stdlib.h"

using namespace std;


basicInfo::basicInfo(){
    chr=-1;
    st=-1;
    ed=-1;
    sv=OTHER;
    cache=NULL;
}



basicInfo::basicInfo(const basicInfo &source) {
    chr=source.chr;
    st=source.st;
    ed=source.ed;
    sv=source.sv;
    if(source.cache!=NULL){
        long long len=strlen(source.cache);
        cache=new char[len+10];
        strcpy(cache,source.cache);
    }
    else{
        cache=NULL;
    }
}

basicInfo &basicInfo::operator =(const basicInfo &source){
    chr=source.chr;
    st=source.st;
    ed=source.ed;
    sv=source.sv;
    if(cache!=NULL){
        delete []cache;
    }
    cache=NULL;
    if(source.cache!=NULL){
        long long len=strlen(source.cache);
        cache=new char[len+10];
        strcpy(cache,source.cache);
    }
    else{
        cache=NULL;
    }
    return *this;
}

basicInfo::basicInfo(int _chr,long long _st,long long _ed,svType _sv){
    chr=_chr;
    st=_st;
    ed=_ed;
    sv=_sv;
    cache=NULL;
}

basicInfo::basicInfo(int _chr, long long _st, long long _ed, const char *qname, svType _sv) {
    new (this) basicInfo(_chr,_st,_ed,_sv);
    //printf("%lld %lld\n",st,ed);
    int qnameSize=strlen(qname);
    cache=new char[qnameSize+1];
    strcpy(cache,qname);
}

bool basicInfo::operator <(const basicInfo &right) const{
    if(chr!=right.chr){
        return chr<right.chr;
    }
    else if(st!=right.st){
        return st<right.st;
    }
    else{
        return ed<right.ed;
    }
}

bool basicInfo::operator ==(basicInfo &right){
    return chr==right.chr&&st==right.st&&ed==right.ed&&sv==right.sv;
}

void basicInfo::printMe(FILE *fp) {
    if(fp==NULL){
        printf("chr:%d st:%lld ed:%lld\n",chr,st,ed);
    }
    else{
        fprintf(fp,"%d\t%lld\t%lld\t",chr,st,ed);
    }
}

void basicInfo::printRange(FILE *fp){
    if(fp==NULL){
        printf("range:%lld\n",ed-st);
    }
    else{
        fprintf(fp,"%lld\t",ed-st);
    }
}

long long basicInfo::getLength() {
    return ed-st+1;
}

char* basicInfo::getAgct(char *ref) {
    deleteCache();
    long long len=getLength();
    cache=new char[len+10];
    strncpy(cache,ref+st,len);
    cache[len]='\0';
    return cache;
}

void basicInfo::deleteCache() {
    if(cache!=NULL){
        delete []cache;
        cache=NULL;
    }
}

basicInfo::~basicInfo() {
    if(cache!=NULL){
        delete []cache;
    }
}

double basicInfo::overlapRate(const basicInfo &tmp) {
    if(chr!=tmp.chr||st>tmp.ed||ed<=tmp.st){
        return 0;
    }
    return ((double)min(ed,tmp.ed)-max(st,tmp.st))/(ed-st+1.0);
}

bool basicInfo::isIllegal() {
    if(chr<0||st<0||ed<0) {
        return true;
    }
    return false;
}

bool basicInfo::myStGreaterThan(const basicInfo &tmp){
    if(chr>tmp.chr){
        return true;
    }
    else if(chr==tmp.chr){
        return st>tmp.ed;
    }
    else{
        return false;
    }
}

bool basicInfo::myEdSmallerThan(const basicInfo &tmp){
    if(chr<tmp.chr){
        return true;
    }
    else if(chr=tmp.chr){
        return ed<tmp.st;
    }
    else{
        return false;
    }
}

bool basicInfo::myEdLargerThanEdOf(const basicInfo &tmp){
    if(chr>=tmp.chr&&ed>=tmp.ed){
        return true;
    }
    return false;
}

bool basicInfo::myEdSmallerThanStOf(const basicInfo &tmp){
    if(chr<=tmp.chr&&ed<=tmp.st){
        return true;
    }
    return false;
}
#endif
