//
// Created by ttbond on 18-10-18.
//
#ifndef SRCFIND_GCSTA_CPP
#define SRCFIND_GCSTA_CPP
#include "gcSta.h"

gcSta::gcSta(){
    numa=numg=numc=numt=0;
}
void gcSta::addStr(char *str,long long len,long long num){
    char *ed=str+len;
    for(char *i=str;i<ed;i++){
        switch(*i){
            case 'a':
            case 'A':
                numa+=num;
                break;
            case 'g':
            case 'G':
                numg+=num;
                break;
            case 'c':
            case 'C':
                numc+=num;
                break;
            case 't':
            case 'T':
                numt+=num;
                break;
        }
    }
}
double gcSta::getGcContent(){
    return ((double)(numg+numc))/(numg+numc+numa+numt);
}
#endif