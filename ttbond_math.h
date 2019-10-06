//
// Created by ttbond on 18-9-14.
//

#ifndef SRCFIND_TTBOND_MATH_H
#define SRCFIND_TTBOND_MATH_H

#include<string.h>
#include<string>
using namespace std;

enum svType{DEL,DUP,INV,INS,OTHER};
enum nonBTyp{HAIRPIN,HDNA,DIRREP};

inline char *svType2str(char *rel,svType tmp){
    switch(tmp){
        case DEL:
            strcpy(rel,"DEL");break;
        case DUP:
            strcpy(rel,"DUP");break;
        case INV:
            strcpy(rel,"INV");break;
        case INS:
            strcpy(rel,"INS");break;
        case OTHER:
            strcpy(rel,"OTHER");break;
    }
    return rel;
}


//translate a small str data to int data
inline int stringChrName2int(const char *str){
    int rel=0;
    const char *i;
    if((*str)=='c'||(*str)=='C'){
        i=str+3;
    }
    else{
        i=str;
    }
    if((*i)=='x'||(*i)=='X'){
        return 23;
    }
    if((*i)=='y'||(*i)=='Y'){
        return 24;
    }
    if((*i)=='m'||(*i)=='M'){
        return 25;
    }
    for(;(*i)!='\0';i++){
        rel*=10;
        rel+=(*i)-'0';
    }
    if(rel>=1&&rel<=24){
        return rel;
    }
    else{
        return -1;
    }
}


inline long long string2ll(char *str){
    long long rel=0;
    char *i=str;
    for(;(*i)!='\0';i++){
        rel*=10;
        rel+=(*i)-'0';
    }
    return rel;
}

inline bool selectFromSepStr(char *rel,char *source,char *tar,char sepSymbol){
    int length=strlen(source);
    char cache[length+10];
    strcpy(cache,source);
    char *strSt=source,*i=source;
    for(;;i++){
        if((*i)==sepSymbol||(*i)=='\0'){
            bool endStrFlag=false;
            if((*i)=='\0'){
                endStrFlag=true;
            }
            (*i)='\0';
            bool findE=false;
            char *j=strSt;
            for(;j!=i;j++){
                if((*j)=='='){
                    findE=true;
                    (*j)='\0';
                    break;
                }
            }
            if(strcmp(strSt,tar)==0){
                if(findE){
                    strcpy(rel,j+1);
                }
                else{
                    rel[0]='\0';
                }
                return true;
            }
            if(endStrFlag){
                return false;
            }
            strSt=i+1;
        }
    }
}

inline char getComBase(char base){
    switch(base){
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'a':
            return 't';
        case 't':
            return 'a';
        case 'g':
            return 'c';
        case 'c':
            return 'g';
    }
}

inline int nuc2int(char nuc){
    switch(nuc){
        case 'A':
        case 'a':
            return 0;
        case 'G':
        case 'g':
            return 1;
        case 'C':
        case 'c':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return -1;
    }
}

inline char int2nuc(int num,bool capital=true){
    if(capital){
        switch(num){
            case 0:
                return 'A';
            case 1:
                return 'G';
            case 2:
                return 'C';
            case 3:
                return 'T';
            default:
                return 'N';
        }
    }
    else{
        switch(num){
            case 0:
                return 'a';
            case 1:
                return 'g';
            case 2:
                return 'c';
            case 3:
                return 't';
        }
    }
}

inline int max(int a,int b){
    if(a>b)return a;
    return b;
}

#endif //SRCFIND_TTBOND_MATH_H
