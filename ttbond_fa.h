//
// Created by ttbond on 19-6-10.
//

#ifndef SRCFIND_COPY_TTBOND_FA_H
#define SRCFIND_COPY_TTBOND_FA_H
#include <stdio.h>
#include <string.h>
#include <ttbond_math.h>
inline bool isRegularChr(char *str);
inline bool isNormalAgctSymbol(const char &tmp);
inline long long loadAgctByChr(int chrNum,char *fileName,char *agct);
inline bool isRegularChr(const char  *str){
    char *ref=">chr";
    char source[10];
    strncpy(source,str,4);
    source[4]='\0';
    if(strcmp(source,ref)==0){
        bool firstNum=false;
        const char *i;
        for(i=str+4;(*i)>='0'&&(*i)<='9';i++){
            firstNum=true;
        }
        if(firstNum==false&&(*(str+4)=='X'||*(str+4)=='Y'||*(str+4)=='x'||*(str+4)=='y')){
            firstNum=true;
            i=str+5;
        }
        if((*i)!=' '||firstNum==false){
            return false;
        }
        else{
            return true;
        }
    }
    else{
        //>1 not >chr1
        int len=strlen(str);
        const char *strEd=str+len;
        for(const char *i=str+1;i<strEd;i++){
            if((*i)==' '){
                break;
            }
            if((*i)<'0'&&(*i)>'9'&&(*i)!='X'&&(*i)!='x'&&(*i)!='Y'&&(*i)!='y'){
                return false;
            }
            if(i-str>2){
                return false;
            }
        }
        return true;
    }
}

//if the symbol is 'A/G/C/T/N' return true
inline bool isNormalAgctSymbol(const char &tmp){
    if(tmp=='A'||tmp=='G'||tmp=='C'||tmp=='T'||tmp=='N'){
        return true;
    }
    else{
        return false;
    }
}


//if the specific chrNum is found then it's sequencing will be loaded to array agct
// and return the length of the chr .... if not exist return -1
inline long long loadAgctByChr(int chrNum,const char *fileName,char *agct){
    printf("loading chr%d from %s...\n",chrNum,fileName) ;
    FILE *fp;
    fp=fopen(fileName,"r");
    char chrName[10];
    bool findChr=false;
    char tmp[1000];
    const char *ctmp=tmp;
    while(fgets(tmp,1000,fp)!=NULL){
        if(tmp[0]=='>'&&isRegularChr(ctmp)){
            if(tmp[1]=='c'||tmp[1]=='C') {
                sscanf(tmp + 4, "%s", chrName);
            }
            else{
                sscanf(tmp+1,"%s",chrName);
            }
            int chrN=stringChrName2int(chrName);
            if(chrN==chrNum){
                findChr=true;
                long long chrLen=0;
                char *tp=agct;
                bool comZero=false;
                while(fgets(tp,1000,fp)!=NULL){
                    if(!isNormalAgctSymbol(*tp)){
                        *tp='\0';
                        comZero=true;
                        break;
                    }
                    int newLen=strlen(tp);
                    tp+=newLen-1;
                    chrLen+=newLen-1;
                }
                if(!comZero){
                    *tp='\0';
                }
                printf("loading successfully %lld agct\n",chrLen);
                return chrLen;
            }
        }
    }
    if(!findChr){
        return -1;
    }
}

#endif //SRCFIND_COPY_TTBOND_FA_H
