//
// Created by ttbond on 18-10-15.
//
#ifndef SRCFIND_FAFILE_CPP
#define SRCFIND_FAFILE_CPP
#include "faFile.h"
faFile::faFile(char *_agct){
    agct=_agct;
}

inline bool faFile::isRegularChr(const char *str){
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
        if(firstNum==false&&(*(str+4)=='X'||*(str+4)=='Y'||*(str+4)=='x'||*(str+4)=='y'||*(str+4)=='m'||*(str+4)=='M')){
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
            if((*i)<'0'&&(*i)>'9'&&(*i)!='X'&&(*i)!='x'&&(*i)!='Y'&&(*i)!='y'&&(*i)!='M'&&(*i)!='m'){
                return false;
            }
            if(i-str>2){
                return false;
            }
        }
        return true;
    }
}

inline bool faFile::isNormalAgctSymbol(const char &tmp){
    if(tmp=='A'||tmp=='G'||tmp=='C'||tmp=='T'||tmp=='N'){
        return true;
    }
    else{
        return false;
    }
}

long long faFile::loadAgctByChr(int chrNum,char *fileName){
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

void faFile::getValidatedRegionByChr(int chrNum,char *fileName){
    validatedRegion.clear();
    agctLen=loadAgctByChr(chrNum,fileName);
    if(agctLen==-1){
        printf("loading chr:%d failed",chrNum);
        exit(-1);
    }
    char *agctEd=agct+agctLen;
    bool noN=false;
    long long tmpSt,tmpEd;
    for(char *i=agct;i<agctEd;i++){
        if((*i)!='N'&&noN==false){
            noN=true;
            tmpSt=i-agct;
        }
        if((*i)=='N'&&noN==true){
            noN=false;
            tmpEd=i-agct-1;
            validatedRegion.push_back(basicInfo(chrNum,tmpSt,tmpEd));
        }
    }
}


void faFile::getRandomBackgroundRegion(long long length,double percent,std::vector<basicInfo> &dest){
    std::vector<long long> validatedNum;
    std::vector<basicInfo>::iterator it;
    std::vector<long long>::iterator lit;
    totalValidatedNum=0;
    for(it=validatedRegion.begin();it!=validatedRegion.end();it++){
        long long tmpNum=(*it).getLength()-length+1;
        validatedNum.push_back(tmpNum);
        totalValidatedNum+=tmpNum;
    }
    long long requestNum=totalValidatedNum*percent;
    std::set<long long> nowSet;
    std::default_random_engine random(time(NULL));
    std::uniform_int_distribution<long long> dis(0,totalValidatedNum);
    long long nowPos;
    printf("Genearating %lld sub backgroud regions...\n",requestNum);
    for(int i=0;i<requestNum;i++){
        nowPos=dis(random);
        while(nowSet.find(nowPos)!=nowSet.end()){
            nowPos=dis(random);
        }
        long long tmpSum=0,tmpId=0;
        for(lit=validatedNum.begin();lit!=validatedNum.end();lit++,tmpId++){
            long long kk=(*lit);
            if(tmpSum+kk>=nowPos){
                basicInfo tmpInfo=validatedRegion[tmpId];
                long long newSt=nowPos-tmpSum-1+tmpInfo.st;
                long long newEd=newSt+length-1;
                dest.push_back(basicInfo(tmpInfo.chr,newSt,newEd));
                break;
            }
            tmpSum+=kk;
        }
    }
}

void faFile::getRandomPos(int posNum, std::vector<basicInfo> &dest) {
    std::vector<long long> validatedNum;
    std::vector<basicInfo>::iterator it;
    std::vector<long long>::iterator lit;
    totalValidatedNum=0;
    for(it=validatedRegion.begin();it!=validatedRegion.end();it++){
        long long tmpNum=(*it).getLength();
        validatedNum.push_back(tmpNum);
        totalValidatedNum+=tmpNum;
    }
    long long requestNum=posNum;
    std::default_random_engine random(time(NULL));
    std::uniform_int_distribution<long long> dis(0,totalValidatedNum);
    long long nowPos;
    printf("Genearating %lld randomPos...\n",requestNum);
    for(int i=0;i<requestNum;i++){
        nowPos=dis(random);
        long long tmpSum=0,tmpId=0;
        for(lit=validatedNum.begin();lit!=validatedNum.end();lit++,tmpId++){
            long long kk=(*lit);
            if(tmpSum+kk>=nowPos){
                basicInfo tmpInfo=validatedRegion[tmpId];
                long long newSt=nowPos-tmpSum-1+tmpInfo.st;
                long long newEd=-1;
                dest.push_back(basicInfo(tmpInfo.chr,newSt,newEd));
                break;
            }
            tmpSum+=kk;
        }
    }
}

#endif