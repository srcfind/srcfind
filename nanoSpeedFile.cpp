//
// Created by ttbond on 19-6-11.
//

#include <string.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <set>
#include <string>
#include "nanoSpeedFile.h"
using namespace std;

char nanoSpeedFile::cache[1000000];
vector<long long> nanoSpeedFile::vecLLCache;

nanoSpeedFile::nanoSpeedFile(const char *fileName) {
    int fileNameLen=strlen(fileName);
    fName=new char[fileNameLen+1];
    strcpy(fName,fileName);
    indexName();
}



nanoSpeedFile::~nanoSpeedFile() {
    delete []fName;
}

void nanoSpeedFile::indexName() {
    readName2Id.clear();
    FILE *fp=fopen(fName,"r");
    int cot=0;
    long long pos=0;
    bool indexRel=true;
    fileStructure *tmpSt;
    string tmpReadName;

    while(fgets(cache,1000000,fp)!=NULL){
        int tmpLen=strlen(cache);
        switch(cot%4){
            case 0:
                if(cache[0]<'0'||cache[0]>'9'){
                    indexRel=false;
                }
                tmpSt=new fileStructure();
                tmpSt->speedPos=pos;
                break;
            case 1:
                if(cache[0]!='@'){
                    indexRel=false;
                }
                for(int i=strlen(cache)-1;i>=0&&
                (cache[i]==' '||cache[i]=='\n');i--){
                    cache[i]=0;
                }
                tmpSt->readNamePos=pos;
                tmpReadName=string(cache+1);
                break;
            case 2:
                if(cache[0]!='A'&&cache[0]!='C'&&cache[0]!='G'&&cache[0]!='T'){
                    indexRel=false;
                }
                tmpSt->sequencePos=pos;
                break;
            case 3:
                tmpSt->qualityPos=pos;
                readName2Id[tmpReadName]=*tmpSt;
                delete tmpSt;
        }
        pos+=tmpLen;
        cot++;
    }
    if(indexRel){
        printf("index %d reads with speed!\n",cot/4);
    }
    else{
        printf("%s index fault!\n",fName);
        exit(0);
    }
    fclose(fp);
}

void nanoSpeedFile::getSpeed(const char *readsName) {
    string tmpName(readsName);
    if(readName2Id.find(tmpName)==readName2Id.end()){
        printf("find %s failed!\n",readsName);
    }
    fileStructure pos=readName2Id[tmpName];
    getLine(pos.speedPos);
}

void nanoSpeedFile::getLine(long long pos) {
    FILE *fp=fopen(fName,"r");
    fseek(fp,pos,SEEK_SET);
    fgets(cache,1000000,fp);
    fclose(fp);
}

long long nanoSpeedFile::cutRegion(long long st,long long ed){
    long long pos=0,sepNum=0,rel=-1;
    for(char &c : cache){
        if(c==0){
            break;
        }
        if(c==' '){
            sepNum++;
        }
        if(rel==-1&&sepNum==st-1){
            rel=pos;
        }
        else if(sepNum==ed){
            c=0;
            return rel+1;
        }
        pos++;
    }
}

long long nanoSpeedFile::getReadLen(){
    if(cache==NULL){
        puts("Cache error while getting read length");
    }
    long long length=0;
    for(auto chr :cache){
        if(chr==0){
            break;
        }
        if(chr==' '){
            length++;
        }
    }
    return length;
}


void nanoSpeedFile::printOneSpeed(const basicInfo &it){
    getSpeed(it.cache);
    auto len=getReadLen();
    //left beside the regions
    //auto cacheSt=cutRegion(it.st-(it.ed-it.st),it.st);
    auto cacheSt=cutRegion(max(0,it.st-200),min(it.ed+200,len));
    //auto cacheSt=cutRegion(it.st,it.ed);
    puts(cache+cacheSt);
}

void nanoSpeedFile::generateSpeSpeedFile(const char *outFileName,vector<basicInfo> list){
    for(const basicInfo &it : list){
        printOneSpeed(it);
    }
}

vector<basicInfo> nanoSpeedFile::staMethod1(vector<long long> &data){
    const int lowerLimit=60;
    const int extendLen=20;
    const int siz=data.size();
    const int stPosOfReads=40;
    const int edPosOfReads=100;
    vector<basicInfo> rel;
    for(int i=stPosOfReads;i<siz-edPosOfReads;i++){
        if(data[i]>=lowerLimit){
            int st=max(1,i-extendLen+1);
            int ed=min(siz,i+extendLen+1);
            if(!rel.empty()&&rel.rbegin()->ed>=st){
                (*rel.rbegin()).ed=ed;
            }
            else{
                rel.push_back(basicInfo(0,st,ed));
            }
        }
    }
    return rel;
}

void nanoSpeedFile::strList2ll(const char *strCache){
    vecLLCache.clear();
    do {
        while (' ' == (*strCache)) {
            strCache++;
            if ('\n' == (*strCache) || 0 == (*strCache)) {
                return;
            }
        }
        long long tmpLL=-1;
        char tmpStr[1000];
        sscanf(strCache,"%lld",&tmpLL);
        if(tmpLL>-1){
            vecLLCache.push_back(tmpLL);
        }
        else{
            puts("error in nanoSpeedFile::strList2ll");
            exit(0);
        }
        sscanf(strCache,"%s",tmpStr);
        strCache+=strlen(tmpStr);
    }while('\n'!=(*strCache)&&0!=(*strCache));
}

void nanoSpeedFile::generateBedFile(vector<basicInfo> (nanoSpeedFile::*staMethod)(vector<long long> &)) {
    readName2readBed.clear();
    map<string,fileStructure>::iterator i;
    long long totalReadsNum=readName2Id.size(),cot=0;
    for(i=readName2Id.begin();i!=readName2Id.end();i++) {
        if(cot%100==0){
            printf("%d/%d now/total\n",cot,totalReadsNum);
        }
        getSpeed((i->first).c_str());
        strList2ll(cache);
        vector<basicInfo> staRel = (this->*staMethod)(vecLLCache);
        if (!staRel.empty()) {
            readName2readBed[i->first] = staRel;
        }
        cot++;
    }
    printf("total %d reads that have regions\n",readName2readBed.size());
}

void nanoSpeedFile::generateBedFile(vector<basicInfo> (nanoSpeedFile::*staMethod)(vector<long long int> &), set<string> &readsNameSet) {
    readName2readBed.clear();
    map<string,fileStructure>::iterator i;
    long long totalReadsNum=readName2Id.size(),cot=0;
    for(i=readName2Id.begin();i!=readName2Id.end();i++) {
        if(cot%100==0){
            printf("%d/%d now/total\n",cot,totalReadsNum);
        }
        cot++;
        if(readsNameSet.find(i->first)==readsNameSet.end()){
            continue;
        }
        getSpeed((i->first).c_str());
        strList2ll(cache);
        vector<basicInfo> staRel = (this->*staMethod)(vecLLCache);
        if (!staRel.empty()) {
            readName2readBed[i->first] = staRel;
        }
    }
    printf("total %d reads that have regions\n",readName2readBed.size());
}
