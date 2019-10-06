//
// Created by ttbond on 19-6-11.
//

#ifndef SRCFIND_COPY_NANOSPEEDFILE_H
#define SRCFIND_COPY_NANOSPEEDFILE_H
#include <iostream>
#include <map>
#include <vector>
#include <string.h>
#include <set>
#include "basicInfo.h"
using namespace std;



class nanoSpeedFile {
public:
    class fileStructure{
    public:
        fileStructure(){
            memset(this,0,sizeof(fileStructure));
        }
        long long speedPos,readNamePos,sequencePos,qualityPos;
    };
    nanoSpeedFile(const char * fileName);
    ~nanoSpeedFile();
    void indexName();
    void getSpeed(const char *readsName);
    void getLine(long long pos);
    void generateSpeSpeedFile(const char *outFileName,vector<basicInfo> list);
    long long cutRegion(long long st,long long ed);
    void printOneSpeed(const basicInfo &it);
    long long getReadLen();
    void strList2ll(const char *strCache);
    vector<basicInfo> staMethod1(vector<long long>&);
    void generateBedFile(vector<basicInfo> (nanoSpeedFile::*staMethod)(vector<long long>&));
    void generateBedFile(vector<basicInfo> (nanoSpeedFile::*staMethod)(vector<long long>&),set<string> &readsNameSet);
    static char cache[1000000];
    static vector<long long> vecLLCache;
    char *fName;
    map<string,fileStructure> readName2Id;
    map<string,vector<basicInfo> >readName2readBed;
};


#endif //SRCFIND_COPY_NANOSPEEDFILE_H
