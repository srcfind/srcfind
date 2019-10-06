//
// Created by ttbond on 19-7-4.
//

#ifndef SRCFIND_FASTQFILE_H
#define SRCFIND_FASTQFILE_H
#include <string>
#include <iostream>
#include <set>
#include <vector>
#include "oneFastq.h"
using namespace std;

class fastqFile {
public:
    fastqFile(const char *sourceFile);
    ~fastqFile();
    void generateSubFile(const char *newFileName,const char *selectRangeFileName);
    void classifyFailAndPass1D();
    void getFast5FileName(const char *newFileName);
    void selectByName(set<string> &nameSet,const char *outF);
    vector<oneFastq> selectByNameL(set<string> &nameSet);
    vector<oneFastq> selectByName(set<string> &nameSet);
    vector<long long> staReadsLen();
    void selectByLen(const long long lowLim,const char *outF);
    set<string> getQnameSet();
    char *fileName;


    class oneFaqstq{
        char *name;

    };
};


#endif //SRCFIND_FASTQFILE_H
