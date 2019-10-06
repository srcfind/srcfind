//
// Created by ttbond on 19-6-10.
//

#ifndef SRCFIND_BAMFILE_H
#define SRCFIND_BAMFILE_H
#include <vector>
#include <string>
#include <map>
#include <set>
#include <htslib/sam.h>
#include "bedFile.h"
#include "basicInfo.h"
#include "alignInfo.h"
using namespace std;

class bamFile {
public:
    bamFile(const char * fileName);
    ~bamFile();
    vector<basicInfo> getReadsPos(bedFile &refRegion);
    vector<basicInfo> getRefPos(bam1_t *aln,const vector<basicInfo> &readBed,vector<bool> &matchNoSoftArea);
    vector<basicInfo> getRefPos(bam1_t *aln,const vector<basicInfo> &readBed);
    void nameMapStPos();
    vector<basicInfo> generateRefBedFromReadBed(map<string,vector<basicInfo> > &readBed,const char *bedFileName);
    void selectSA(const char *outBamFileName);
    set<string> getReadsSet();
    vector<alignInfo> simpleQuery(const char *qname=NULL);
    vector<alignInfo> simpleQueryAlign(const char *qname,long long readLen);
    void generateDotplot(const char *outP,const char *refFile,char *refCache,char *_qname);
    long long getRefMatchLen(bam1_t *aln);
    long long getReadLen(bam1_t *aln);
    void selectByRefRange(const char *inputBedName,const char *outPutBamName);
    void selectByRefRange(vector<basicInfo> inputBed,const char *outPutBamName);
    void selectByName(set<string> nameSet,const char *outPutBamName,bool selected=true);
    void selectByComAlign(const char *outPutBamName);
    void selectByRefOverlap(const char *outPutBamName);
    void selectByModel1(const char *outPutBamName);
    void sortByPos();
    bool isSorted();
    void modifyFileName(string &fn);
    void selectByReadLen(int minLen,const char *outPutBamName);
    void copyHeader(const char *tmpName,const char *outPutBamName);
    set<string> getQnameSet();
    vector<basicInfo> simpleQueryBasic(const char *qname);
    long long getHardClipLen(bam1_t *aln);
    char *fName;
    bool nameIsMapStPos;
};


#endif //SRCFIND_BAMFILE_H
