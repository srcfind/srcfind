//
// Created by ttbond on 19-8-8.
//

#ifndef SRCFIND_ONEFASTQ_H
#define SRCFIND_ONEFASTQ_H
#include<stdio.h>
#include<string>
using namespace std;

class oneFastq {
public:
    oneFastq();
    oneFastq(const oneFastq &tmp);
    oneFastq& operator =(const oneFastq &tmp);
    bool getNext(FILE * &fp,bool storeInfo=false,bool storeSeq=false,bool storeQual=false);
    void write(FILE * &fp);
    bool checkCacheContent(char *tmpCache);
    void deleteSelf(char *add);
    ~oneFastq();
    char *seqCache,*qualCache,*infoCache;
    string qname;
    long long len;
};


#endif //SRCFIND_ONEFASTQ_H
