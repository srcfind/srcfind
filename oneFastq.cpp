//
// Created by ttbond on 19-8-8.
//

#include <string.h>
#include "oneFastq.h"

oneFastq::oneFastq() {
    seqCache=NULL;
    qualCache=NULL;
    infoCache=NULL;
}

bool oneFastq::getNext(FILE *&fp,bool storeInfo,bool storeSeq,bool storeQual) {
    char lineCache[8000005];
    int lineCacheSize=8000000;
    if(fgets(lineCache,lineCacheSize,fp)!=NULL){
        checkCacheContent(lineCache);
        qname=string("");
        for(char *i=lineCache+1;(*i)!=' ';i++){
            qname+=*i;
        }
        if(seqCache!=NULL){
            delete []seqCache;
        }
        if(qualCache!=NULL){
            delete []qualCache;
        }
        if(infoCache!=NULL){
            delete []infoCache;
        }
        if(storeInfo){
            long long infoLen=strlen(lineCache);
            lineCache[infoLen-1]=0;
            infoCache=new char[infoLen+5];
            strcpy(infoCache,lineCache);
        }
        fgets(lineCache,lineCacheSize,fp);
        checkCacheContent(lineCache);
        len=strlen(lineCache)-1;
        if(storeSeq){
            lineCache[len]=0;
            seqCache=new char[len+5];
            strcpy(seqCache,lineCache);
        }
        fgets(lineCache,lineCacheSize,fp);
        checkCacheContent(lineCache);
        fgets(lineCache,lineCacheSize,fp);
        checkCacheContent(lineCache);
        if(storeQual){
            lineCache[len]=0;
            qualCache=new char[len+5];
            strcpy(qualCache,lineCache);
        }
        return true;
    }
    else{
        return false;
    }
}

bool oneFastq::checkCacheContent(char *tmpCache){
    long long len=strlen(tmpCache);
    if(len>=8000000-1){
        puts("not enough memory for one line!");
        exit(-1);
    }
}

void oneFastq::write(FILE *&fp) {
    fprintf(fp,"%s\n",infoCache);
    fprintf(fp,"%s\n",seqCache);
    fprintf(fp,"+\n");
    fprintf(fp,"%s\n",qualCache);
}

oneFastq::~oneFastq() {
    if(infoCache==NULL){
        delete []infoCache;
    }
    if(seqCache==NULL){
        delete []seqCache;
    }
    if(qualCache==NULL){
        delete []qualCache;
    }
}

oneFastq::oneFastq(const oneFastq &tmp) {
    qname=tmp.qname;
    len=tmp.len;
    if(tmp.seqCache!=NULL){
        long long len=strlen(tmp.seqCache);
        seqCache=new char[len+10];
        strcpy(seqCache,tmp.seqCache);
    }
    else{
        seqCache=NULL;
    }
    if(tmp.qualCache!=NULL){
        long long len=strlen(tmp.qualCache);
        qualCache=new char[len+10];
        strcpy(qualCache,tmp.qualCache);
    }
    else{
        qualCache=NULL;
    }
    if(tmp.infoCache!=NULL){
        long long len=strlen(tmp.infoCache);
        infoCache=new char[len+10];
        strcpy(infoCache,tmp.infoCache);
    }
    else{
        infoCache=NULL;
    }
}

oneFastq &oneFastq::operator=(const oneFastq &tmp) {
    qname=tmp.qname;
    len=tmp.len;
    if(tmp.seqCache!=NULL){
        long long len=strlen(tmp.seqCache);
        seqCache=new char[len+10];
        strcpy(seqCache,tmp.seqCache);
    }
    else{
        seqCache=NULL;
    }
    if(tmp.qualCache!=NULL){
        long long len=strlen(tmp.qualCache);
        qualCache=new char[len+10];
        strcpy(qualCache,tmp.qualCache);
    }
    else{
        qualCache=NULL;
    }
    if(tmp.infoCache!=NULL){
        long long len=strlen(tmp.infoCache);
        infoCache=new char[len+10];
        strcpy(infoCache,tmp.infoCache);
    }
    else{
        infoCache=NULL;
    }
    return *this;
}

void oneFastq::deleteSelf(char *add) {
    if(add!=NULL){
        delete []add;
    }
}

