//
// Created by ttbond on 19-7-4.
//


#include <string.h>
#include <stdio.h>
#include <math.h>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include "fastqFile.h"
using namespace std;

fastqFile::fastqFile(const char *sourceFile) {
    int len=strlen(sourceFile);
    fileName=new char[len+5];
    strcpy(fileName,sourceFile);
}

fastqFile::~fastqFile() {
    delete []fileName;
}

void fastqFile::generateSubFile(const char *newFileName, const char *selectRangeFileName) {
    FILE *ofp=fopen(newFileName,"a");
    FILE *ifp=fopen(fileName,"r");
    FILE *rangeFp=fopen(selectRangeFileName,"r");
    set<string> range;
    char lineCache[1000000];
    int lineCacheSize=1000000;
    while(fgets(lineCache,lineCacheSize,rangeFp)!=NULL){
        int len=strlen(lineCache);
        lineCache[len-1]=0;
        range.insert(string(lineCache));
    }
    fclose(rangeFp);
    while(fgets(lineCache,lineCacheSize,ifp)!=NULL) {
        if('@'==lineCache[0]&&'-'==lineCache[9]&&'-'==lineCache[14]&&'-'==lineCache[19]){
            int len=strlen(lineCache);
            for(int i=0;i<len;i++){
                if(lineCache[i]=='_'){
                    lineCache[i]=0;
                    break;
                }
            }
            if(range.find(string(lineCache+1))!=range.end()){
                lineCache[strlen(lineCache)]='_';
                fprintf(ofp,"%s",lineCache);
                char tmpCharCache;
                int enterNum=0;
                while(tmpCharCache=getc(ifp)){
                    if('\n'==tmpCharCache){
                        enterNum++;
                    }
                    fprintf(ofp,"%c",tmpCharCache);
                    if(3==enterNum){
                        break;
                    }
                }
            }
        }
    }

}

void fastqFile::classifyFailAndPass1D() {
    ifstream sourceF;
    ofstream passOutF,failOutF;
    sourceF.open(fileName);
    passOutF.open(string(fileName)+".pass.fastq");
    failOutF.open(string(fileName)+".fail.fastq");
    string qname,qname2,seq,plus,qua;
    int cot=0;
    while(sourceF>>qname>>qname2>>seq>>plus>>qua){
        cot++;
        if(qname[0]!='@'){
            cout << "cot:" << cot << endl;
            cout << qname << endl;
            puts("error in classifyFailAndPass1D()");
            exit(0);
        }
        double totalErp=0;
        for(int i:qua){
            totalErp+=pow(10,-(i-33)/10.0);
        }
        totalErp/=qua.length();
        double meanQ=-10*log10(totalErp);
        if(meanQ>=7){
            passOutF << qname << " " << qname2 << endl << seq << endl << plus << endl << qua << endl;
        }
        else{
            failOutF << qname << " " << qname2 << endl << seq << endl << plus << endl << qua << endl;
        }
    }
    sourceF.close();
    passOutF.close();
    failOutF.close();
}

void fastqFile::getFast5FileName(const char *newFileName) {

}

void fastqFile::selectByName(set<string> &nameSet,const char *outF){
    char lineCache[8000005];
    int lineCacheSize=8000000;
    FILE *fp=fopen(fileName,"r");
    FILE *ofp=fopen(outF,"w");
    long long nowLine=0,lineFlag=-100;
    int num=0;
    while(fgets(lineCache,lineCacheSize,fp)!=NULL) {
        nowLine++;
        long long len=strlen(lineCache);
        if(len>=lineCacheSize-1){
            puts("not enough memory for one line!");
            exit(-1);
        }
        lineCache[len-1]=0;
        if('@'==lineCache[0]&&'-'==lineCache[9]&&'-'==lineCache[14]&&'-'==lineCache[19]){
            string nowName;
            for(char *i=lineCache+1;(*i)!=' ';i++){
                nowName+=(*i);
            }
            if(nameSet.find(nowName)!=nameSet.end()){
                lineFlag=nowLine;
            }
        }
        if(nowLine<lineFlag+4){
            num++;
            printf("num:%d\n",num);
            fprintf(ofp,"%s\n",lineCache);
        }
    }
    fclose(fp);
    fclose(ofp);
}

vector<long long> fastqFile::staReadsLen(){
    char lineCache[1000005];
    int lineCacheSize=1000000;
    FILE *fp=fopen(fileName,"r");
    long long nowLine=0;
    vector<long long>rel;
    while(fgets(lineCache,lineCacheSize,fp)!=NULL) {
        nowLine++;
        if(nowLine%100==0){
            cout << nowLine << endl;
        }
        if(nowLine%4==2){
            rel.push_back(strlen(lineCache)-1);
        }
    }
    return rel;
}

set<string> fastqFile::getQnameSet(){
    char lineCache[8000005];
    int lineCacheSize=8000000;
    FILE *fp=fopen(fileName,"r");
    long long nowLine=0;
    set<string>rel;
    while(fgets(lineCache,lineCacheSize,fp)!=NULL) {
        long long len=strlen(lineCache);
        if(len>=lineCacheSize-1){
            puts("not enough memory for one line!");
            exit(-1);
        }
        lineCache[len-1]=0;
        if('@'==lineCache[0]&&'-'==lineCache[9]&&'-'==lineCache[14]&&'-'==lineCache[19]) {
            string nowName;
            for(char *i=lineCache+1;(*i)!=' ';i++){
                nowName+=(*i);
            }
            rel.insert(nowName);
        }
    }
    return rel;
}

void fastqFile::selectByLen(const long long lowLim, const char *outF) {
    FILE *ifp=fopen(fileName,"r");
    FILE *ofp=fopen(outF,"w");
    oneFastq nowFastq;
    while(nowFastq.getNext(ifp,true,true,true)){
        if(nowFastq.len>=1000){
            nowFastq.write(ofp);
        }
    }
}

vector<oneFastq> fastqFile::selectByName(set<string> &nameSet) {
    vector<oneFastq> rel;
    FILE *ifp = fopen(fileName, "r");
    oneFastq nowFastq;
    int num=0,totalNum=0;
    while (nowFastq.getNext(ifp, true, true, true)){
        totalNum++;
        if(totalNum%100==0){
            printf("total:%d\n",totalNum);
        }
        if(nameSet.find(nowFastq.qname)!=nameSet.end()){
            rel.push_back(nowFastq);
            num++;
            if(num%100==0){
                printf("%d\n",num);
            }
        }
    }
    return rel;
}

vector<oneFastq> fastqFile::selectByNameL(set<string> &nameSet) {
    vector<oneFastq> rel;
    FILE *ifp = fopen(fileName, "r");
    oneFastq nowFastq;
    int num=0,totalNum=0;
    while (nowFastq.getNext(ifp, false, false, false)){
        totalNum++;
        if(totalNum%100==0){
            printf("total:%d\n",totalNum);
        }
        if(nameSet.find(nowFastq.qname)!=nameSet.end()){
            rel.push_back(nowFastq);
            num++;
            if(num%100==0){
                printf("%d\n",num);
            }
        }
    }
    return rel;
}

