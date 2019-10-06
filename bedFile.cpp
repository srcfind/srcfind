//
// Created by ttbond on 18-10-11.
//
#ifndef SRCFIND_BEDFILE_CPP
#define SRCFIND_BEDFILE_CPP
#include "bedFile.h"
using namespace std;


//to initialize a void bedFile object
bedFile::bedFile() {
    agct=NULL;
    nowPos=0;
    nowRefChr=0;
}

//read bed file and store basicInfo into myRegion
bedFile::bedFile(const char *fileName,char *_agct) {
    agct=_agct;
    nowPos=0;
    nowRefChr=0;
    FILE *fp=fopen(fileName,"r");
    char cache[100];
    int chr;
    long long st,ed;
    while(fgets(cache,100,fp)!=NULL){
        char *cacheSt=cache;
        if((*cache)=='c'||(*cache)=='C'){
            cacheSt=cache+3;
        }
        sscanf(cacheSt,"%d%lld%lld",&chr,&st,&ed);
        myRegion.push_back(basicInfo(chr,st,ed));
    }
    fclose(fp);
}

bedFile::bedFile(std::vector <basicInfo> _source,char *_agct) {
    myRegion=_source;
    agct=_agct;
    nowPos=0;
    nowRefChr=0;
}

bool bedFile::load(const char *fileName,fileT fileType){
    FILE *fp=fopen(fileName,"r");
    if(fp==NULL){
        printf("bedFile::load file error\n");
        exit(0);
    }
    bool loadStatus=false;
    switch(fileType){
        case BEDF: break;
        case NONBF:loadStatus=loadNonBF(fp);break;
    }
}

bool bedFile::loadNonBF(FILE * fp) {
    const int lineCacheSize=100000,partCacheSize=20000;
    char lineCache[lineCacheSize],partCache[partCacheSize],chr[10];
    int chrNum;
    while(fgets(lineCache,lineCacheSize,fp)!=NULL){
        long long st,ed;
        if(lineCache[0]=='#' || lineCache[0]=='\n'){
            continue;
        }
        sscanf(lineCache,"%s%s%s%lld%lld",partCache,partCache,partCache,&st,&ed);
        sscanf(lineCache,"%s",partCache);

        if(partCache[0]=='c'||partCache[0]=='C'){
            sscanf(partCache+3,"%s",chr);
            if(chr[0]=='X'||chr[0]=='x'){
                chrNum=23;
            }
            else if(chr[0]=='Y'||chr[0]=='y'){
                chrNum=24;
            }
            else{
                sscanf(chr,"%d",&chrNum);
            }
        }
        myRegion.push_back(basicInfo(chrNum,st,ed,OTHER));
    }

}

void bedFile::solAndPrint(char *fileName,char *writeType) {


    std::vector<basicInfo>::iterator it;
    FILE *fp=fopen(fileName,writeType);
    int refChr=0;
    //for debug
    /*
    for(it=myRegion.begin();it!=myRegion.end();it++){
        if(refChr!=(*it).chr){
            refChr=(*it).chr;
            loadAgctByChr(refChr,"GRCh38.d1.vd1.fa",agct);
        }
        detectRegion region(agct,(*it).chr,(*it).st,(*it).ed,(*it).sv);
        region.getReverseComScore(fp);
        region.getDirectRepeatScore(fp);
        region.getMirrorRepeatScore(fp);
    }
     */

    //for background set
    int infoNum=myRegion.size();
    for(int i=0;i<infoNum;i++){
        if(myRegion[i].chr!=refChr){
            refChr=myRegion[i].chr;
            loadAgctByChr(refChr,"GRCh38.d1.vd1.fa",agct);
        }
        detectRegion region(agct,myRegion[i].chr,myRegion[i].st,myRegion[i].ed,myRegion[i].sv);
        region.getMirrorRepeatScore();
        region.getDirectRepeatScore();
        region.getReverseComScore();
        myRel.push_back(solRel(region));
    }
    for(int i=0;i<infoNum;i++){
        if(myRegion[i].chr==24){
            continue;
        }
        myRel[i].printMe(fp);
        fprintf(fp,"%lld\n",myRel[i].info->getLength());
    }
    fclose(fp);
}
void bedFile::printSingleRcdScore(char *fileName){
    char fileNameCache[50];
    strcpy(fileNameCache,fileName);
    FILE *mirFp=fopen(strcat(fileNameCache,"Mir.scr"),"w");
    strcpy(fileNameCache,fileName);
    FILE *revFp=fopen(strcat(fileNameCache,"Rev.scr"),"w");
    strcpy(fileNameCache,fileName);
    FILE *dirFp=fopen(strcat(fileNameCache,"Dir.scr"),"w");
    int infoNum=myRegion.size();
    int refChr=0;
    for(int i=0;i<infoNum;i++){
        if(myRegion[i].chr!=refChr){
            refChr=myRegion[i].chr;
            loadAgctByChr(refChr,"GRCh38.d1.vd1.fa",agct);
        }
        detectRegion region(agct,myRegion[i].chr,myRegion[i].st,myRegion[i].ed,myRegion[i].sv);
        region.getMirrorRepeatScore();
        region.printScore(region.mirrorRepScore,region.mirrorRepScore+region.length,mirFp);
        region.getReverseComScore();
        region.printScore(region.reverseComScore,region.reverseComScore+region.length,revFp);
        region.getDirectRepeatScore();
        region.printScore(region.directRepScore,region.directRepScore+region.length,dirFp);
    }
    fclose(mirFp);
    fclose(revFp);
    fclose(dirFp);
}

char *bedFile::getChrAgct(basicInfo &tar){
    if(nowRefChr!=tar.chr){
        nowRefChr=tar.chr;
        loadAgctByChr(nowRefChr,"GRCh38.d1.vd1.fa",agct);
    }
    return agct;
}

basicInfo bedFile::getNowRegion() {
    if(atRegionEnd()){
        puts("At illeagle index");
        exit(0);
    }
    return myRegion[nowPos++];
}

void bedFile::init(){
    nowPos=0;
    nowRefChr=0;
}

bool bedFile::atRegionEnd() {
    if(nowPos>=(int)myRegion.size()){
        return true;
    }
    else{
        return false;
    }
}

void bedFile::sortMyRegion(){
    sort(myRegion.begin(),myRegion.end());
}
#endif
