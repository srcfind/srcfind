//
// Created by ttbond on 18-9-16.
//

#ifndef SRCFIND_DETECTREGION_H
#define SRCFIND_DETECTREGION_H
#include<stdio.h>
#include<string.h>
#include<vector>
#include<cstring>
#include<stdlib.h>
#include"ttbond_math.h"
#include "basicInfo.h"


class detectRegion {
public:
    detectRegion(char *fileName);
    detectRegion(char *str,const basicInfo &source);
    detectRegion(char *str,int _chr,long long _st,long long _ed,svType _mysv);
    detectRegion(char *str1,int _chr,long long _st1,long long _ed1,char *str2,int _chr2,long long _st2,long long _ed2,svType _mysv);
    ~detectRegion();
    double *getReverseComScore(FILE *fp=NULL);
    double *getDirectRepeatScore(FILE *fp=NULL);
    double *getMirrorRepeatScore(FILE *fp=NULL);
    int printMaxScore(int *score,int *scoreEd,FILE *fp);
    std::vector< std::vector<basicInfo> > getPotDirPos();
    void localizeFirstBasePos();
    void releaseFirstBasePos(int **&_firstBasePos);
    void release(double *&toRelease);
    void release(char *&toRelease);
    void release(int **&toRelease);
    void releaseDetectRel();
    void release(int *&toRelease);
    void revComDfs(int *st,int *ed,int *misMatchSt,int *misMatchEd,int *st2,int *ed2,int *misMatchSt2,int *misMatchEd2,int depth);
    void mirRepDfs(int *st,int *ed,int *misMatchSt,int *misMatchEd,int *st2,int *ed2,int *misMatchSt2,int *misMatchEd2,int depth);
    void dirRepDfs(int *st,int *ed,int *misMatchSt,int *misMatchEd,int depth);
    void dirRepDfsRec(int *st,int *ed,int *misMatchSt,int *misMatchEd,int depth);
    int *getInitMisSt(int iteNum,int **_firstBasePos,int *_cacheForDfs);
    int *getInitMisEd(int iteNum,int **_firstBasePos,int *_cacheForDfs);
    void copyFirstBasePos();
    void printDetectRel(FILE *fp=NULL,bool test=false);
    void printScore(double *score,double *scoreEd,FILE *fp);
    double *statisticScore(double *score,double &sumScore,FILE *fp=NULL);
    bool operator <(detectRegion &right);

    static bool firstRegion;
    static bool firstLineOfFile;
    static char *base;
    static int maxMismatchNum;
    static int maxDetectLen;
    static int minDetectLen;
    int chr,st,ed,length;
    svType mySv;
    //if the region is constructed by two regions in genome
    //this will be used
    int chr2,st2,ed2;
    char *agct,*agctEnd;
    double *reverseComScore,sumRevComScore;
    double *directRepScore,sumDirRepScore;
    double *mirrorRepScore,sumMirRepScore;
    int **firstBasePos,**firstBasePos2;
    int *indEndOfFirstBasePos[4],*indEndOfFirstBasePos2[4];
    int *cacheForDfs,*cacheForDfs2;
    long long cacheLen;
    int *detectRel;
    std::vector< std::vector<basicInfo> > recRel;
};


#endif //SRCFIND_DETECTREGION_H
