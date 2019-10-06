//
// Created by ttbond on 18-10-17.
//
#ifndef SRCFIND_AGCTTREE_CPP
#define SRCFIND_AGCTTREE_CPP
#include <stdio.h>
#include <stdlib.h>
#include "agctTree.h"
#include "agctNode.h"
agctTree::agctTree(){
    root=new agctNode();
    myGcSta=new gcSta();
    maxLevel=100;
    charCache=NULL;
}
void agctTree::release(agctNode *tar,int level){
    for(int i=0;i<4;i++){
        if(tar->son[i]!=NULL){
            release(tar->son[i],level+1);
        }
    }
    delete tar;
}
void agctTree::release(char *&tmp){
    if(tmp!=NULL){
        delete[] tmp;
    }
    tmp=NULL;
}
inline int agctTree::nu2int(char *from){
    switch(*from){
        case 'a':
        case 'A':
            return 0;
        case 'g':
        case 'G':
            return 1;
        case 'c':
        case 'C':
            return 2;
        case 't':
        case 'T':
            return 3;
        default:
            return -1;
    }
}
inline char agctTree::int2nu(int from){
    switch(from){
        case 0:
            return 'A';
        case 1:
            return 'G';
        case 2:
            return 'C';
        case 3:
            return 'T';
        default:
            printf("Wrong parameter in agctTree::int2nu value:%d\n",from);
            exit(-1);
    }
}
void agctTree::addSeq(agctNode *nowNode,char *from,int len,int level){
    if(level<maxLevel&&len>0){
        agctNode *nextNode=nowNode->enterSon(nu2int(from));
        addSeq(nextNode,from+1,len-1,level+1);
    }
    else{
        nowNode->num++;
    }
}
void agctTree::getTotalNum(std::vector<int> &rel,agctNode *nowNode,int level){
    if(nowNode->num>0&&level>=leastLen){
        rel.push_back(-(nowNode->num));
    }
    for(int i=0;i<4;i++){
        if(nowNode->son[i]!=NULL){
            getTotalNum(rel,nowNode->son[i],level+1);
        }
    }
}
void agctTree::printMostSeq(int rank,int _leastLen,FILE *fp){
    leastLen=_leastLen;
    std::vector<int> nums;
    getTotalNum(nums,root,0);
    sort(nums.begin(),nums.end());
    std::vector<int>::iterator it;
    for(it=nums.begin();it!=nums.end();it++){
        fprintf(fp,"%d\t",(*it));
    }
    minSelectFre=-nums[rank-1];
    release(charCache);
    charCache=new char[maxLevel+10];
    minSelectLen=leastLen;
    fprintf(fp,"\n%d\t%d\t\n",minSelectFre,minSelectLen);
    fpCache=fp;
    selectByFre(root,0);
    fprintf(fp,"%.5lf\n",myGcSta->getGcContent());
}
void agctTree::selectByFre(agctNode *nowNode,int level){
    for(int i=0;i<4;i++){
        if(nowNode->son[i]!=NULL){
            charCache[level]=int2nu(i);
            selectByFre(nowNode->son[i],level+1);
        }
    }
    if(nowNode->num>=minSelectFre&&level>=minSelectLen){
        charCache[level]='\0';
        myGcSta->addStr(charCache,strlen(charCache),nowNode->num);
        fprintf(fpCache,"%s\n",charCache);
    }
}
void agctTree::addSeq(int *lenSt,int lenEd,char *seq){
    for(int i=0;i<lenEd;i++){
        if(lenSt[i]<5){
            continue;
        }
        addSeq(root,seq+i,lenSt[i],0);
    }
}
agctTree::~agctTree() {
    delete myGcSta;
    release(root,0);
}

#endif