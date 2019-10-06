//
// Created by ttbond on 19-8-21.
//

#include "alignManager.h"


RInside alignManager::R;

int alignManager::code(const char tmp,bool com) {
    if(com){
        switch(tmp){
            case 'A':return 4;
            case 'G':return 2;
            case 'C':return 1;
            case 'T':return 0;
        }
    }
    switch(tmp){
        case 'A':return 0;
        case 'G':return 1;
        case 'C':return 2;
        case 'T':return 4;
    }
    return 0;
}

void alignManager::generateKmer(map<long long, vector<long long> > &rel, const char *seq,bool com) {
    long long siz=strlen(seq);
    int kmerSize=10;
    long long unit=1;
    long long cut=(unit<<(kmerSize*2))-1;
    long long nowk=0;
    if(!com) {
        for (long long i = 0; i < siz; i++) {
            nowk <<= 2;
            nowk &= cut;
            nowk |= code(seq[i]);
            if (i >= kmerSize - 1) {
                rel[nowk].push_back(i);
            }
        }
    }
    else{
        for(long long i=0;i<siz;i++) {
            nowk <<= 2;
            nowk &= cut;
            nowk |= code(seq[siz - 1 - i], true);
            if (i >= kmerSize - 1) {
                rel[nowk].push_back(siz-1-i);
            }
        }
    }
}

vector<pair<int, int> > alignManager::getXY(bool refCom,bool readCom) {
    map1.clear();
    map2.clear();
    vector<pair <int,int> >rel;
    generateKmer(map1,seq1,refCom);
    generateKmer(map2,seq2,readCom);
    for(auto it=map1.begin();it!=map1.end();it++){
        long long kmer=it->first;
        if(map2.find(kmer)!=map2.end()){
            vector<long long>&vec1=it->second,&vec2=map2[kmer];
            for(auto i : vec1){
                for(auto j:vec2){
                    rel.push_back(make_pair(i,j));
                }
            }
        }
    }
    return rel;
}
/*
//modelB 的函数主体
alignManager::alignManager(oneFastq &read, char *refFile, const char *bamFileName,char *refCache,const char *outF,bool plot) {
    if(plot){
        R.parseEvalQ("library('ggplot2')");
    }
    cache=refCache;
    posId=NULL;
    idNum=NULL;
    posId2=NULL;
    idNum2=NULL;
    secPosId=NULL;
    secIdNum=NULL;
    segTreeL=segTreeR=segNodesHave=segNodesFull=NULL;
    memset(segPosId,0,sizeof(segPosId));
    memset(segIdNum,0,sizeof(segIdNum));
    //截取ref，并准备好readSeq
    bamFile bamF(bamFileName);
    vector<alignInfo> posReg;
    posReg=bamF.simpleQuery(read.qname.c_str());
    long long refMinPos=1e18,refMaxPos=-1;
    int chr=-1;
    map<int,vector<alignInfo> > pos;
    for(auto info:posReg){
        if(info.refInfo.chr<0||info.refInfo.chr>24){
            continue;
        }
        pos[info.refInfo.chr].push_back(info);
    }
    for(auto info:pos){
        refMinPos=1e18,refMaxPos=-1;
        int siz=info.second.size();
        bool refOverlap=false;
        long long readDis[2];
        memset(readDis,0,sizeof(readDis));
        int longerDir;
        for(int i=0;i<siz;i++){
            refMinPos=min(refMinPos,info.second[i].refInfo.st);
            refMaxPos=max(refMaxPos,info.second[i].refInfo.ed);
            readDis[(bool)(info.second[i].flag&16)]+=info.second[i].readInfo.ed-info.second[i].readInfo.st;
            for(int j=i+1;j<siz;j++){
                if((info.second)[i].refInfo.overlapRate((info.second)[j].refInfo)>0 && (info.second[i].flag&16)!=(info.second[j].flag&16)){
                    refOverlap=true;
                    break;
                }
            }
        }
        if(readDis[0]>=readDis[1]){
            longerDir=0;
        }
        else{
            longerDir=1;
        }
        printf("longerDir:%d\n",longerDir);
        bool oneInvPass=false;
        for(int j=0;j<siz;j++){
            if(((bool)(info.second[j].flag&16))!=longerDir){
                bool thisInvPass=true;
                for(int i=0;i<siz;i++){
                    if(((bool)(info.second[i].flag&16))==longerDir) {

                        basicInfo refPos=getRefPosFromAlign(info.second[i],info.second[j],longerDir);
                        refPos.printMe();
                        faFile ref(refCache);
                        ref.loadAgctByChr(refPos.chr, refFile);
                        char *refReg1,*refReg2;
                        refReg1=new char[refPos.ed-refPos.st+50];
                        refReg2=new char[info.second[j].refInfo.ed-info.second[j].refInfo.st+50];
                        char *refPint=refReg1;
                        for (char *it = refCache + refPos.st; it <= refCache + refPos.ed; it++, refPint++) {
                            *refPint = *it;
                        }
                        *refPint=0;
                        refPint=refReg2;
                        for (char *it = refCache + info.second[j].refInfo.st; it <= refCache + info.second[j].refInfo.ed; it++, refPint++) {
                            *refPint = *it;
                        }
                        *refPint=0;
                        seq1=refReg1;
                        seq2=refReg2;
                        auto rel=getXY(true,false);
                        sort(rel.begin(),rel.end());
                        getCluster(rel);
                        long long refLen=refPos.ed-refPos.st;
                        auto negKLen=getClusterRefLen(posId,idNum,rel,refLen);
                        auto posKLen=getClusterRefLen(posId2,idNum2,rel,refLen);
                        long long comRefLen=max(negKLen,posKLen);
                        double overPercentage=comRefLen*1.0/refLen;
                        printf("pppppppercentage:%.2lf\n",overPercentage);
                        if(overPercentage>=0.7){
                            thisInvPass=false;
                        }
                        delete []refReg1;
                        delete []refReg2;
                    }
                    if(!thisInvPass){
                        break;
                    }
                }
                if(thisInvPass) {
                    oneInvPass=true;
                    break;
                }
            }
        }
        if(!oneInvPass){
            FILE *fp=fopen(outF,"a");
            fprintf(fp,"%s\t%d\tnoValidInvSig\n",read.qname.c_str(),omp_get_thread_num());
            fclose(fp);
            return;
        }
        return;
        if(refOverlap){
            if(refMaxPos-refMinPos>=1.5*read.len||read.len>80000){
                FILE *fp=fopen(outF,"a");
                fprintf(fp,"%s\t%d\tlargeRegion\n",read.qname.c_str(),omp_get_thread_num());
                fclose(fp);
                return;
            }
            char *refSeq = new char[refMaxPos - refMinPos + 10];
            char *refPint = refSeq;
            #pragma omp critical
            {
                faFile ref(refCache);
                ref.loadAgctByChr(chr, refFile);
                for (char *i = refCache + refMinPos; i <= refCache + refMaxPos; i++, refPint++) {
                    *refPint = *i;
                }
            }
            *refPint=0;
            long long refLen=strlen(refSeq);
            seq1=refSeq;
            seq2=read.seqCache;
            int dirK=0,comK=0;
            vector<pair<int, int> >  rel[3],mainRel[3];
            //首先进行直接比对，根据直接比对中两种斜率的分布结果，确定最适合直接比对的斜率
            rel[1]=getXY(true,true);
            sort(rel[1].begin(),rel[1].end());
            getCluster(rel[1]);
            long long negKLen=getClusterRefLen(posId,idNum,rel[1],refLen);
            long long posKLen=getClusterRefLen(posId2,idNum2,rel[1],refLen);
            long long dirRefLen=max(negKLen,posKLen);
            pair<long long,long long> dirSamLen,comSamLen;

            if(negKLen>posKLen){
                dirK=-1;
                storeMainRel(posId,idNum,rel[1],mainRel[1]);
                dirSamLen=getClusterSamLen(posId,idNum,rel[1],read.len);
                copyCluster(1,posId,idNum,segPosId[1],segIdNum[1],rel[1].size());
            }
            else{
                dirK=1;
                storeMainRel(posId2,idNum2,rel[1],mainRel[1]);
                dirSamLen=getClusterSamLen(posId2,idNum2,rel[1],read.len);
                copyCluster(1,posId2,idNum2,segPosId[1],segIdNum[1],rel[1].size());
            }

            //进行反向互补比对，并确定其对应的斜率
            rel[2]=getXY(true,false);
            sort(rel[2].begin(),rel[2].end());
            getCluster(rel[2]);
            negKLen=getClusterRefLen(posId,idNum,rel[2],refLen);
            posKLen=getClusterRefLen(posId2,idNum2,rel[2],refLen);
            long long comRefLen=max(negKLen,posKLen);
            if(negKLen>posKLen){
                comK=-1;
                storeMainRel(posId,idNum,rel[2],mainRel[2]);
                comSamLen=getClusterSamLen(posId,idNum,rel[2],read.len);
                copyCluster(2,posId,idNum,segPosId[2],segIdNum[2],rel[2].size());
            }
            else{
                comK=1;
                storeMainRel(posId2,idNum2,rel[2],mainRel[2]);
                comSamLen=getClusterSamLen(posId2,idNum2,rel[2],read.len);
                copyCluster(2,posId2,idNum2,segPosId[2],segIdNum[2],rel[2].size());
            }

            //如果两种比对对应的斜率是相同的，那么这条read不符合条件，退出
            if(comK==dirK){
                FILE *fp=fopen(outF,"a");
                fprintf(fp,"%s\t%d\tdirectoryError\n",read.qname.c_str(),omp_get_thread_num());
                fclose(fp);
                delete[] refSeq;
                return;
            }

            //根据直接和反向互补比对的在sam上的长度来确定哪方面是主要比对
            int mainK=0,mainId=0;
            if(dirSamLen.second>=comSamLen.second) {
                mainK=dirK;
                mainId=1;
            }else{
                mainK=comK;
                mainId=2;
            }
            int secK=mainK*(-1),secId=3-mainId;
            refBoundaryPos=refMaxPos;
            //getEndBoundaryInfo(chr,refBoundaryPos-100,refBoundaryPos+100,read);

            //保证主要比对的斜率是1，进行标准化工作
            if(-1==mainK){
                refBoundaryPos=refMinPos;
                int siz=rel[mainId].size();
                for(int i=0;i<siz;i++){
                    rel[mainId][i].first=int(refLen-rel[mainId][i].first);
                }
                siz=rel[secId].size();

                for(int i=0;i<siz;i++){
                    rel[secId][i].first=int(refLen-rel[secId][i].first);
                }
                mainK*=-1;
                secK*=-1;
            }

            //在主要比对和次要比对中找最长的比对，用来进行线性回归，以确定端点处的距离
            vector< pair<int,int> > maxPos[3],polePos[3];
            maxPos[mainId]=getMaxClusterOnSam(segPosId[mainId],segIdNum[mainId],rel[mainId]);
            polePos[mainId]=getPolePos(maxPos[mainId]);
            pair<double,double> linearPara[3];
            if(maxPos[mainId].size()<5){
                FILE *fp=fopen(outF,"a");
                fprintf(fp,"%s\t%d\tnoMainMax\n",read.qname.c_str(),omp_get_thread_num());
                fclose(fp);
                return;
            }
            #pragma omp critical
            {
                linearPara[mainId] = getLinearRegressionPara(maxPos[mainId]);
            }
            pair<double,double>mainEndPos=make_pair(refLen,linearPara[mainId].second*refLen+linearPara[mainId].first);
            vector<pair<int,int> >  secPos;
            for(auto i:rel[secId]){
                if(i.second>=mainEndPos.second){
                    secPos.push_back(i);
                }
            }
            sort(secPos.begin(),secPos.end());
            getCluster(secPos,secK,secPosId,secIdNum);
            maxPos[secId]=getMaxClusterOnSam(secPosId,secIdNum,secPos);
            polePos[secId]=getPolePos(maxPos[secId]);
            if(maxPos[secId].size()<5){
                FILE *fp=fopen(outF,"a");
                fprintf(fp,"%s\t%d\tnoSecMax\n",read.qname.c_str(),omp_get_thread_num());
                fclose(fp);
                return;
            }
            #pragma omp critical
            {
                linearPara[secId] = getLinearRegressionPara(maxPos[secId]);
            }
            pair<double,double>secStaPos=make_pair(refLen,linearPara[secId].second*refLen+linearPara[secId].first);
            FILE *fp=fopen(outF,"a");
            fprintf(fp,"%s\t%d\tdis\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",\
                                  read.qname.c_str(),omp_get_thread_num(),secStaPos.second-mainEndPos.second,(secStaPos.second-mainEndPos.second)/read.len,\
                                  (polePos[secId][1].second-polePos[secId][0].second)/read.len,\
                                  (refMaxPos-refMinPos-polePos[secId][1].first)/(refMaxPos-refMinPos),\
                                  (polePos[mainId][1].second-polePos[mainId][0].second)/read.len,\
                                  (refMaxPos-refMinPos-polePos[mainId][1].first)/(refMaxPos-refMinPos),\
                                  linearPara[secId].second,linearPara[mainId].second);
            fclose(fp);
            if(plot){
                //plotDotplot(rel[mainId],1,rel[secId],2,read.qname,refLen,read.len,0);
                //plotDotplot(maxPos[mainId],1,maxPos[secId],2,read.qname,refLen,read.len,1);
                //plotRefBoundary(refBoundaryPos,refCache,read.qname);
            }
        }
    }
}
*/
alignManager::alignManager(oneFastq &read){
    posId=NULL;
    idNum=NULL;
    posId2=NULL;
    idNum2=NULL;
    secPosId=NULL;
    secIdNum=NULL;
    segTreeL=segTreeR=segNodesHave=segNodesFull=NULL;
    memset(segPosId,0,sizeof(segPosId));
    memset(segIdNum,0,sizeof(segIdNum));
    myFastq=read;
}
void alignManager::getCluster(vector<pair<int, int>> &pos) {
    if(posId!=NULL){
        delete[] posId;
        delete[] posId2;
    }
    if(idNum!=NULL){
        delete[] idNum;
        delete[] idNum2;
    }
    int siz=pos.size();
    posId=new int[siz+5];
    idNum=new int[siz+5];
    posId2=new int[siz+5];
    idNum2=new int[siz+5];
    for(int i=0;i<siz;i++){
        posId[i]=i;
        idNum[i]=0;
        posId2[i]=i;
        idNum2[i]=0;
    }
    for(int i=0;i<siz;i++){
        for(int j=i+1;j<siz;j++){
            if(pos[j].second<pos[i].second&&pos[i].second-pos[j].second<=50&&pos[j].first-pos[i].first<=50){
                posId[getClusterId(j,posId)]=getClusterId(i,posId);
            }
            if(pos[j].second>pos[i].second&&pos[j].second-pos[i].second<=50&&pos[j].first-pos[i].first<=50){
                posId2[getClusterId(j,posId2)]=getClusterId(i,posId2);
            }
        }
    }
    for(int i=0;i<siz;i++){
        posId[i]=getClusterId(i,posId);
        idNum[posId[i]]++;
        posId2[i]=getClusterId(i,posId2);
        idNum2[posId2[i]]++;
    }
}

void alignManager::getCluster(vector< pair<int,int> > &pos,int k,int *&index,int *&num){
    if(index!=NULL){
        delete[] index;
    }
    if(num!=NULL){
        delete[] num;
    }
    int siz=pos.size();
    index=new int[siz+5];
    num=new int[siz+5];
    for(int i=0;i<siz;i++){
        index[i]=i;
        num[i]=0;
    }
    if(-1==k){
        for(int i=0;i<siz;i++){
            for(int j=i+1;j<siz;j++){
                if(pos[j].second<pos[i].second&&pos[i].second-pos[j].second<=50&&pos[j].first-pos[i].first<=50){
                    index[getClusterId(j,index)]=getClusterId(i,index);
                }
            }
        }
    }
    else{
        for(int i=0;i<siz;i++){
            for(int j=i+1;j<siz;j++){
                if(pos[j].second>pos[i].second&&pos[j].second-pos[i].second<=50&&pos[j].first-pos[i].first<=50){
                    index[getClusterId(j,index)]=getClusterId(i,index);
                }
            }
        }
    }
    for(int i=0;i<siz;i++){
        index[i]=getClusterId(i,index);
        num[index[i]]++;
    }
}

alignManager::~alignManager() {
    if(posId!=NULL){
        delete[] posId;
        delete[] posId2;
    }
    if(idNum!=NULL){
        delete[] idNum;
        delete[] idNum2;
    }
    if(segTreeL!=NULL){
        delete[] segTreeL;
        delete[] segTreeR;
        delete[] segNodesFull;
        delete[] segNodesHave;
    }
    for(int i=0;i<3;i++){
        if(segIdNum[i]!=NULL){
            delete[] segIdNum[i];
        }
        if(segPosId[i]!=NULL){
            delete[] segPosId[i];
        }
    }
    if(secPosId!=NULL){
        delete[] secPosId;
        delete[] secIdNum;
    }
}

int alignManager::getClusterId(int tar,int *fa) {
    if(fa[tar]!=tar){
        return fa[tar]=getClusterId(fa[tar],fa);
    }
    return tar;
}

long long alignManager::getClusterRefLen(int *Index,int *numIndex,vector<pair<int, int> > &rel,int refLen) {
    initSegTree(refLen+5);
    int *refLeftBound,*refRightBound;
    int siz=rel.size();
    refLeftBound=new int[siz+10];
    refRightBound=new int[siz+10];
    int init_max=1e9,init_min=-1;
    for(int i=0;i<siz;i++){
        refLeftBound[i]=init_max;
        refRightBound[i]=init_min;
    }
    int p=0;
    for(auto i:rel){
        if(numIndex[Index[p]]>=10) {
            refLeftBound[Index[p]] = min(refLeftBound[Index[p]], i.first);
            refRightBound[Index[p]] = max(refRightBound[Index[p]], i.first);
        }
        p++;
    }
    long long ans=0;
    for(int i=0;i<siz;i++){
        if(refRightBound[i]-refLeftBound[i]>0){
            dfsInsSeg(refLeftBound[i]+1,refRightBound[i]+1,1);
        }
    }
    delete[] refLeftBound;
    delete[] refRightBound;
    return getTotalSegLen(1,refLen+1,1);
}

pair<long long,long long> alignManager::getClusterSamLen(int *Index,int *numIndex,vector<pair<int, int> > &rel,int samLen) {
    initSegTree(samLen+5);
    int *samLeftBound,*samRightBound;
    int siz=rel.size();
    samLeftBound=new int[siz+10];
    samRightBound=new int[siz+10];
    int init_max=1e9,init_min=-1;
    for(int i=0;i<siz;i++){
        samLeftBound[i]=init_max;
        samRightBound[i]=init_min;
    }
    int p=0;
    for(auto i:rel){
        if(numIndex[Index[p]]>=10) {
            samLeftBound[Index[p]] = min(samLeftBound[Index[p]], i.second);
            samRightBound[Index[p]] = max(samRightBound[Index[p]], i.second);
        }
        p++;
    }
    long long ans=0;
    for(int i=0;i<siz;i++){
        if(samRightBound[i]-samLeftBound[i]>0){
            dfsInsSeg(samLeftBound[i]+1,samRightBound[i]+1,1);
        }
    }
    delete[] samLeftBound;
    delete[] samRightBound;
    return make_pair(getTotalSegLen(1,samLen/2,1),getTotalSegLen(1,samLen,1));
}

void alignManager::printPos(int *Index, int *numIndex,vector<pair<int, int> > &pos,int type) {
    int p=0;
    for(auto i:pos){
        if(numIndex[Index[p]]>=10){
            printf("%d\t%d\t%d\n",i.first,i.second,type);
        }
        p++;
    }
}

void alignManager::initSegTree(int range) {
    if(segTreeL!=NULL){
        delete[] segTreeL;
        delete[] segTreeR;
        delete[] segNodesHave;
        delete[] segNodesFull;
    }
    int treeSiz=(range<<2);
    segTreeL=new int[treeSiz];
    segTreeR=new int[treeSiz];
    segNodesFull=new int[treeSiz];
    segNodesHave=new int[treeSiz];
    dfsInitSegTree(1,range,1);
}

void alignManager::dfsInitSegTree(int l, int r, int id) {
    segTreeL[id]=l;
    segTreeR[id]=r;
    segNodesHave[id]=0;
    segNodesFull[id]=0;
    if(l!=r) {
        int mid=((l+r)>>1);
        dfsInitSegTree(l,mid,id<<1);
        dfsInitSegTree(mid+1,r,(id<<1)|1);
    }
}

void alignManager::dfsInsSeg(int l, int r, int id) {
    if(segNodesFull[id]==1){
        return;
    }
    if(segTreeL[id]==l&&segTreeR[id]==r){
        segNodesFull[id]=1;
        segNodesHave[id]=1;
        return;
    }
    segNodesHave[id]=1;
    int mid=(segTreeL[id]+segTreeR[id])>>1;
    if(r<=mid){
        dfsInsSeg(l,r,id<<1);
    }
    else if(l>mid){
        dfsInsSeg(l,r,(id<<1)|1);
    }
    else{
        dfsInsSeg(l,mid,id<<1);
        dfsInsSeg(mid+1,r,(id<<1)|1);
    }
}

int alignManager::getTotalSegLen(int l,int r,int id) {
    if(segNodesHave[id]==0){
        return 0;
    }
    if(segNodesFull[id]==1){
        return r-l+1;
    }
    int mid=(segTreeL[id]+segTreeR[id])>>1;
    if(r<=mid){
        return getTotalSegLen(l,r,id<<1);
    }
    else if(l>mid){
        return getTotalSegLen(l,r,(id<<1)|1);
    }
    else{
        return getTotalSegLen(l,mid,id<<1)+getTotalSegLen(mid+1,r,(id<<1)|1);
    }
}

void alignManager::storeMainRel(int *Index, int *numIndex, vector<pair<int, int>> &rel, vector<pair<int, int> > &mainRel) {
    int p=0;
    for(auto i:rel){
        if(numIndex[Index[p]]>=10){
            mainRel.push_back(i);
        }
        p++;
    }
}

void alignManager::copyCluster(int segId, int *index, int *num,int * &des_index,int * &des_num, int relLen) {
    if(des_index!=NULL){
        delete[] des_index;
        delete[] des_num;
    }
    des_index=new int[relLen+5];
    des_num=new int[relLen+5];
    for(int i=0;i<relLen;i++) {
        des_index[i] = index[i];
        des_num[i] = num[i];
    }
}

void alignManager::printPos(vector<pair<int, int>> &pos,int type) {
    for(auto i:pos){
        printf("%d\t%d\t%d\n",i.first,i.second,type);
    }
}

vector<pair<int, int> > alignManager::getMaxClusterOnSam(int *Index, int *numIndex, vector<pair<int, int>> &pos) {
    int *samLeftBound,*samRightBound;
    int siz=pos.size();
    samLeftBound=new int[siz+10];
    samRightBound=new int[siz+10];
    int init_max=1e9,init_min=-1;
    for(int i=0;i<siz;i++){
        samLeftBound[i]=init_max;
        samRightBound[i]=init_min;
    }
    int p=0;
    for(auto i:pos){
        if(numIndex[Index[p]]>=10) {
            samLeftBound[Index[p]] = min(samLeftBound[Index[p]], i.second);
            samRightBound[Index[p]] = max(samRightBound[Index[p]], i.second);
        }
        p++;
    }
    long long ans=-1;
    int maxRange=-1;
    for(int i=0;i<siz;i++){
        if(samRightBound[i]-samLeftBound[i]>maxRange){
            maxRange=samRightBound[i]-samLeftBound[i];
            ans=i;
        }
    }
    vector< pair<int,int>> selectedRel;
    for(int i=0;i<siz;i++){
        if(Index[i]==ans){
            selectedRel.push_back(pos[i]);
        }
    }
    delete[] samLeftBound;
    delete[] samRightBound;
    return selectedRel;
}

pair<double, double> alignManager::getLinearRegressionPara(vector<pair<int, int> > &pos) {
    int siz=pos.size();
    Rcpp::NumericMatrix r_pos(pos.size(),2);
    for(int i=0;i<siz;i++){
        r_pos(i,0)=pos[i].first;
        r_pos(i,1)=pos[i].second;
    }
    R["pos"]=r_pos;
    R.parseEvalQ("dat<-data.frame(pos)");
    R.parseEvalQ("colnames(dat)<-c('a','b')");
    R.parseEvalQ("model<-lm(b~a,data=dat)");
    vector<double> rel=Rcpp::as<vector<double> >(R.parseEval("model$coefficient"));
    //printf("%.2lf\t%.2lf\n",rel[0],rel[1]);
    return make_pair(rel[0],rel[1]);
}

void alignManager::plotDotplot(vector<pair<int, int> > pointSet1, int colTyp1, vector<pair<int, int> > pointSet2, int colTyp2,string &qname,int refLen,int readLen,int id) {
    int siz1=pointSet1.size();
    int siz2=pointSet2.size();
    Rcpp::NumericMatrix r_pos(siz1+siz2,3);
    for(int i=0;i<siz1;i++){
        r_pos(i,0)=pointSet1[i].first;
        r_pos(i,1)=pointSet1[i].second;
        r_pos(i,2)=colTyp1;
    }
    for(int i=0;i<siz2;i++){
        r_pos(i+siz1,0)=pointSet2[i].first;
        r_pos(i+siz1,1)=pointSet2[i].second;
        r_pos(i+siz1,2)=colTyp2;
    }
    R["pos"]=r_pos;
    char cacheId[20];
    sprintf(cacheId,"_%d.png",id);
    string plotFileName=string("/media/ttbond/642A8ADA2A8AA91E/fast5Test/scr/custom_dotplot/")+qname+cacheId;
    char cache[500];
    int width=refLen,height=readLen;
    if(refLen>=10000||readLen>=10000){
        width=refLen/10;
        height=readLen/10;
    }
    sprintf(cache,"png(file='%s',width=%d,height=%d,res=300)\n\
                   print(p)\n\
                   dev.off()",plotFileName.c_str(),width,height);
    char scale_limit_cache[1000];
    sprintf(scale_limit_cache,"scale_x_continuous(expand = c(0, 0),limits=c(%d,%d))+\
                               scale_y_continuous(expand = c(0, 0),limits=c(%d,%d))",0,refLen,0,readLen);
    R.parseEvalQ("dat<-data.frame(pos)");
    R.parseEvalQ("colnames(dat)<-c('x','y','typ')");
    R.parseEvalQ("library(Hmisc)");
    R.parseEvalQ("p <- ggplot(data=dat,aes(x=x,y=y))+\
                    geom_point(aes(color=factor(typ)),size=0.4)+\
                    geom_smooth(method='lm',formula=y~x,aes(color=factor(typ)),fullrange=TRUE,size=0.2)+\
                    theme(legend.position=\"none\")+\
                    theme_bw()+\
                    theme(legend.position=\"none\")+\
                    theme(axis.title.x=element_blank(),\
                          axis.text.x=element_blank(),\
                          axis.ticks.x=element_blank())+\
                    theme(axis.title.y=element_blank(),\
                          axis.text.y=element_blank(),\
                          axis.ticks.y=element_blank())+\
                    "+string(scale_limit_cache));
    R.parseEvalQ(cache);
}

void alignManager::plotRefBoundary(long long refBoundaryPos, char *refCache,string &qname) {
    FILE *fp[2];
    string dir("/media/ttbond/642A8ADA2A8AA91E/fast5Test/scr/custom_dotplot");
    string outP=dir+"/"+qname+"_2.png";
    fp[0]=fopen((dir+"/tmp0.fq").c_str(),"w");
    fp[1]=fopen((dir+"/tmp1.fq").c_str(),"w");
    for(int i=0;i<2;i++){
        fprintf(fp[i],">%d\n",i);
    }
    for(char *i=refCache+refBoundaryPos-100;i<refCache+refBoundaryPos+100;i++){
        for(int j=0;j<2;j++){
            fprintf(fp[j],"%c",*i);
        }
    }
    for(int i=0;i<2;i++){
        fprintf(fp[i],"\n");
        fclose(fp[i]);
    }
    string dotplotCommand=string("java -cp /home/ttbond/softwares/gepard/dist/Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 ")+
                                  dir+"/tmp0.fq"+string(" -seq2 ")+dir+"/tmp1.fq"+string(" -matrix /home/ttbond/softwares/gepard/resources/matrices/edna.mat")+
                                  string(" -maxwidth 2000 -maxheight 2000 -outfile ")+outP+string(" -silent");
    system(dotplotCommand.c_str());
}

void alignManager::getEndBoundaryInfo(int chr,long long st,long long ed,oneFastq &read) {
    detectRegion region(cache,chr,st,ed,OTHER);
    FILE *fp=fopen("/media/ttbond/642A8ADA2A8AA91E/fast5Test/scr/repLengthOfCom.dat","a");
    fprintf(fp,"%s\t",read.qname.c_str());
    region.getDirectRepeatScore(fp);
    region.getReverseComScore(fp);
    region.getMirrorRepeatScore(fp);
    fprintf(fp,"\n");
    fclose(fp);
}

vector<pair<int, int> > alignManager::getPolePos(vector<pair<int, int> > &posSet) {
    int minx=1e9,miny=1e9,maxx=-1,maxy=-1;
    for(auto i:posSet){
        minx=min(minx,i.first);
        maxx=max(maxx,i.first);
        miny=min(miny,i.second);
        maxy=max(maxy,i.second);
    }
    vector<pair<int,int> >rel;
    rel.push_back(make_pair(minx,miny));
    rel.push_back(make_pair(maxx,maxy));
    return rel;
}

void alignManager::getDotplot(char *refFile, const char *bamFileName, char *refCache, const char *outD,const char *wd) {
    bamFile bamF(bamFileName);
    vector<alignInfo> posReg=bamF.simpleQuery(myFastq.qname.c_str());
    map<int,vector<alignInfo> > posRel;
    for(auto i:posReg){
        if(i.refInfo.chr<=0||i.refInfo.chr>23){
            continue;
        }
        posRel[i.refInfo.chr].push_back(i);
    }
    for(auto i:posRel){
        int minn=1e9,maxx=-1;
        for(auto j:i.second){
            minn=min((int)j.refInfo.st,minn);
            maxx=max((int)j.refInfo.ed,maxx);
        }
        minn-=15000;
        maxx+=15000;
        char dotplotF[1000];
        sprintf(dotplotF,"%s_%d_%d_%d.png",myFastq.qname.c_str(),i.first,minn,maxx);
        puts((string(outD)+"/tmpRef.fa").c_str());
        FILE *fp=fopen((string(outD)+"/tmpRef.fa").c_str(),"w");
        fprintf(fp,">%d_%d_%d\n",i.first,minn,maxx);
        faFile ref(refCache);
        ref.loadAgctByChr(i.first, refFile);
        for (char *ch = refCache + minn; ch<= refCache + maxx; ch++) {
            fprintf(fp,"%c",*ch);
        }
        fprintf(fp,"\n");
        fclose(fp);
        FILE *fp2=fopen((string(outD)+"/tmpRead.fa").c_str(),"w");
        fprintf(fp2,">%s\n%s\n",myFastq.qname.c_str(),myFastq.seqCache);
        fclose(fp2);
        string dir(wd);
        string outP=string(outD)+"/"+string(dotplotF);
        string dotplotCommand=string("java -cp /home/xutun/gepard-master/dist/Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 ")+
                                  dir+"/tmpRef.fa"+string(" -seq2 ")+dir+"/tmpRead.fa"+string(" -matrix /home/xutun/gepard-master/resources/matrices/edna.mat")+
                                  string(" -maxwidth 2000 -maxheight 2000 -outfile ")+outP+string(" -silent");
        system(dotplotCommand.c_str());
    }
}

basicInfo alignManager::getRefPosFromAlign(alignInfo source,alignInfo tar,bool rev) {
    if(rev){
        swap(source.refInfo.st,source.refInfo.ed);
    }
    long double k=((long double)source.readInfo.st-source.readInfo.ed)/(source.refInfo.st-source.refInfo.ed);
    long double b=((long double)source.readInfo.st)-k*(source.refInfo.st);
    long double tarRefSt=((long double)tar.readInfo.st-b)/k;
    long double tarRefEd=((long double)tar.readInfo.ed-b)/k;
    if(rev){
        swap(tarRefSt,tarRefEd);
    }
    return basicInfo(source.refInfo.chr,(long long)tarRefSt,(long long)tarRefEd);
}
///*
//旧版的主函数体
alignManager::alignManager(oneFastq &read, char *refFile, const char *bamFileName,char *refCache,const char *outF,bool plot) {
    if(plot){
        R.parseEvalQ("library('ggplot2')");
    }
    myFastq=read;
    cache=refCache;
    posId=NULL;
    idNum=NULL;
    posId2=NULL;
    idNum2=NULL;
    secPosId=NULL;
    secIdNum=NULL;
    segTreeL=segTreeR=segNodesHave=segNodesFull=NULL;
    memset(segPosId,0,sizeof(segPosId));
    memset(segIdNum,0,sizeof(segIdNum));
    //截取ref，并准备好readSeq
    bamFile bamF(bamFileName);
    int readLen=strlen(read.seqCache);
    vector<basicInfo> posReg=bamF.simpleQueryBasic(read.qname.c_str());
    long long refMinPos=1e18,refMaxPos=-1;
    int chr=-1;
    map<int,vector<basicInfo> > posInfo;
    for(auto info:posReg){
        if(info.chr<1||info.chr>25){
            continue;
        }
        posInfo[info.chr].push_back(info);
    }
    for(auto info:posInfo){
        refMinPos=1e18,refMaxPos=-1;
        for(auto me:info.second){
            chr=me.chr;
            refMinPos=min(min(refMinPos,me.st),me.ed);
            refMaxPos=max(max(refMaxPos,me.st),me.ed);
        }
        if(refMaxPos-refMinPos>=1.5*read.len||read.len>80000){
            FILE *fp=fopen(outF,"a");
            fprintf(fp,"%s\t%d\tlargeRegion\n",read.qname.c_str(),omp_get_thread_num());
            fclose(fp);
            return;
        }
        char *refSeq = new char[refMaxPos - refMinPos + 10];
        char *refPint = refSeq;
        #pragma omp critical
        {
            faFile ref(refCache);
            ref.loadAgctByChr(chr, refFile);
            for (char *i = refCache + refMinPos; i <= refCache + refMaxPos; i++, refPint++) {
                *refPint = *i;
            }
        }
        bornT=clock();
        *refPint=0;
        long long refLen=strlen(refSeq);
        seq1=refSeq;
        seq2=read.seqCache;
        int dirK=0,comK=0;
        vector<pair<int, int> >  rel[3],mainRel[3];
        //首先进行直接比对，根据直接比对中两种斜率的分布结果，确定最适合直接比对的斜率
        rel[1]=getXY(true,true);
        if(rel[1].size()>5e8){
            FILE *fp=fopen(outF,"a");
            fprintf(fp,"%s\t%d\tTimeError\n",read.qname.c_str(),omp_get_thread_num());
            fclose(fp);
        }
        sort(rel[1].begin(),rel[1].end());
        getCluster(rel[1]);
        long long negKLen=getClusterRefLen(posId,idNum,rel[1],refLen);
        long long posKLen=getClusterRefLen(posId2,idNum2,rel[1],refLen);
        long long dirRefLen=max(negKLen,posKLen);
        pair<long long,long long> dirSamLen,comSamLen;

        if(negKLen>posKLen){
            dirK=-1;
            storeMainRel(posId,idNum,rel[1],mainRel[1]);
            dirSamLen=getClusterSamLen(posId,idNum,rel[1],read.len);
            copyCluster(1,posId,idNum,segPosId[1],segIdNum[1],rel[1].size());
        }
        else{
            dirK=1;
            storeMainRel(posId2,idNum2,rel[1],mainRel[1]);
            dirSamLen=getClusterSamLen(posId2,idNum2,rel[1],read.len);
            copyCluster(1,posId2,idNum2,segPosId[1],segIdNum[1],rel[1].size());
        }

        //进行反向互补比对，并确定其对应的斜率
        rel[2]=getXY(true,false);
        if(rel[2].size()>5e8){
            FILE *fp=fopen(outF,"a");
            fprintf(fp,"%s\t%d\tTimeError\n",read.qname.c_str(),omp_get_thread_num());
            fclose(fp);
            return;
        }
        sort(rel[2].begin(),rel[2].end());
        getCluster(rel[2]);
        negKLen=getClusterRefLen(posId,idNum,rel[2],refLen);
        posKLen=getClusterRefLen(posId2,idNum2,rel[2],refLen);
        long long comRefLen=max(negKLen,posKLen);
        if(negKLen>posKLen){
            comK=-1;
            storeMainRel(posId,idNum,rel[2],mainRel[2]);
            comSamLen=getClusterSamLen(posId,idNum,rel[2],read.len);
            copyCluster(2,posId,idNum,segPosId[2],segIdNum[2],rel[2].size());
        }
        else{
            comK=1;
            storeMainRel(posId2,idNum2,rel[2],mainRel[2]);
            comSamLen=getClusterSamLen(posId2,idNum2,rel[2],read.len);
            copyCluster(2,posId2,idNum2,segPosId[2],segIdNum[2],rel[2].size());
        }
        getDotplot(basicInfo(chr,refMinPos,refMaxPos),refSeq,"/home/xutun/nanoporeBase/NA12878/minimap2_chrM/chrM_dotplot");
        //如果两种比对对应的斜率是相同的，那么这条read不符合条件，退出
        if(comK==dirK){
            FILE *fp=fopen(outF,"a");
            fprintf(fp,"%s\t%d\tdirectoryError\n",read.qname.c_str(),omp_get_thread_num());
            fclose(fp);
            delete[] refSeq;
            return;
        }
        //根据直接和反向互补比对的在sam上的长度来确定哪方面是主要比对
        int mainK=0,mainId=0;
        if(dirSamLen.second>=comSamLen.second) {
            mainK=dirK;
            mainId=1;
        }else{
            mainK=comK;
            mainId=2;
        }
        int secK=mainK*(-1),secId=3-mainId;
        refBoundaryPos=refMaxPos;
        //getEndBoundaryInfo(chr,refBoundaryPos-100,refBoundaryPos+100,read);

        //保证主要比对的斜率是1，进行标准化工作
        if(-1==mainK){
            refBoundaryPos=refMinPos;
            int siz=rel[mainId].size();
            for(int i=0;i<siz;i++){
                rel[mainId][i].first=int(refLen-rel[mainId][i].first);
            }
            siz=rel[secId].size();
            for(int i=0;i<siz;i++){
                rel[secId][i].first=int(refLen-rel[secId][i].first);
            }
            mainK*=-1;
            secK*=-1;
        }

        //在主要比对和次要比对中找最长的比对，用来进行线性回归，以确定端点处的距离
        vector< pair<int,int> > maxPos[3],polePos[3];
        maxPos[mainId]=getMaxClusterOnSam(segPosId[mainId],segIdNum[mainId],rel[mainId]);
        polePos[mainId]=getPolePos(maxPos[mainId]);
        pair<double,double> linearPara[3];
        if(maxPos[mainId].size()<5){
            FILE *fp=fopen(outF,"a");
            fprintf(fp,"%s\t%d\tnoMainMax\n",read.qname.c_str(),omp_get_thread_num());
            fclose(fp);
            return;
        }
        #pragma omp critical
        {
            linearPara[mainId] = getLinearRegressionPara(maxPos[mainId]);
        }
        pair<double,double>mainEndPos=make_pair(refLen,linearPara[mainId].second*refLen+linearPara[mainId].first);
        vector<pair<int,int> >  secPos;
        for(auto i:rel[secId]){
            if(i.second>=mainEndPos.second){
                secPos.push_back(i);
            }
        }
        sort(secPos.begin(),secPos.end());
        getCluster(secPos,secK,secPosId,secIdNum);
        maxPos[secId]=getMaxClusterOnSam(secPosId,secIdNum,secPos);
        polePos[secId]=getPolePos(maxPos[secId]);
        if(maxPos[secId].size()<5){
            FILE *fp=fopen(outF,"a");
            fprintf(fp,"%s\t%d\tnoSecMax\n",read.qname.c_str(),omp_get_thread_num());
            fclose(fp);
            return;
        }
        #pragma omp critical
        {
            linearPara[secId] = getLinearRegressionPara(maxPos[secId]);
        }
        pair<double,double>secStaPos=make_pair(refLen,linearPara[secId].second*refLen+linearPara[secId].first);

        #pragma omp critical
        {
            FILE *fp = fopen(outF, "a");
            long long boundaryPos;
            if (mainK == 1) {
                boundaryPos = refMaxPos;
            } else {
                boundaryPos = refMinPos;
            }
            fprintf(fp, "%s\t%d\tdis\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t", \
                                  read.qname.c_str(), omp_get_thread_num(), secStaPos.second - mainEndPos.second,
                                  (secStaPos.second - mainEndPos.second) / read.len, \
                                  (polePos[secId][1].second - polePos[secId][0].second) / read.len, \
                                  (refMaxPos - refMinPos - polePos[secId][1].first) / (refMaxPos - refMinPos), \
                                  (polePos[mainId][1].second - polePos[mainId][0].second) / read.len, \
                                  (refMaxPos - refMinPos - polePos[mainId][1].first) / (refMaxPos - refMinPos), \
                                  linearPara[secId].second, linearPara[mainId].second);
            fprintf(fp, "%d\t%lld\n", chr, boundaryPos);
            fclose(fp);
        }
        if(plot){
            //plotDotplot(rel[mainId],1,rel[secId],2,read.qname,refLen,read.len,0);
            //plotDotplot(maxPos[mainId],1,maxPos[secId],2,read.qname,refLen,read.len,1);
            //plotRefBoundary(refBoundaryPos,refCache,read.qname);
        }
        delete[] refSeq;
    }

    return;
}
//*/
///*
//筛选符合model1的qname
alignManager::alignManager(oneFastq &read,const char *bamFileName,const char *outF) {
    posId=NULL;
    idNum=NULL;
    posId2=NULL;
    idNum2=NULL;
    secPosId=NULL;
    secIdNum=NULL;
    segTreeL=segTreeR=segNodesHave=segNodesFull=NULL;
    memset(segPosId,0,sizeof(segPosId));
    memset(segIdNum,0,sizeof(segIdNum));
    //截取ref，并准备好readSeq
    bamFile bamF(bamFileName);
    vector<alignInfo> posReg;
    posReg=bamF.simpleQueryAlign(read.qname.c_str(),read.len);
    map<int,vector<alignInfo> > pos;
    for(auto info:posReg){
        if(info.refInfo.chr<0||info.refInfo.chr>24){
            continue;
        }
        pos[info.refInfo.chr].push_back(info);
    }
    for(auto info:pos){
        int siz=info.second.size();
        bool refOverlap=false;
        for(int i=0;i<siz;i++){
            for(int j=i+1;j<siz;j++){
                if((info.second)[i].refInfo.overlapRate((info.second)[j].refInfo)>0 && (info.second[i].flag&16)!=(info.second[j].flag&16)){
                    refOverlap=true;
                    break;
                }
            }
        }
        if(refOverlap){
            FILE *fp=fopen(outF,"a");
            fprintf(fp,"%s\n",read.qname.c_str());
            fclose(fp);
        }
    }
}

void alignManager::getDotplot(basicInfo refPos, const char *refCache, const char *wd) {
    char dotplotF[1000];
    char coreId[100];
    sprintf(coreId,"%d",omp_get_thread_num());
    sprintf(dotplotF,"%s_%d_%lld_%lld.38.png",myFastq.qname.c_str(),refPos.chr,refPos.st,refPos.ed);
    string tmpRefFa=string(wd)+"/"+string(coreId)+"tmpRef.fa";
    string tmpReadFa=string(wd)+"/"+string(coreId)+"tmpRead.fa";
    FILE *fp=fopen(tmpRefFa.c_str(),"w");
    fprintf(fp,">%d_%lld_%lld\n",refPos.chr,refPos.st,refPos.ed);
    fprintf(fp,"%s\n",refCache);
    fclose(fp);
    FILE *fp2=fopen(tmpReadFa.c_str(),"w");
    fprintf(fp2,">%s\n%s\n",myFastq.qname.c_str(),myFastq.seqCache);
    fclose(fp2);
    string dir(wd);
    string outP=string(wd)+"/"+string(dotplotF);
    string dotplotCommand=string("java -cp /home/xutun/gepard-master/dist/Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 ")+
                                  tmpReadFa+string(" -seq2 ")+tmpRefFa+string(" -matrix /home/xutun/gepard-master/resources/matrices/edna.mat")+
                                  string(" -maxwidth 2000 -maxheight 2000 -outfile ")+outP+string(" -silent");
    system(dotplotCommand.c_str());
}


