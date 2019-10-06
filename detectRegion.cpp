//
// Created by ttbond on 18-9-16.
//
#ifndef SRCFIND_DETECTREGION_CPP
#define SRCFIND_DETECTREGION_CPP
#include "detectRegion.h"

char *detectRegion::base=new char[5];
bool detectRegion::firstRegion=true;
int detectRegion::maxMismatchNum=1;
int detectRegion::maxDetectLen=200;
int detectRegion::minDetectLen=3;
bool detectRegion::firstLineOfFile=true;


detectRegion::detectRegion(char *str,int _chr,long long _st,long long _ed,svType _mySv){
    if(firstRegion){
        base[0]='A';
        base[1]='G';
        base[2]='C';
        base[3]='T';
        base[4]='\0';
        firstRegion=false;
    }
    st=_st;
    ed=_ed;
    chr=_chr;
    mySv=_mySv;
    length=ed-st+1;
    cacheLen=max(20000,length*30+10);
    agct=new char[length+10];
    char *i=agct,*j=str+st-1,*charEd=str+ed-1;
    for(;j<=charEd;j++,i++){
        (*i)=(*j);
    }
    (*i)='\0';
    agctEnd=agct+length;
    firstBasePos=NULL;
    cacheForDfs=new int[cacheLen];
    cacheForDfs2=NULL;
    reverseComScore=NULL;
    directRepScore=NULL;
    mirrorRepScore=NULL;
    firstBasePos2=NULL;
    detectRel=new int[length+10];
}

//region2 should be at the letf region of region1
detectRegion::detectRegion(char *str1,int _chr,long long _st1,long long _ed1,char *str2,int _chr2,long long _st2,long long _ed2,svType _mySv){
    if(firstRegion){
        base[0]='A';
        base[1]='G';
        base[2]='C';
        base[3]='T';
        base[4]='\0';
        firstRegion=false;
    }
    st=_st1;
    ed=_ed1;
    chr=_chr;
    st2=_st2;
    ed2=_ed2;
    chr2=_chr2;
    mySv=_mySv;
    length=(ed-st+1)+(ed2-st2+1);
    cacheLen=length*4+10;
    agct=new char[length+10];
    char *i=agct,*j=str1+st-1,*charEd=str1+ed-1;
    for(;j<=charEd;j++,i++){
        (*i)=(*j);
    }
    j=str2+st2-1;
    charEd=str2+ed2-1;
    for(;j<=charEd;j++,i++){
        (*i)=(*j);
    }
    (*i)='\0';
    agctEnd=agct+length;
    firstBasePos=NULL;
    cacheForDfs=new int[cacheLen];
    reverseComScore=NULL;
    directRepScore=NULL;
    mirrorRepScore=NULL;
    firstBasePos2=NULL;
    detectRel=new int[length+10];
}

detectRegion::detectRegion(char *str,const basicInfo &source) {
    new (this) detectRegion(str,source.chr,source.st,source.ed,source.sv);
}

detectRegion::detectRegion(char *_agct) {
    if(firstRegion){
        base[0]='A';
        base[1]='G';
        base[2]='C';
        base[3]='T';
        base[4]='\0';
        firstRegion=false;
    }
    st=1;
    ed=strlen(_agct);
    chr=-1;
    mySv=OTHER;
    length=ed-st+1;
    cacheLen=max(20000,length*10+10);
    agct=new char[length+10];
    strcpy(agct,_agct);
    agctEnd=agct+length;
    firstBasePos=NULL;
    cacheForDfs=new int[cacheLen];
    reverseComScore=NULL;
    directRepScore=NULL;
    mirrorRepScore=NULL;
    firstBasePos2=NULL;
    detectRel=new int[length+10];
}

void detectRegion::localizeFirstBasePos(){
    releaseFirstBasePos(firstBasePos);
    firstBasePos=new int*[4];
    for(int **i=firstBasePos;i<firstBasePos+4;i++){
        (*i)=new int[cacheLen];
    }
    for(int i=0;i<4;i++){
        indEndOfFirstBasePos[i]=firstBasePos[i];
    }
    for(char *i=agct;i<agctEnd;i++){
        for(char *j=base;j<base+4;j++){
            if((*i)==(*j)){
                int id=j-base;
                *(indEndOfFirstBasePos[id]++)=i-agct;
                break;
            }
        }
    }
}

void detectRegion::releaseFirstBasePos(int **&_firstBasePos){
    if(_firstBasePos!=NULL){
        for(int i=0;i<4;i++){
            delete[] _firstBasePos[i];
        }
    }
    delete[] _firstBasePos;
    _firstBasePos=NULL;
}

void detectRegion::release(double *&toRelease){
    if(toRelease!=NULL){
        delete[] toRelease;
    }
    toRelease=NULL;
}
void detectRegion::release(int *&toRelease){
    if(toRelease!=NULL){
        delete[] toRelease;
    }
    toRelease=NULL;
}
void detectRegion::release(char *&toRelease){
    if(toRelease!=NULL){
        delete[] toRelease;
    }
    toRelease=NULL;
}
void detectRegion::release(int **&toRelease){
    if(toRelease!=NULL) {
        for (int i = 0; i < 4; i++) {
            delete[] *(toRelease + i);
        }
    }
    delete[] toRelease;
    toRelease=NULL;
}

void detectRegion::releaseDetectRel(){
    memset(detectRel,0,(length+10)*sizeof(int));
}

int *detectRegion::getInitMisSt(int iteNum,int **_firstBasePos,int *_cacheForDfs){
    if(_firstBasePos==NULL){
        puts("firstBasePos hasn't been loacalized");
        exit(-1);
    }
    return _cacheForDfs;
}
int *detectRegion::getInitMisEd(int iteNum,int **_firstBasePos,int *_cacheForDfs){
    if(_firstBasePos==NULL){
        puts("firstBasePos hasn't been loacalized");
        exit(-1);
    }
    return _cacheForDfs+(indEndOfFirstBasePos[iteNum]-firstBasePos[iteNum]);
}

void detectRegion::copyFirstBasePos(){
    if(firstBasePos==NULL){
        puts("No firstBasePos to copy");
        exit(-1);
    }
    releaseFirstBasePos(firstBasePos2);
    firstBasePos2=new int*[4];
    for(int **i=firstBasePos2;i<firstBasePos2+4;i++){
        (*i)=new int[cacheLen];
    }
    for(int i=0;i<4;i++){
        indEndOfFirstBasePos2[i]=firstBasePos2[i];
        for(int *j=firstBasePos[i];j<indEndOfFirstBasePos[i];j++){
            *(indEndOfFirstBasePos2[i]++)=(*j);
        }
    }
    release(cacheForDfs2);
    cacheForDfs2=new int[cacheLen];
    memset(cacheForDfs2,0,cacheLen*sizeof(int));
}

void detectRegion::printDetectRel(FILE *fp,bool test){
    int tMax=0;
    if(test) {
        for (int i = 0; i < length; i++) {
            tMax = max(tMax, detectRel[i]);
            printf("%d ", detectRel[i]);
        }
        printf("Total lenght:%d MAX:%d\n",length,tMax);
        return;
    }
    std::vector<int>::iterator it;
    std::vector<int> tmp[maxDetectLen+10];
    for(int i=0;i<length;i++){
        tmp[detectRel[i]].push_back(i);
    }
    int *relEd=detectRel+length;
    for(int *i=detectRel;i<relEd;i++){
        if(tMax<(*i)){
            tMax=*i;
        }
    }
    if(fp==NULL){
        printf("#CHROM CPST CPED ALT ST ED AGCT\n");
        char tmpCache[10];
        for(int *i=detectRel;i<relEd;i++){
            if(*i==tMax){
                int inPos=i-detectRel;
                printf("chr%d\t%lld\t%lld\t%s\t%lld\t%lld\t",chr,st,ed,svType2str(tmpCache,mySv),st+inPos,st+inPos+(*i));
                char *tst=agct+inPos,*ted=agct+inPos+(*i);
                for(;tst<ted;tst++){
                    printf("%c",*tst);
                }
                printf("\n");
            }
        }
    }
    else{
        if(firstLineOfFile) {
            fprintf(fp, "#CHROM CPST CPED ALT ST ED AGCT\n");
            firstLineOfFile = false;
        }
        char tmpCache[10];
        for(int *i=detectRel;i<relEd;i++){
            if(*i==tMax){
                int inPos=i-detectRel;
                fprintf(fp,"chr%d\t%lld\t%lld\t%s\t%lld\t%lld\t",chr,st,ed,svType2str(tmpCache,mySv),st+inPos,st+inPos+(*i));
                char *tst=agct+inPos,*ted=agct+inPos+(*i);
                for(;tst<ted;tst++){
                    fprintf(fp,"%c",*tst);
                }
                fprintf(fp,"\n");
            }
        }

    }
}

void detectRegion::printScore(double *score,double *scoreEd,FILE *fp){
    char chrCache[10];
    if(fp==NULL){
        return;
        printf("%d %lld %lld %s\n",chr,st,ed,svType2str(chrCache,mySv));
        for(double *i=score;i<scoreEd;i++){
            printf("%.0lf ",(*i));
        }
        printf("\n");
        return;
    }
    if(fp!=NULL){
        fprintf(fp,"%d %lld %lld %s\n",chr,st,ed,svType2str(chrCache,mySv));
        for(double *i=score;i<scoreEd;i++){
            fprintf(fp,"%.0lf ",(*i));
        }
        fprintf(fp,"\n");
    }
}



double *detectRegion::statisticScore(double *score,double &sumScore,FILE *fp){
    memset(score,0,length*sizeof(double));
    double *edOfScore=score+length;
    int *tdetectRel=detectRel;
    sumScore=0;
    for(double *i=score;i<edOfScore;i++,tdetectRel++) {
        double *ti = i, *iEd = i + (*tdetectRel);
        for (; ti < iEd; ti++) {
            (*ti) += 1;
        }
        sumScore+=(*i);
    }
    return score;
}

double  *detectRegion::getReverseComScore(FILE *fp){
    release(reverseComScore);
    releaseDetectRel();
    reverseComScore=new double[length+10];
    memset(cacheForDfs,0,cacheLen*sizeof(int));
    localizeFirstBasePos();
    copyFirstBasePos();
    for(int i=0;i<4;i++){
        int comBaseNum=nuc2int(getComBase(int2nuc(i)));
        revComDfs(firstBasePos[i],indEndOfFirstBasePos[i],getInitMisSt(i,firstBasePos,cacheForDfs),getInitMisEd(i,firstBasePos,cacheForDfs),firstBasePos2[comBaseNum],indEndOfFirstBasePos2[comBaseNum],getInitMisSt(comBaseNum,firstBasePos2,cacheForDfs2),getInitMisEd(comBaseNum,firstBasePos2,cacheForDfs2),1);
    }
    //statisticScore(reverseComScore,sumRevComScore);
    //printScore(reverseComScore,reverseComScore+length,fp);
    printMaxScore(detectRel,detectRel+length,fp);
}


double *detectRegion::getDirectRepeatScore(FILE *fp) {
    release(directRepScore);
    releaseDetectRel();
    directRepScore=new double[length+10];
    memset(cacheForDfs,0,cacheLen*sizeof(int));
    localizeFirstBasePos();
    for(int i=0;i<4;i++){
        dirRepDfs(firstBasePos[i],indEndOfFirstBasePos[i],getInitMisSt(i,firstBasePos,cacheForDfs),getInitMisEd(i,firstBasePos,cacheForDfs),1);
    }
    //statisticScore(directRepScore,sumDirRepScore);
    //printScore(directRepScore,directRepScore+length,fp);
    printMaxScore(detectRel,detectRel+length,fp);
    //release(directRepScore);
}

double  *detectRegion::getMirrorRepeatScore(FILE *fp){
    release(mirrorRepScore);
    releaseDetectRel();
    mirrorRepScore=new double[length+10];
    memset(cacheForDfs,0,cacheLen*sizeof(int));
    localizeFirstBasePos();
    copyFirstBasePos();
    for(int i=0;i<4;i++){
        mirRepDfs(firstBasePos[i],indEndOfFirstBasePos[i],getInitMisSt(i,firstBasePos,cacheForDfs),getInitMisEd(i,firstBasePos,cacheForDfs),firstBasePos2[i],indEndOfFirstBasePos2[i],getInitMisSt(i,firstBasePos2,cacheForDfs2),getInitMisEd(i,firstBasePos2,cacheForDfs2),1);
    }
    //statisticScore(mirrorRepScore,sumMirRepScore);
    //printScore(mirrorRepScore,mirrorRepScore+length,fp);
    printMaxScore(detectRel,detectRel+length,fp);
    //release(directRepScore);
}


void detectRegion::dirRepDfs(int *st,int *ed,int *misMatchSt,int *misMatchEd,int depth){
    for(char *nuc=base;nuc<base+4;nuc++) {
        int *nst=ed,*ned=ed;
        int *nmst = misMatchEd, *nmed = misMatchEd;
        for(int *tst=st,*tmst=misMatchSt;tst<ed;tst++,tmst++){
            int partial=(*tst)+depth;
            char *nowBase=agct+partial;
            if(nowBase<agctEnd&&*nowBase!='N'){
                if(*nowBase==*nuc){
                    *(ned++)=*tst;
                    *(nmed++)=*tmst;
                }
                else{
                    if(*tmst<maxMismatchNum){
                        *(ned++)=*tst;
                        *(nmed++)=(*tmst)+1;
                    }
                    else{
                        if(depth>=minDetectLen){
                            detectRel[*tst]=max(depth,detectRel[*tst]);
                        }
                    }
                }
            }
            else{
                if(depth>=minDetectLen){
                    detectRel[*tst]=max(depth,detectRel[*tst]);
                }
            }
        }
        if(ned-nst<=1||depth>maxDetectLen){
            for(;nst<ned;nst++){
                detectRel[*nst]=max(detectRel[*nst],depth);
            }
        }
        else{
            dirRepDfs(nst,ned,nmst,nmed,depth+1);
        }
    }
}

void detectRegion::dirRepDfsRec(int *st,int *ed,int *misMatchSt,int *misMatchEd,int depth){
    for(char *nuc=base;nuc<base+4;nuc++) {
        int *nst=ed,*ned=ed;
        int *nmst = misMatchEd, *nmed = misMatchEd;
        for(int *tst=st,*tmst=misMatchSt;tst<ed;tst++,tmst++){
            int partial=(*tst)+depth;
            char *nowBase=agct+partial;
            if(nowBase<agctEnd&&*nowBase!='N'){
                if(*nowBase==*nuc){
                    *(ned++)=*tst;
                    *(nmed++)=*tmst;
                }
                else{
                    if(*tmst<maxMismatchNum){
                        *(ned++)=*tst;
                        *(nmed++)=(*tmst)+1;
                    }
                }
            }
        }
        if(ned-nst>=2&&depth+1>=maxDetectLen){
            std::vector<basicInfo> tmpGrp;
            for(;nst<ned;nst++){
                tmpGrp.push_back(basicInfo(chr,(*nst)+this->st-1,(*nst)+depth+this->st-1));
                //printf("%d %lld %lld\n",chr,*nst,(*nst)+depth);
            }
            //exit(0);
            recRel.push_back(tmpGrp);
        }
        if(ned-nst>=2&&depth+1<maxDetectLen){
            dirRepDfsRec(nst,ned,nmst,nmed,depth+1);
        }
    }
}

void detectRegion:: revComDfs(int *st,int *ed,int *misMatchSt,int *misMatchEd,int *st2,int *ed2,int *misMatchSt2,int *misMatchEd2,int depth){
    for(char *nuc=base;nuc<base+4;nuc++) {
        int *nst=ed,*ned=ed;
        int *nmst = misMatchEd, *nmed = misMatchEd;
        for(int *tst=st,*tmst=misMatchSt;tst<ed;tst++,tmst++){
            int partial=(*tst)+depth;
            char *nowBase=agct+partial;
            if(nowBase<agctEnd&&*nowBase!='N'){
                if(*nowBase==*nuc){
                    *(ned++)=*tst;
                    *(nmed++)=*tmst;
                }
                else{
                    if(*tmst<maxMismatchNum){
                        *(ned++)=*tst;
                        *(nmed++)=(*tmst)+1;
                    }
                    else{
                        if(depth>=minDetectLen){
                            detectRel[*tst]=max(depth,detectRel[*tst]);
                        }
                    }
                }
            }
            else{
                if(depth>=minDetectLen){
                    detectRel[*tst]=max(depth,detectRel[*tst]);
                }
            }
        }

        char tarBase=getComBase(*nuc);
        int *nst2=ed2,*ned2=ed2;
        int *nmst2 = misMatchEd2, *nmed2 = misMatchEd2;
        for(int *tst=st2,*tmst=misMatchSt2;tst<ed2;tst++,tmst++){
            int partial=(*tst)-depth;
            char *nowBase=agct+partial;
            if(nowBase>=agct&&*nowBase!='N'){
                if(tarBase==*nowBase){
                    *(ned2++)=*tst;
                    *(nmed2++)=*tmst;
                }
                else{
                    if(*tmst<maxMismatchNum){
                        *(ned2++)=*tst;
                        *(nmed2++)=(*tmst)+1;
                    }
                }
            }
        }

        if(nst==ned||nst2==ned2||depth>maxDetectLen){
            for(;nst<ned;nst++){
                detectRel[*nst]=max(detectRel[*nst],depth);
            }
        }
        else{
            revComDfs(nst,ned,nmst,nmed,nst2,ned2,nmst2,nmed2,depth+1);
        }
    }
}

void detectRegion:: mirRepDfs(int *st,int *ed,int *misMatchSt,int *misMatchEd,int *st2,int *ed2,int *misMatchSt2,int *misMatchEd2,int depth){
    for(char *nuc=base;nuc<base+4;nuc++) {
        int *nst=ed,*ned=ed;
        int *nmst = misMatchEd, *nmed = misMatchEd;
        for(int *tst=st,*tmst=misMatchSt;tst<ed;tst++,tmst++){
            int partial=(*tst)+depth;
            char *nowBase=agct+partial;
            if(nowBase<agctEnd&&*nowBase!='N'){
                if(*nowBase==*nuc){
                    *(ned++)=*tst;
                    *(nmed++)=*tmst;
                }
                else{
                    if(*tmst<maxMismatchNum){
                        *(ned++)=*tst;
                        *(nmed++)=(*tmst)+1;
                    }
                    else{
                        if(depth>=minDetectLen){
                            detectRel[*tst]=max(depth,detectRel[*tst]);
                        }
                    }
                }
            }
            else{
                if(depth>=minDetectLen){
                    detectRel[*tst]=max(depth,detectRel[*tst]);
                }
            }
        }
        int *nst2=ed2,*ned2=ed2;
        int *nmst2 = misMatchEd2, *nmed2 = misMatchEd2;
        for(int *tst=st2,*tmst=misMatchSt2;tst<ed2;tst++,tmst++){
            int partial=(*tst)-depth;
            char *nowBase=agct+partial;
            if(nowBase>=agct&&*nowBase!='N'){
                if(*nuc==*nowBase){
                    *(ned2++)=*tst;
                    *(nmed2++)=*tmst;
                }
                else{
                    if(*tmst<maxMismatchNum){
                        *(ned2++)=*tst;
                        *(nmed2++)=(*tmst)+1;
                    }
                }
            }
        }

        if(nst==ned||nst2==ned2||depth>maxDetectLen){
            for(;nst<ned;nst++){
                detectRel[*nst]=max(detectRel[*nst],depth);
            }
        }
        else{
            mirRepDfs(nst,ned,nmst,nmed,nst2,ned2,nmst2,nmed2,depth+1);
        }
    }
}
bool detectRegion::operator <(detectRegion &right){
    if(chr!=right.chr){
        return chr<right.chr;
    }
    else if(st!=right.st){
        return st<right.st;
    }
    else{
        return ed<right.ed;
    }
}

detectRegion::~detectRegion(){
    release(reverseComScore);
    release(directRepScore);
    release(mirrorRepScore);
    release(agct);
    release(firstBasePos);
    release(firstBasePos2);
    release(cacheForDfs);
    release(cacheForDfs2);
    release(detectRel);
}

std::vector< std::vector<basicInfo> > detectRegion::getPotDirPos(){
    recRel.clear();
    memset(cacheForDfs,0,cacheLen*sizeof(int));
    localizeFirstBasePos();
    for(int i=0;i<4;i++){
        dirRepDfsRec(firstBasePos[i],indEndOfFirstBasePos[i],getInitMisSt(i,firstBasePos,cacheForDfs),getInitMisEd(i,firstBasePos,cacheForDfs),1);
    }
    return recRel;
}

int detectRegion::printMaxScore(int *score,int *scoreEd, FILE *fp) {
    int rel=-1;
    for(int *i=score;i<scoreEd;i++){
        rel=max(rel,(*i));
    }
    fprintf(fp,"%d\t",rel);
    printf("%d\t",rel);
}

#endif