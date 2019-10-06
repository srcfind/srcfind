//
// Created by ttbond on 18-10-9.
//

#include "geneRegion.h"

int geneRegion::firstGeneRegion=1;
int geneRegion::refChr=0;
std::map<int,std::string> geneRegion::id2name;
std::set<int> geneRegion::allGeneId;
geneRegion::geneRegion(int init) {
    if(init==1) {
        firstGeneRegion = 1;
        id2name.clear();
        allGeneId.clear();
    }
}
geneRegion::geneRegion(char *source,char *_agct) {
    if(firstGeneRegion==1) {
        char cache[1000], name[100];
        firstGeneRegion = 0;
        int _id;
        FILE *fp = fopen("./geneExons/name2id.dat", "r");
        while (fgets(cache, 1000, fp) != NULL) {
            sscanf(cache, "%d%s", &_id, name);
            id2name[_id] = std::string(name);
        }
    }
    agct=_agct;
    int _chr;
    long long _st,_ed;
    char cache[10];
    sscanf(source,"%s%lld%lld%d",cache,&_st,&_ed,&id);
    _chr=stringChrName2int(cache);
    allGeneId.insert(id);
    addRegion(_chr,_st,_ed);
}

void geneRegion::addRegion(char *source){
    char _chr[10];
    long long _st,_ed;
    sscanf(source,"%s%lld%lld",_chr,&_st,&_ed);
    addRegion(stringChrName2int(_chr),_st,_ed);
}
void geneRegion::addRegion(int _chr,long long _st,long long _ed){
    myInfo.push_back(basicInfo(_chr,_st,_ed,OTHER));
}
void geneRegion::solAndPrintRel(FILE *fp){
    sort(myInfo.begin(),myInfo.end());
    std::vector<basicInfo>::iterator newEnd;
    newEnd=unique(myInfo.begin(),myInfo.end());
    myInfo.erase(newEnd,myInfo.end());
    int exonNum=myInfo.size();
    for(int i=0;i<exonNum;i++){
        if(myInfo[i].chr==24){
            continue;
        }
        if(myInfo[i].chr!=refChr){
            refChr=myInfo[i].chr;
            loadAgctByChr(refChr,"GRCh38.d1.vd1.fa",agct);
        }
        printf("gene:%d\texonNum:%d\t",id,i);
        printf("%lld %lld\n",myInfo[i].st,myInfo[i].ed);
        detectRegion region(agct,myInfo[i].chr,myInfo[i].st,myInfo[i].ed,myInfo[i].sv);
        region.getMirrorRepeatScore();
        region.getDirectRepeatScore();
        region.getReverseComScore();
        myRel.push_back(solRel(region));
    }
    for(int i=0;i<exonNum;i++){
        if(myInfo[i].chr==24){
            continue;
        }
        fprintf(fp,"%s\t%d\t",id2name[id].c_str(),i+1);
        //fprintf(fp,"%d\t%d\t",id,i+1);
        myRel[i].printMe(fp);
        myRel[i].info->printRange(fp);
        fprintf(fp,"\n");
    }
}

void geneRegion::getSeqSta(agctTree *sta) {
    sort(myInfo.begin(),myInfo.end());
    std::vector<basicInfo>::iterator newEnd;
    newEnd=unique(myInfo.begin(),myInfo.end());
    myInfo.erase(newEnd,myInfo.end());
    int exonNum=myInfo.size();
    for(int i=0;i<exonNum;i++){
        if(myInfo[i].chr==24){
            continue;
        }
        if(myInfo[i].chr!=refChr){
            refChr=myInfo[i].chr;
            loadAgctByChr(refChr,"GRCh38.d1.vd1.fa",agct);
        }
        printf("gene:%d\texonNum:%d\t",id,i);
        printf("%lld %lld\n",myInfo[i].st,myInfo[i].ed);
        detectRegion region(agct,myInfo[i].chr,myInfo[i].st,myInfo[i].ed,myInfo[i].sv);
        region.getMirrorRepeatScore();
        sta->addSeq(region.detectRel,region.length,region.agct);
    }

}