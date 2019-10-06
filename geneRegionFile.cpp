//
// Created by ttbond on 18-10-9.
//

#include "geneRegionFile.h"
geneRegionFile::geneRegionFile(char *file,char *agct) {
    FILE *fp=fopen(file,"r");
    char cache[1000],cache2[10];
    long long tmp;
    int tmpId;
    geneRegion init(1);
    std::vector<geneRegion>::iterator it;
    int totalGeneNum=0;
    myT=NULL;
    while(fgets(cache,1000,fp)!=NULL){
        sscanf(cache,"%s%lld%lld%d",cache2,&tmp,&tmp,&tmpId);
        if(init.allGeneId.find(tmpId)==init.allGeneId.end()){
            myGeneRegion.push_back(geneRegion(cache,agct));
            it=myGeneRegion.end()-1;
            id2ind[tmpId]=totalGeneNum;
            totalGeneNum++;
        }
        else{
            myGeneRegion[id2ind[tmpId]].addRegion(cache);
        }
    }
}

void geneRegionFile::solAndPrintRel() {
    std::vector<geneRegion>::iterator it;
    FILE *fp=fopen("./geneExons/geneExonNoCancer.gms","w");
    for(it=myGeneRegion.begin();it!=myGeneRegion.end();it++){
        (*it).solAndPrintRel(fp);
    }
    fclose(fp);
}

void geneRegionFile::getSeqSta() {
    std::vector<geneRegion>::iterator it;
    FILE *fp=fopen("./geneExons/geneExonNoCancer.sta","w");
    myT=new agctTree();
    for(it=myGeneRegion.begin();it!=myGeneRegion.end();it++){
        (*it).getSeqSta(myT);
    }
    myT->printMostSeq(200,15,fp);
    fprintf(fp,"\n");
    fclose(fp);
}
geneRegionFile::~geneRegionFile() {
    if(myT!=NULL){
        myT->~agctTree();
    }
}