//
// Created by ttbond on 18-9-29.
//

#include "bpSolRelSet.h"

bpSolRelSet::bpSolRelSet() {

}
bpSolRelSet::bpSolRelSet(std::vector<bpSolRel> &source){
    bpSols=source;
}
void bpSolRelSet::addBpSol(bpSolRel &tmp){
    bpSols.push_back(tmp);
}
void bpSolRelSet::printMax(char side,char type,FILE *fp){
    if(fp==NULL){
        return;
    }
    else{
        std::vector<bpSolRel>::iterator it;
        for(it=bpSols.begin();it!=bpSols.end();it++){
            fprintf(fp,"%.0lf ",it->getMax(side,type));
        }
        fprintf(fp,"\n");
    }
}
void bpSolRelSet::printMaxDis(char side,char type,FILE *fp){
    if(fp==NULL){
        return;
    }
    else{
        std::vector<bpSolRel>::iterator it;
        for(it=bpSols.begin();it!=bpSols.end();it++){
            fprintf(fp,"%lld ",it->getMaxDis(side,type));
        }
        fprintf(fp,"\n");
    }
}
void bpSolRelSet::printDis2bp(char side,char type,FILE *fp){
    if(fp==NULL){
        return;
    }
    else{
        std::vector<bpSolRel>::iterator it;
        for(it=bpSols.begin();it!=bpSols.end();it++){
            fprintf(fp,"%lld ",it->getDis2bp(side,type));
        }
        fprintf(fp,"\n");
    }
}