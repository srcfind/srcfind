//
// Created by ttbond on 18-10-9.
//

#ifndef SRCFIND_GENEREIONFILE_H
#define SRCFIND_GENEREIONFILE_H
#include<vector>
#include<map>
#include"geneRegion.h"
#include"agctTree.h"


class geneRegionFile {
public:
    geneRegionFile(char *file,char *agct);
    void solAndPrintRel();
    void getSeqSta();
    ~geneRegionFile();
    std::vector<geneRegion> myGeneRegion;
    std::map<int,int> id2ind;
    agctTree *myT;
};


#endif //SRCFIND_GENEREIONFILE_H
