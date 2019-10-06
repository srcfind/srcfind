//
// Created by ttbond on 18-10-9.
//

#ifndef SRCFIND_GENEREGION_H
#define SRCFIND_GENEREGION_H
#include<map>
#include<vector>
#include<set>
#include<algorithm>
#include<string>
#include"detectRegion.h"
#include"solRel.h"
#include"ttbond_math.h"
#include"ttbond_fa.h"
#include"basicInfo.h"
#include"agctTree.h"


class geneRegion {
public:
    geneRegion(int init);
    geneRegion(char *source,char *_agct);
    void addRegion(char *source);
    void addRegion(int _chr,long long _st,long long _ed);
    void solAndPrintRel(FILE *fp);
    void getSeqSta(agctTree *myT);


    static int firstGeneRegion,refChr;
    static std::map<int,std::string> id2name;
    static std::set<int> allGeneId;

    int id,num;
    char *agct;
    std::vector<basicInfo> myInfo;
    std::vector<solRel>myRel;
};


#endif //SRCFIND_GENEREGION_H
