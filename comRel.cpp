//
// Created by ttbond on 19-9-4.
//

#include "comRel.h"

comRel::comRel(const char *fileName) {
    fName=string(fileName);
}

vector<comRel::oneRel> comRel::selectByDis(double low, double high) {
    FILE *fp=fopen(fName.c_str(),"r");
    char cache[1000];
    double dis;
    vector<comRel::oneRel> rel;
    while(fscanf(fp,"%s%lf",cache,&dis)!=EOF){
        if(dis>=low&&dis<=high) {
            rel.push_back(oneRel(string(cache), dis));
        }
    }
    return rel;
}


comRel::oneRel::oneRel(string _qname, double _dis) {
    qname=_qname;
    dis=_dis;
}

comRel::oneRel::oneRel() {
    dis=-1.11;
    return;
}
