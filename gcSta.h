//
// Created by ttbond on 18-10-18.
//

#ifndef SRCFIND_GCSTA_H
#define SRCFIND_GCSTA_H


class gcSta {
public:
    gcSta();
    void addStr(char *str,long long len,long long num=1);
    double getGcContent();
    long long numa,numg,numc,numt;
    double gcContent;
};


#endif //SRCFIND_GCSTA_H
