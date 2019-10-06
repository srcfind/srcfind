//
// Created by ttbond on 18-11-15.
//

#ifndef SRCFIND_NONB_H
#define SRCFIND_NONB_H
#include<vector>
#include<utility>
#include<algorithm>
#include<string.h>
#include"basicInfo.h"
#include"ttbond_math.h"

class nonB {
public:
    nonB(std::vector<basicInfo> _myRegion,char *agct,nonBTyp _myNonBTyp);
    int getMisNum();
    bool isBrinkDiverge();
    bool operator ==(const nonB &source);
    bool operator <(const nonB &source);
    void printMe(FILE *relFp);
    nonBTyp myNonBTyp;
    int regionNum;
    std::vector<basicInfo> myRegion;
    std::vector<basicInfo> misMachPoint;
    std::vector<std::pair<int,int> > relativePos;
    std::vector<std::vector<char> > misMachAgct;


};


#endif //SRCFIND_NONB_H
