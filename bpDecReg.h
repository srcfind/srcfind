//
// Created by ttbond on 18-9-28.
//

#ifndef SRCFIND_BPDECREG_H
#define SRCFIND_BPDECREG_H
#include<stdlib.h>
#include"detectRegion.h"

class bpDecReg: public detectRegion{
public:
    bpDecReg(char *str,int _chr,long long _st,long long _ed,svType _mySv,long long _bpPos,bpDecReg *_pairReg=NULL);
    bpDecReg *pairReg;
    long long bpPos;
};


#endif //SRCFIND_BPDECREG_H
