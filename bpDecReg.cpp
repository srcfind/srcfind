//
// Created by ttbond on 18-9-28.
//

#include "bpDecReg.h"
bpDecReg::bpDecReg(char *str,int _chr,long long _st,long long _ed,svType _mySv,long long _bpPos,bpDecReg *_pairReg):
                detectRegion(str,_chr,_st,_ed,_mySv){
    bpPos=_bpPos;
    pairReg=_pairReg;
    if(pairReg!=NULL){
        pairReg->pairReg=this;
    }
}