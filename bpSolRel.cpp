//
// Created by ttbond on 18-9-29.
//

#include "bpSolRel.h"

bpSolRel::bpSolRel(bpSideSolRel &_left, bpSideSolRel &_right) {
    left= new bpSideSolRel(_left);
    right= new bpSideSolRel(_right);
}

long long bpSolRel::getDis2bp(char side,char type){
    bpSideSolRel *tmp;
    switch(side){
        case 'L':
            tmp=left;break;
        case 'R':
            tmp=right;break;
        default:
            printf("Wrong director in bpSolRel::getDis2bp\n");
            exit(1);
    }
    switch(type){
        case 'R':
            return tmp->dis2bpRev;
        case 'D':
            return tmp->dis2bpDir;
        case 'M':
            return tmp->dis2bpMir;
        default:
            printf("Wrong director in bpSolRel::getDis2bp\n");
            exit(1);
    }
}
long long bpSolRel::getMaxDis(char side,char type){
    bpSideSolRel *tmp;
    switch(side){
        case 'L':
            tmp=left;break;
        case 'R':
            tmp=right;break;
        default:
            printf("Wrong director in bpSolRel::getDis2bp\n");
            exit(1);
    }
    switch(type){
        case 'R':
            return tmp->maxDisRev;
        case 'D':
            return tmp->maxDisDir;
        case 'M':
            return tmp->maxDisMir;
        default:
            printf("Wrong director in bpSolRel::getDis2bp\n");
            exit(1);
    }
}
double bpSolRel::getMax(char side,char type){
    bpSideSolRel *tmp;
    switch(side){
        case 'L':
            tmp=left;break;
        case 'R':
            tmp=right;break;
        default:
            printf("Wrong director in bpSolRel::getDis2bp\n");
            exit(1);
    }
    switch(type){
        case 'R':
            return tmp->maxRevComScore;
        case 'D':
            return tmp->maxDirRepScore;
        case 'M':
            return tmp->maxMirRepScore;
        default:
            printf("Wrong director in bpSolRel::getDis2bp\n");
            exit(1);
    }
}
