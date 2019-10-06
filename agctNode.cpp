//
// Created by ttbond on 18-10-17.
//
#ifndef SRCFIND_AGCTNODE_CPP
#define SRCFIND_AGCTNODE_CPP
#include<stdio.h>
#include<stdlib.h>
#include "agctNode.h"
agctNode::agctNode(){
    num=0;
    for(int i=0;i<4;i++){
        son[i]=NULL;
    }
}



#endif