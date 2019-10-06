//
// Created by ttbond on 19-7-8.
//

#include "fastqReads.h"



bool fastqReads::load(fstream &in) {
    if(in>>qname>>fast5Name>>seq>>pn>>qual){
        return true;
    }
    return false;
}
