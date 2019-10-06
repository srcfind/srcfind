//
// Created by ttbond on 19-7-8.
//

#ifndef SRCFIND_FASTQREADS_H
#define SRCFIND_FASTQREADS_H
#include <string>
#include <fstream>
using namespace std;


class fastqReads {
public:
    bool load(fstream &in);
    string qname,fast5Name,seq,pn,qual;
};


#endif //SRCFIND_FASTQREADS_H
