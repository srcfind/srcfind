//
// Created by ttbond on 19-9-4.
//

#ifndef SRCFIND_COMREL_H
#define SRCFIND_COMREL_H

#include<string>
#include<vector>
using namespace std;

class comRel {
public:
    class oneRel{
    public:
        oneRel(string qname,double dis);
        oneRel();
        string qname;
        double dis;
    };
    comRel(const char *fileName);
    vector<oneRel> selectByDis(double low,double high);
    string fName;
};


#endif //SRCFIND_COMREL_H
