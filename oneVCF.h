//
// Created by ttbond on 18-9-14.
//

#ifndef SRCFIND_TTBOND_VCF_H
#define SRCFIND_TTBOND_VCF_H

#include<stdio.h>
#include<map>
#include<string>
#include<vector>

#include"ttbond_math.h"
#include"stdlib.h"

class oneVCF{
public:

    oneVCF(char *str);
    oneVCF();
    void initFormat(char *str);



    static int knowFormat;
    static int infoNum;
    static std::vector<std::string>infoName;
    char *str;
    int chr;
    long long st;
    long long ed;
    long long length;
    svType sv;
    std::vector<std::string> info;
};



#endif //SRCFIND_TTBOND_VCF_H
