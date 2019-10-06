//
// Created by ttbond on 18-9-15.
//

#include "oneVCF.h"

int oneVCF::knowFormat=-1;
int oneVCF::infoNum=0;
std::vector<std::string> oneVCF::infoName;

oneVCF::oneVCF(char *str) {
    if(knowFormat==-1){
         initFormat(str);
         return;
    }
    char *j=str,*strSt=str;
    for(int i=0;i<infoNum;i++){
        //error will happen if the info is not formal
        for(;(*j)!='\t'&&(*j)!='\n';j++){
            if((*j)=='\0'){
                printf("format is not formal\n");
                puts(strSt);
                exit(-1);
            }
        };
        (*j)='\0';
        info.push_back(std::string(strSt));
        if(infoName[i]==std::string("CHROM")){
            chr=stringChrName2int(strSt);
        }
        else if(infoName[i]==std::string("POS")){
            st=string2ll(strSt);
        }
        else if(infoName[i]==std::string("ALT")){
            if(strlen(strSt)==5){
                if(strcmp(strSt,"<DEL>")==0){
                    sv=DEL;
                }
                else if(strcmp(strSt,"<DUP>")==0){
                    sv=DUP;
                }
                else if(strcmp(strSt,"<INV>")==0){
                    sv=INV;
                }
                else if(strcmp(strSt,"<INS>")==0){
                    sv=INS;
                }
                else{
                    sv=OTHER;
                }
            }
        }
        else if(infoName[i]==std::string("INFO")){
            char cache[100];
            if(selectFromSepStr(cache,strSt,"END",';')){
                if((*cache)){
                    ed=string2ll(cache);
                    length=ed-st+1;
                }
                else{
                    puts("no end info");
                    exit(-1);
                }
            }
            else{
                puts("no end info");
                exit(-1);
            }
        }
        strSt=++j;
    }
    /*
    //test weather load successfully
    for(int i=0;i<infoNum;i++){
        printf("%s\n",info[i].c_str());
    }
     */
}

oneVCF::oneVCF(){

}

void oneVCF::initFormat(char *str) {
    FILE *fp=fopen(str,"r");
    char cache[10000];
    while(fgets(cache,10000,fp)!=NULL){
        if(cache[0]==cache[1]&&cache[0]=='#'){
            continue;
        }
        else if(cache[0]=='#'&&cache[1]!=cache[0]){
            char *j=cache+1;
            knowFormat=1;
            for(char *i=cache+1;(*i)!='\0';i++){
                if((*i)=='\t'||(*i)=='\n'){
                    (*i)='\0';
                    infoName.push_back(std::string(j));
                    j=i+1;
                }
            }
            infoNum=infoName.size();
            break;
        }
        else{
            printf("Can't build format of vcf.\n");
            fclose(fp);
            exit(-1);
        }
    }
    fclose(fp);
}

