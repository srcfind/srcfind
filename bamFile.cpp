//
// Created by ttbond on 19-6-10.
//

#ifndef SRCFIND_BAMFILE_CPP
#define SRCFIND_BAMFILE_CPP

#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <htslib/sam.h>
#include <set>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "bamFile.h"
#include "alignInfo.h"
#include "ttbond_math.h"
#include "ttbond_fa.h"
using namespace std;


bamFile::bamFile(const char *fileName) {
    int siz=strlen(fileName);
    fName=new char[siz+1];
    strcpy(fName,fileName);
    nameIsMapStPos=false;
}
bamFile::~bamFile() {
    delete []fName;
}


//from ref pos get readspos
vector<basicInfo> bamFile::getReadsPos(bedFile &refRegion){
    vector<basicInfo> &bedRegion=refRegion.myRegion;
    vector<basicInfo> rel;
    sort(bedRegion.begin(),bedRegion.end());
    int bedRegionSiz=bedRegion.size();
    int bedRegionStP=0;
    char *qname=NULL;
    bam_hdr_t *header=NULL;
    bam1_t *aln=bam_init1();
    samFile *in=hts_open(fName,"r");
    uint32_t *cigar;
    header=sam_hdr_read(in);
    int cot=0;
    //遍历每一条比对好的reads
    while(sam_read1(in,header,aln)>=0){
        qname=bam_get_qname(aln);
        int32_t chr_id= aln->core.tid;
        if(chr_id<0){
            continue;
        }
        const char *chr=header->target_name[chr_id];
        if(strcmp(chr,"1")!=0){
            continue;
        }
        long long refSt=aln->core.pos+1;
        cigar=bam_get_cigar(aln);
        int matchLen=0;
        for(int i=0;i<aln->core.n_cigar;i++){
            int icigar=cigar[i];
            int op=bam_cigar_op(icigar);
            int len=bam_cigar_oplen(icigar);
            if(op==0 || op==2){
                matchLen+=len;
            }
        }
        long long refEd=refSt+matchLen;
        while(bedRegionStP<bedRegionSiz && bedRegion[bedRegionStP].st<refSt){
            bedRegionStP++;
        }
        for(int i=bedRegionStP;i<bedRegionSiz;i++){
            if(bedRegion[i].st>refEd){
                break;
            }
            else if(bedRegion[i].ed>refEd){
                continue;
            }
            long long offSam=0,offSt=0,samSt=-1,samEd=-1;
            //printf("bed: %lld %lld\n",bedRegion[i].st,bedRegion[i].ed);
            //printf("ref: %lld %lld\n",refSt,refEd);
            for(int j=0;j<aln->core.n_cigar;j++){
                int icigar=cigar[j];
                int op=bam_cigar_op(icigar);
                int len=bam_cigar_oplen(icigar);
                if(op==0 || op==2){
                    offSt+=len;
                }
                if(op!=2){
                    offSam+=len;
                }
                if(samSt==-1 && refSt+offSt>=bedRegion[i].st){
                    if(op==0){
                        samSt=offSam-(refSt+offSt-bedRegion[i].st);
                    }
                    else if(op==2){
                        samSt=offSam;
                    }
                }
                if(refSt+offSt>=bedRegion[i].ed){
                    if(op==0){
                        samEd=offSam-(refSt+offSt-bedRegion[i].ed);
                    }
                    else if(op==2){
                        samEd=offSam;
                    }
                    break;
                }
            }
            rel.push_back(basicInfo(1,samSt,samEd,qname));
        }
    }
    hts_close(in);
    return rel;
}

vector<basicInfo> bamFile::getRefPos(bam1_t *aln,const vector<basicInfo> &readBed,vector<bool> &matchNoSoftArea) {
    vector<basicInfo> rel;
    //It should promised that readBed is sorted and no overlap
    int32_t chr = aln->core.tid;
    long long refSt = aln->core.pos + 1;
    uint32_t *cigar = bam_get_cigar(aln);
    long long refOff=refSt,preRefOff=refSt;
    long long readOff=0,preReadOff=0;
    int readBedSiz=readBed.size(),nowReadBedP=0;
    for(int i=0;i<readBedSiz;i++){
        rel.push_back(basicInfo());
    }
    int op,len;
    for(int j=0;j<aln->core.n_cigar;j++) {
        int icigar = cigar[j];
        op = bam_cigar_op(icigar);
        len = bam_cigar_oplen(icigar);
        switch(op){
            case 0:
                refOff+=len;
                readOff+=len;
                break;
            case 1:
                readOff+=len;
                break;
            case 2:
                refOff+=len;
                break;
            case 4:
                readOff+=len;
                break;
            case 5:
                readOff+=len;
                break;
            default:
                puts("innormal cigar!");
                exit(0);
        }
        for(int i=nowReadBedP;i<readBedSiz;i++){
            long long residual=0;
            if(-1==rel[i].st && readBed[i].st<=readOff){
                switch(op){
                    case 0:
                        residual=readOff-readBed[i].st;
                        rel[i].st=refOff-residual;
                        break;
                    case 1:
                        rel[i].st=refOff;
                        break;
                    case 2:
                        break;
                    case 4:
                        rel[i].st=-2;
                        break;
                    case 5:
                        rel[i].st=-2;
                        break;
                }
            }
            if(-1==rel[i].ed && readBed[i].ed<=readOff){
                switch(op){
                    case 0:
                        residual=readOff-readBed[i].ed;
                        rel[i].ed=refOff-residual;
                        break;
                    case 1:
                        rel[i].ed=refOff;
                    case 2:
                        break;
                    case 4:
                        rel[i].ed=-2;
                        break;
                    case 5:
                        rel[i].ed=-2;
                        break;
                }
                if(-2!=rel[i].st && -2!=rel[i].ed){
                    matchNoSoftArea[i]=true;
                }
                nowReadBedP=i+1;
            }
            if(readBed[i].st>readOff){
                break;
            }
        }
    }
    return rel;
}
vector<basicInfo> bamFile::getRefPos(bam1_t *aln,const vector<basicInfo> &readBed) {
    vector<basicInfo> rel;
    //It should promised that readBed is sorted and no overlap
    int32_t chr = aln->core.tid;
    long long refSt = aln->core.pos + 1;
    uint32_t *cigar = bam_get_cigar(aln);
    long long refOff=refSt,preRefOff=refSt;
    long long readOff=0,preReadOff=0;
    int readBedSiz=readBed.size(),nowReadBedP=0;
    for(int i=0;i<readBedSiz;i++){
        rel.push_back(basicInfo());
    }
    int op,len;
    for(int j=0;j<aln->core.n_cigar;j++) {
        int icigar = cigar[j];
        op = bam_cigar_op(icigar);
        len = bam_cigar_oplen(icigar);
        switch(op){
            case 0:
                refOff+=len;
                readOff+=len;
                break;
            case 1:
                readOff+=len;
                break;
            case 2:
                refOff+=len;
                break;
            case 4:
                readOff+=len;
                break;
            case 5:
                readOff+=len;
                break;
            default:
                puts("innormal cigar!");
                exit(0);
        }
        for(int i=nowReadBedP;i<readBedSiz;i++){
            long long residual=0;
            if(-1==rel[i].st && readBed[i].st<=readOff){
                switch(op){
                    case 0:
                        residual=readOff-readBed[i].st;
                        rel[i].st=refOff-residual;
                        break;
                    case 1:
                        rel[i].st=refOff;
                        break;
                    case 2:
                        break;
                    case 4:
                        rel[i].st=-2;
                        break;
                    case 5:
                        rel[i].st=-2;
                        break;
                }
            }
            if(-1==rel[i].ed && readBed[i].ed<=readOff){
                switch(op){
                    case 0:
                        residual=readOff-readBed[i].ed;
                        rel[i].ed=refOff-residual;
                        break;
                    case 1:
                        rel[i].ed=refOff;
                    case 2:
                        break;
                    case 4:
                        rel[i].ed=-2;
                        break;
                    case 5:
                        rel[i].ed=-2;
                        break;
                }
                nowReadBedP=i+1;
            }
            if(readBed[i].st>readOff){
                break;
            }
        }
    }
    return rel;
}

vector<basicInfo> bamFile::generateRefBedFromReadBed(map<string,vector<basicInfo> > &readBed,const char *bedFileName){
    vector<basicInfo> rel;
    map<string,vector<bool> > machNoSoftArea;
    for(auto it=readBed.begin();it!=readBed.end();it++){
        vector<bool> tmp((it->second).size(),false);
        machNoSoftArea[it->first]=tmp;
    }
    int bedRegionStP=0;
    char *qname=NULL;
    bam_hdr_t *header=NULL;
    bam1_t *aln=bam_init1();
    samFile *in=hts_open(fName,"r");
    uint32_t *cigar;
    header=sam_hdr_read(in);
    int cot=0;
    //遍历每一条比对好的reads
    while(sam_read1(in,header,aln)>=0) {
        if(cot%100==0){
            printf("%d bam align have been processed!\n",cot);
        }
        cot++;
        qname = bam_get_qname(aln);
        if(aln->core.tid < 0){
            continue;
        }
        if(readBed.find(string(qname))!=readBed.end()) {
            vector<basicInfo> tmp=getRefPos(aln,readBed[string(qname)],machNoSoftArea[string(qname)]);
            for(basicInfo &info:tmp){
                info.cache=new char[strlen(qname)+5];
                info.chr=stringChrName2int(header->target_name[aln->core.tid]);
                strcpy(info.cache,qname);
            }
            rel.insert(rel.end(),tmp.begin(),tmp.end());
        }
    }
    hts_close(in);
    FILE *fp=fopen(bedFileName,"w");
    char *machNoSoftAreaBedFileName=new char[strlen(bedFileName)+20];
    strcpy(machNoSoftAreaBedFileName,bedFileName);
    strcat(machNoSoftAreaBedFileName,"NoSoftMatch");
    FILE *fp2=fopen(machNoSoftAreaBedFileName,"w");
    for(auto it=machNoSoftArea.begin();it!=machNoSoftArea.end();it++){
        vector<bool> &tmpVector=(it->second);
        vector<basicInfo> &tmpBasicInfoV=readBed[(it->first)];
        int siz=tmpVector.size();
        for(int i=0;i<siz;i++){
            if(true||(!tmpVector[i])){
                tmpBasicInfoV[i].printMe(fp2);
                fprintf(fp2,"%s\n",(it->first).c_str());
            }
        }
    }
    for(const basicInfo &info:rel){
        fprintf(fp,"%d\t%lld\t%lld\t%s\n",info.chr,info.st,info.ed,info.cache);
    }
    fclose(fp);
    fclose(fp2);
    delete []machNoSoftAreaBedFileName;
    return rel;
}

void bamFile::selectSA(const char *outBamFileName) {
    set<string> SAReadName;
    set<string> allReadName;

    {
        char *qname = NULL;
        bam_hdr_t *header = NULL;
        bam1_t *aln = bam_init1();
        samFile *in = hts_open(fName, "r");
        uint32_t *cigar;
        header = sam_hdr_read(in);
        while (sam_read1(in, header, aln) >= 0) {
            qname = bam_get_qname(aln);
            allReadName.insert(string(qname));
            if ((aln->core.flag) & BAM_FSUPPLEMENTARY) {
                SAReadName.insert(string(qname));
            }
        }
    }
    {
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        samFile *out = hts_open(outBamFileName, "w");
        while (sam_read1(in2, header2, aln2) >= 0) {
            qname = bam_get_qname(aln2);
            if (SAReadName.find(string(qname)) != SAReadName.end()) {
                sam_write1(out, header2, aln2);
            }
        }
        hts_close(out);
    }
    printf("Total %d reads in the bam file\n",allReadName.size());
    printf("and %d(%.2lf) SA reads were found\n",SAReadName.size(),SAReadName.size()*1.0/allReadName.size());


}

void bamFile::selectByName(set<string> qnameSet,const char *outBamFileName,bool selected) {
    string tmpName(outBamFileName);
    int num=0;
    {
        tmpName=tmpName+".tmp";
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        samFile *out = hts_open(tmpName.c_str(), "w");
        while (sam_read1(in2, header2, aln2) >= 0) {
            num++;
            if(num%10000==0){
                printf("num:%d\n",num);
            }
            qname = bam_get_qname(aln2);
            if (( qnameSet.find(string(qname)) != qnameSet.end() )==selected) {
                sam_write1(out, header2, aln2);
            }
        }
        hts_close(out);
        hts_close(in2);
    }
    string com_samheader=string("samtools view -H ")+fName+string(" >")+outBamFileName;
    system(com_samheader.c_str());
    string com_alignAdd=string("cat ")+tmpName+string(" >>")+string(outBamFileName);
    system(com_alignAdd.c_str());
    string com_rmTmp=string("rm ")+tmpName;
    system(com_rmTmp.c_str());
    //printf("Total %d reads in the bam file\n",allReadName.size());
    //printf("and %d(%.2lf) SA reads were found\n",SAReadName.size(),SAReadName.size()*1.0/allReadName.size());
}

set<string> bamFile::getReadsSet() {
    char *qname=NULL;
    bam_hdr_t *header2 = NULL;
    bam1_t *aln2 = bam_init1();
    samFile *in2 = hts_open(fName, "r");
    header2 = sam_hdr_read(in2);
    set<string> readsNameSet;
    while (sam_read1(in2, header2, aln2) >= 0) {
        qname = bam_get_qname(aln2);
        readsNameSet.insert(string(qname));
    }
    hts_close(in2);
    return readsNameSet;
}
vector<basicInfo> bamFile::simpleQueryBasic(const char *qname) {
    vector<basicInfo> totalRel;
    char readName[100];
    if(qname==NULL) {
        scanf("%s", readName);
    }
    else{
        strcpy(readName,qname);
        puts(readName);
    }
    do{
        if(strcmp(readName,"q")==0){
            break;
        }
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        set<string> readsNameSet;
        long long readSt,readEd;
        while (sam_read1(in2, header2, aln2) >= 0) {
            qname=bam_get_qname(aln2);
            auto hardLen=getHardClipLen(aln2);
            int totLen=hardLen+aln2->core.l_qseq;
            if(strcmp(qname,readName)==0){
                uint32_t *cigar;
                cigar=bam_get_cigar(aln2);
                long long readMatchLen=0;
                for(int i=0;i<aln2->core.n_cigar;i++){
                    int icigar=cigar[i];
                    int op=bam_cigar_op(icigar);
                    int len=bam_cigar_oplen(icigar);
                    if(i==0){
                        readSt=len;
                    }
                    if(i==aln2->core.n_cigar-1){
                        readEd=readMatchLen;
                    }
                    if(op!=2){
                        readMatchLen+=len;
                    }
                }
                vector<basicInfo>rel,tmp;
                tmp.push_back(basicInfo(0,readSt+1,readEd));
                rel=getRefPos(aln2,tmp);
                if(((aln2->core.flag)&16)){
                    long long tmpReadSt=readSt;
                    readSt=totLen-readEd;
                    readEd=totLen-tmpReadSt;

                }
                int chr_id=aln2->core.tid;
                const char *chr;
                if(chr_id<0){
                    chr="-1";
                }
                else {
                    chr = header2->target_name[chr_id];
                }
                totalRel.push_back(basicInfo(stringChrName2int(chr),rel[0].st,rel[0].ed));
                printf("ref chr:%s\t st:%lld\ted:%lld\t",chr,rel[0].st,rel[0].ed);
                printf("read st:%d\t",readSt+1);
                printf("ed:%lld\ttotalLen:%lld\tflag:%d\n",readEd,totLen,aln2->core.flag);
            }
        }
        hts_close(in2);
    }while(qname==NULL && scanf("%s",readName));
    return totalRel;
}
vector<alignInfo> bamFile::simpleQuery(const char *qname) {
    vector<alignInfo> totalRel;
    char readName[100];
    if(qname==NULL) {
        scanf("%s", readName);
    }
    else{
        strcpy(readName,qname);
        puts(readName);
    }
    do{
        if(strcmp(readName,"q")==0){
            break;
        }
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        set<string> readsNameSet;
        long long readSt,readEd;
        long long readLen=0;
        while (sam_read1(in2, header2, aln2) >= 0) {
            qname=bam_get_qname(aln2);
            if(strcmp(qname,readName)==0){
                uint32_t *cigar;
                cigar=bam_get_cigar(aln2);
                long long readMatchLen=0;
                for(int i=0;i<aln2->core.n_cigar;i++){
                    int icigar=cigar[i];
                    int op=bam_cigar_op(icigar);
                    int len=bam_cigar_oplen(icigar);
                    if(i==0){
                        readSt=len;
                    }
                    if(i==aln2->core.n_cigar-1){
                        readEd=readMatchLen;
                    }
                    if(op!=2){
                        readMatchLen+=len;
                    }
                }
                vector<basicInfo>rel,tmp;
                tmp.push_back(basicInfo(0,readSt+1,readEd));
                rel=getRefPos(aln2,tmp);
                if(((aln2->core.flag)&16)){
                    long long tmpReadSt=readSt;
                    readSt=readLen-readEd;
                    readEd=readLen-tmpReadSt;

                }
                int chr_id=aln2->core.tid;
                const char *chr;
                if(chr_id<0){
                    chr="-1";
                }
                else {
                    chr = header2->target_name[chr_id];
                }
                totalRel.push_back(alignInfo(stringChrName2int(chr),rel[0].st,rel[0].ed,readSt+1,readEd,qname,(int)(aln2->core.flag)));
                printf("ref chr:%s\t st:%lld\ted:%lld\t",chr,rel[0].st,rel[0].ed);
                printf("read st:%d\t",readSt+1);
                printf("ed:%lld\ttotalLen:%lld\tflag:%d\n",readEd,readLen,aln2->core.flag);
            }
        }
        hts_close(in2);
    }while(qname==NULL && scanf("%s",readName));
    return totalRel;
}

long long bamFile::getReadLen(bam1_t *aln) {
    long long rel=0;
    for(int i=0;i<aln->core.n_cigar;i++){
        uint32_t *cigar=bam_get_cigar(aln);
        int icigar=cigar[i];
        int op=bam_cigar_op(icigar);
        int len=bam_cigar_oplen(icigar);
        if(op!=2){
            rel+=len;
        }
    }
    return rel;
}

long long bamFile::getRefMatchLen(bam1_t *aln){
    long long rel=0;
    for(int i=0;i<aln->core.n_cigar;i++){
        uint32_t *cigar=bam_get_cigar(aln);
        int icigar=cigar[i];
        int op=bam_cigar_op(icigar);
        int len=bam_cigar_oplen(icigar);
        if(op==0||op==2||op==7||op==8){
            rel+=len;
        }
    }
    return rel;
}



void bamFile::generateDotplot(const char *outP,const char *refFile,char *refCache,char *_qname){
    //原来这个函数的设计是要自己产生readfq和reffq，但是考虑到整条reads的完整性高，画dotplot的
    //时候应该在fastq文件中提取read
    char *qname = NULL;
    bam_hdr_t *header = NULL;
    bam1_t *aln = bam_init1();
    samFile *in = hts_open(fName, "r");
    uint32_t *cigar;
    header = sam_hdr_read(in);
    char cache[100];
    char command[200];
    long long extendLen;
    string workDir=string(outP)+string("/tmp");
    string tmpRefBed=workDir+string("/tmpRef.bed");
    string tmpRefFasta=workDir+string("/tmpRef.fq");
    string tmpReadFasta=workDir+string("/tmpRead.fq");
    int refChr=-1;
    int readId=0;
    while (sam_read1(in, header, aln) >= 0) {
        readId++;
        string qname(bam_get_qname(aln));
        //puts(qname.c_str());
        if(qname!=string(_qname)){
            continue;
        }
        //puts(tmpRefFasta.c_str());
        int chr_id=aln->core.tid;
        long long stPos=aln->core.pos+1;
        sprintf(cache,"%d_%s_%lld",readId,header->target_name[chr_id],stPos);
        string pos(cache);
        string chrName(header->target_name[chr_id]);
        string fileName=qname+string("_")+string(cache)+string(".png");
        modifyFileName(fileName);
        long long refMatchLen=getRefMatchLen(aln);
        extendLen=max(2000.0,0.2*refMatchLen);
        long long refCutSt=max(1,stPos-extendLen);
        long long refCutEd=stPos+refMatchLen+extendLen;
        if(isRegularChr((string(">")+chrName).c_str())){
            int chr = stringChrName2int(chrName.c_str());
            if(refChr!=chr){
                long long sta=loadAgctByChr(chr,refFile,refCache);
                refChr=chr;
            }
            FILE *fp=fopen(tmpRefFasta.c_str(),"w");
            fprintf(fp,">%s\n",chrName.c_str());
            refCutEd=min((long long)strlen(refCache),refCutEd);
            for(const char *i=refCache+refCutSt-1;i<refCache+refCutEd;i++){
                fprintf(fp,"%c",(*i));
            }
            fprintf(fp,"\n");
            fclose(fp);
            //FILE *readFp=fopen(tmpReadFasta.c_str(),"w");
            //fprintf(readFp,">%s\n",qname.c_str());
            //uint8_t *seq=bam_get_seq(aln);
            //for(int i=0;i<aln->core.l_qseq;i++){
            //    fprintf(readFp,"%c",seq_nt16_str[bam_seqi(seq,i)]);
            //}
            //fprintf(readFp,"\n");
            //fclose(readFp);
            string dotplotCommand=string("java -cp /home/ttbond/softwares/gepard/dist/Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 ")+
                                  tmpRefFasta+string(" -seq2 ")+tmpReadFasta+string(" -matrix /home/ttbond/softwares/gepard/resources/matrices/edna.mat")+
                                  string(" -maxwidth 2000 -maxheight 2000 -outfile ")+string(outP)+string("/")+fileName+string(" -silent");
            system(dotplotCommand.c_str());
            //exit(0);
        }
        /*
        FILE *fp=fopen(tmpRefBedF.c_str(),"w");
        fprintf(fp,"%s\t%lld\t%lld\n",chrName.c_str(),refCutSt,refCurEd);
        */

    }
    hts_close(in);
}

void bamFile::selectByRefRange(const char *inputBedName,const char *outPutBamName){
    bedFile tmpBed(inputBedName,NULL);
    tmpBed.sortMyRegion();
    selectByRefRange(tmpBed.myRegion,outPutBamName);
}

void bamFile::selectByRefRange(vector<basicInfo> inputBed,const char *outPutBamName){
    int bedSt=0,bedSiz=inputBed.size();
    {
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        samFile *out = hts_open(outPutBamName, "w");
        while (sam_read1(in2, header2, aln2) >= 0 && bedSt<bedSiz) {
            qname = bam_get_qname(aln2);
            long long refSt=aln2->core.pos+1;
            long long refEd=refSt+getRefMatchLen(aln2);
            int chr_num=stringChrName2int(header2->target_name[aln2->core.tid]);
            if(chr_num==-1){
                continue;
            }
            basicInfo tmp(chr_num,refSt,refEd);
            basicInfo maxEd=basicInfo();
            for(int i=bedSt;i<bedSiz;i++){
                if(inputBed[i].myStGreaterThan(tmp)){
                    break;
                }
                if(inputBed[i].myEdLargerThanEdOf(maxEd)){
                    maxEd=inputBed[i];
                }
                if(inputBed[i].overlapRate(tmp)>0){
                    sam_write1(out, header2, aln2);
                }
                if(maxEd.myEdSmallerThanStOf(tmp)){
                    bedSt=i+1;
                    continue;
                }
            }
        }
        hts_close(out);
        hts_close(in2);
    }
}

void bamFile::sortByPos(){
    //can not be used and to be update
}

bool bamFile::isSorted() {
    //can not be used and to be update
    char *qname = NULL;
    bam_hdr_t *header = NULL;
    bam1_t *aln = bam_init1();
    samFile *in = hts_open(fName, "r");
    uint32_t *cigar;
    header = sam_hdr_read(in);
    return false;
}

void bamFile::modifyFileName(string &fn){
    for(char &tmp : fn){
        if('/'==tmp){
            tmp='_';
        }
    }
}

void bamFile::selectByComAlign(const char *outPutBamName) {
    set<string> posiStrandQname,negaStrandQname,bothStrandQname;
    {
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        while (sam_read1(in2, header2, aln2) >= 0) {
            qname = bam_get_qname(aln2);
            if ((aln2->core.flag) & BAM_FREVERSE) {
                negaStrandQname.insert(string(qname));
            }
            else{
                posiStrandQname.insert(string(qname));
            }
        }
        hts_close(in2);
    }
    for(auto it=posiStrandQname.begin();it!=posiStrandQname.end();it++){
        if(negaStrandQname.find(*it)!=negaStrandQname.end()){
            bothStrandQname.insert(*it);
        }
    }
    selectByName(bothStrandQname,outPutBamName);
}

void bamFile::selectByRefOverlap(const char *outPutBamName) {
    map<string,vector<basicInfo> > info;
    {
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        while (sam_read1(in2, header2, aln2) >= 0) {
            qname = bam_get_qname(aln2);
            int chr_id=aln2->core.tid;
            const char *chrName=header2->target_name[chr_id];
            long long readLen=getReadLen(aln2);
            uint32_t *cigar=bam_get_cigar(aln2);
            long long readSt=bam_cigar_oplen(cigar[0]);
            long long readEd=readLen-bam_cigar_oplen(cigar[(aln2->core.n_cigar)-1]);
            vector<basicInfo>rel,tmp;
            tmp.push_back(basicInfo(0,readSt+1,readEd));
            rel=getRefPos(aln2,tmp);
            basicInfo tmpInfo(0,rel[0].st,rel[0].ed,chrName);
            info[string(qname)].push_back(tmpInfo);
        }
        hts_close(in2);
    }
    set<string> filteredQname;
    for(auto it=info.begin();it!=info.end();it++){
        vector<basicInfo> &tmp=it->second;
        int siz=tmp.size();
        for(int i=0;i<siz;i++){
            for(int j=i+1;j<siz;j++){
                if(strcmp(tmp[i].cache,tmp[j].cache)==0){
                    if(tmp[i].overlapRate(tmp[j])>0){
                        filteredQname.insert(it->first);
                    }
                }
            }
        }
    }
    selectByName(filteredQname,outPutBamName);
}

void bamFile::selectByReadLen(int minLen,const char *outPutBamName){
    copyHeader(NULL,outPutBamName);
    string tmpName=string(outPutBamName)+".tmp";
    {
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        samFile *out = hts_open(outPutBamName, "a");
        while (sam_read1(in2, header2, aln2) >= 0) {
            qname = bam_get_qname(aln2);
            long long refSt=aln2->core.pos+1;
            long long refEd=refSt+getRefMatchLen(aln2);
            if(aln2->core.tid==-1){
                continue;
            }
            int chr_num=stringChrName2int(header2->target_name[aln2->core.tid]);
            if(chr_num==-1){
                continue;
            }
            basicInfo tmp(chr_num,refSt,refEd);
            basicInfo maxEd=basicInfo();
            if(aln2->core.l_qseq>=minLen){
                sam_write1(out, header2, aln2);
            }
        }
        hts_close(out);
        hts_close(in2);
    }
}

void bamFile::copyHeader(const char *tmpName,const char *outPutBamName) {
    string com_samheader=string("samtools view -H ")+fName+string(" >")+outPutBamName;
    system(com_samheader.c_str());
}

set<string> bamFile::getQnameSet() {
    set<string> rel;
    {
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        while (sam_read1(in2, header2, aln2) >= 0) {
            qname = bam_get_qname(aln2);
            if(aln2->core.tid==-1){
                continue;
            }
            int chr_num=stringChrName2int(header2->target_name[aln2->core.tid]);
            if(chr_num==-1){
                //continue;
            }
            rel.insert(string(qname));
        }
        hts_close(in2);
    }
    return rel;
}

long long bamFile::getHardClipLen(bam1_t *aln) {
    long long rel=0;
    uint32_t *cigar=bam_get_cigar(aln);
    for(int j=0;j<aln->core.n_cigar;j++) {
        int icigar = cigar[j];
        int op = bam_cigar_op(icigar);
        int len = bam_cigar_oplen(icigar);
        switch(op){
            case 0:
            case 1:
            case 2:
            case 4:
                break;
            case 5:
                rel+=len;
                break;
            default:
                puts("innormal cigar!");
                exit(0);
        }
    }
    return rel;
}

vector<alignInfo> bamFile::simpleQueryAlign(const char *qname, long long int readLen) {
    vector<alignInfo> totalRel;
    char readName[100];
    if(qname==NULL) {
        scanf("%s", readName);
    }
    else{
        strcpy(readName,qname);
        puts(readName);
    }
    do{
        if(strcmp(readName,"q")==0){
            break;
        }
        char *qname=NULL;
        bam_hdr_t *header2 = NULL;
        bam1_t *aln2 = bam_init1();
        samFile *in2 = hts_open(fName, "r");
        header2 = sam_hdr_read(in2);
        set<string> readsNameSet;
        long long readSt,readEd;
        while (sam_read1(in2, header2, aln2) >= 0) {
            qname=bam_get_qname(aln2);
            if(strcmp(qname,readName)==0){
                uint32_t *cigar;
                cigar=bam_get_cigar(aln2);
                long long readMatchLen=0;
                for(int i=0;i<aln2->core.n_cigar;i++){
                    int icigar=cigar[i];
                    int op=bam_cigar_op(icigar);
                    int len=bam_cigar_oplen(icigar);
                    if(i==0){
                        readSt=len;
                    }
                    if(i==aln2->core.n_cigar-1){
                        readEd=readMatchLen;
                    }
                    if(op!=2){
                        readMatchLen+=len;
                    }
                }
                vector<basicInfo>rel,tmp;
                tmp.push_back(basicInfo(0,readSt+1,readEd));
                rel=getRefPos(aln2,tmp);
                if(((aln2->core.flag)&16)){
                    long long tmpReadSt=readSt;
                    readSt=readLen-readEd;
                    readEd=readLen-tmpReadSt;

                }
                int chr_id=aln2->core.tid;
                const char *chr;
                if(chr_id<0){
                    chr="-1";
                }
                else {
                    chr = header2->target_name[chr_id];
                }
                totalRel.push_back(alignInfo(stringChrName2int(chr),rel[0].st,rel[0].ed,readSt+1,readEd,qname,(int)(aln2->core.flag)));
                printf("ref chr:%s\t st:%lld\ted:%lld\t",chr,rel[0].st,rel[0].ed);
                printf("read st:%d\t",readSt+1);
                printf("ed:%lld\ttotalLen:%lld\tflag:%d\n",readEd,readLen,aln2->core.flag);
            }
        }
        hts_close(in2);
    }while(qname==NULL && scanf("%s",readName));
    return totalRel;
}

void bamFile::selectByModel1(const char *outF) {
    vector<alignInfo> totalRel;
    char readName[100];
    *readName=0;
    char *qname=NULL;
    bam_hdr_t *header2 = NULL;
    bam1_t *aln2 = bam_init1();
    samFile *in2 = hts_open(fName, "r");
    header2 = sam_hdr_read(in2);
    set<string> readsNameSet;
    long long readSt,readEd;
    while (sam_read1(in2, header2, aln2) >= 0) {
        qname=bam_get_qname(aln2);
        if(strcmp(qname,readName)!=0) {
            //判断有没有相反方向的并且在ref上的overlap
            map<int, vector<alignInfo> > pos;
            for (auto info:totalRel) {
                if (info.refInfo.chr < 0 || info.refInfo.chr > 24) {
                    continue;
                }
                pos[info.refInfo.chr].push_back(info);
            }
            for (auto info:pos) {
                int siz = info.second.size();
                bool refOverlap = false;
                for (int i = 0; i < siz; i++) {
                    for (int j = i + 1; j < siz; j++) {
                        if ((info.second)[i].refInfo.overlapRate((info.second)[j].refInfo) > 0 &&
                        (info.second[i].flag & 16) != (info.second[j].flag & 16)) {
                            refOverlap = true;
                            break;
                        }
                    }
                }
                if (refOverlap) {
                    FILE *fp = fopen(outF, "a");
                    fprintf(fp, "%s\n", qname);
                    fclose(fp);
                }
            }
            strcpy(readName, qname);
            totalRel.clear();
        }
        uint32_t *cigar;
        cigar=bam_get_cigar(aln2);
        long long readMatchLen=0;
        for(int i=0;i<aln2->core.n_cigar;i++){
            int icigar=cigar[i];
            int op=bam_cigar_op(icigar);
            int len=bam_cigar_oplen(icigar);
            if(i==0){
                readSt=len;
            }
            if(i==aln2->core.n_cigar-1){
                readEd=readMatchLen;
            }
            if(op!=2){
                readMatchLen+=len;
            }
        }
        vector<basicInfo>rel,tmp;
        tmp.push_back(basicInfo(0,readSt+1,readEd));
        rel=getRefPos(aln2,tmp);
        if(((aln2->core.flag)&16)){
            long long tmpReadSt=readSt;
            readSt=-readEd;
            readEd=-tmpReadSt;
        }
        int chr_id=aln2->core.tid;
        const char *chr;
        if(chr_id<0){
            chr="-1";
        }
        else {
            chr = header2->target_name[chr_id];
        }
        totalRel.push_back(alignInfo(stringChrName2int(chr),rel[0].st,rel[0].ed,readSt+1,readEd,qname,(int)(aln2->core.flag)));
        //printf("ref chr:%s\t st:%lld\ted:%lld\t",chr,rel[0].st,rel[0].ed);
        //printf("read st:%d\t",readSt+1);
        //printf("ed:%lld\ttotalLen:%lld\tflag:%d\n",readEd,0,aln2->core.flag);
    }
    hts_close(in2);
}


#endif