
//
// Created by ttbond on 19-6-9.
//

#include <htslib/sam.h>
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <stdlib.h>
#include <iostream>
#include "omp.h"
#include "basicInfo.h"
#include "bedFile.h"
#include "bamFile.h"
#include "nanoSpeedFile.h"
#include "nanoSpeedBF.h"
#include "fastqFile.h"
#include "faFile.h"
#include "ttbond_math.h"
#include "alignManager.h"
#include "comRel.h"
using namespace std;


char refCache[350000000];

int main(int argc, char *argv[]){

/*
    // From non-b region cut duration data of reads.
    bedFile bedFromAno=bedFile();
    //Get ref regions
    bedFromAno.load("/media/ttbond/642A8ADA2A8AA91E/non-b/mirror_repeats/gff/chr1_MR.gff",NONBF);
    //Load the bam file
    bamFile bf=bamFile("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.sam");
    //From bam file get reads pos of non-b region
    vector<basicInfo> list=bf.getReadsPos(bedFromAno);
    nanoSpeedFile tmp("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/speed/tmp.dat");
    //Generate speedFile format to visualize using R
    tmp.generateSpeSpeedFile("test",list);
*/


/*
    bamFile bf=bamFile("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.SA.sam");
    set<string> readsNameSet=bf.getReadsSet();
    //Load speedFile format file
    nanoSpeedFile tmp("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/speed/tmp.dat");
    //Use staMethod1 to find imnormal regions of the reads and get pos on reads
    tmp.generateBedFile(&nanoSpeedFile::staMethod1,readsNameSet);
    const char *bedFileName="/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/relSpeed/imnormalRegion.bed";
    //From pos on reads to pos on ref, and output the bed style file
    vector<basicInfo> tmpBf=bf.generateRefBedFromReadBed(tmp.readName2readBed,bedFileName);

    //nanoSpeedBF tmpNanoSpeedBF;
    //tmpNanoSpeedBF.loadNanoSpeedBF(tmpBf);
    //tmpNanoSpeedBF.staRepRegions();

*/

/*
    //bamFile bf=bamFile("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.sam");
    bamFile bf=bamFile("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/NA12878_chr1_p2_0_2_SA_38.sam");
    bf.simpleQuery(argv[1]);

    //bf.selectSA("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.SA.sam");
*/

/*
    //select fastq
    fastqFile tmp(argv[1]);
    tmp.generateSubFile("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fastq/selectedFastq/selectedReads.fastq","/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/relSpeed/imnormalRegion.bedNoSoftMatch.readsName");
*/

/*
    //classify the fastq file into fail and pass two files
    fastqFile tmp(argv[1]);
    tmp.classifyFailAndPass1D();
*/

/*
    //generateDotplot
    //system("samtools view -h /media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_pacbio/HG001.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.trio.bam 1:22902056-22915321 \
    //>/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_pacbio/tmpDotplot/61/tmp.sam");
    //bamFile bf=bamFile("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp3.sam");
    //argv[1] is the bamFile; argv[2] is the doplotDir; argv[3] is the refFile
    //argv[4] is qname of the read
    bamFile bf=bamFile(argv[1]);
    puts(argv[4]);
    bf.generateDotplot(argv[2],argv[3],refCache,argv[4]);
*/

/*
    //to select specific regions but samtools itself can do the identic
    string workDir("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_pacbio");
    bamFile bf=bamFile("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_pacbio/HG001.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.trio.bam");
    string tmpBamName=workDir+string("/tmp.bam");
    vector<basicInfo>tmp;
    tmp.push_back(basicInfo(1,17955227,17956556));

    bf.selectByRefRange(tmp,tmpBamName.c_str());
    puts("aaaa");
    bamFile bf2=bamFile(tmpBamName.c_str());
    puts("bbbb");
    string dotplotPName=workDir+string("/tmpDotplot");
    string refFileName=string("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fastq/selectedFastq/human_g1k_v37.fasta");
    bf2.generateDotplot(dotplotPName.c_str(),refFileName.c_str(),refCache);
    puts("cccc");
    return 0;
    system((string("rm ")+tmpBamName).c_str());
*/

/* generate cut ref to fastq
    faFile refFile(refCache);
    refFile.loadAgctByChr(stringChrName2int(argv[2]),argv[1]);
    FILE *fp=fopen(argv[5],"w");
    long long stP=string2ll(argv[3]),edP=string2ll(argv[4]);
    fprintf(fp,">%s_%s_%s\n",argv[2],argv[3],argv[4]);
    for(char *i=refCache+stP-1;i<refCache+edP;i++){
        fprintf(fp,"%c",*i);
    }
    fprintf(fp,"\n");
*/

/* 将在某bam文件中的qname提取出来，然后在fastq文件中筛选，并且创建新的的fastq文件
    bamFile bf("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.SA.sam");
    set<string> qnameSet=bf.getReadsSet();
    fastqFile ff("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/total.fq");
    ff.selectByName(qnameSet,"/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/total_SA.fq");
*/

/*
    fastqFile ff("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/total_SA.fq");
    ff.selectByLen(1000,"tmp.dat");
    fastqFile ff2("tmp.dat");
    set<string> rel=ff2.getQnameSet();
    FILE *fp=fopen("tmp2.dat","w");
    for(auto i : rel) {
        fprintf(fp,"%s\n",i.c_str());
    }

*/


/*
    //识别有complementary alignment的reads
    //这部分是通过k-mer的方法将比对成功的k-mer坐标提取的程序
    //首先过滤有com比对的reads
    fastqFile fastqF("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/total_SA.fq");
    fastqF.selectByLen(1000,"/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/tmp.fastq");
    fastqFile fastqF1000("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/tmp.fastq");
    set<string> fastqF1000Name=fastqF1000.getQnameSet();
    bamFile bamF("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.SA.sam");
    bamF.selectByName(fastqF1000Name,"/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sam");
    bamFile bamF2("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sam");
    bamF2.selectByComAlign("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp2.sam");
    bamFile bamF3("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp2.sam");
    bamF3.selectByRefOverlap("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp3.sam");
    bamFile bamF4("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp3.sam");
    set<string> qnameSet=bamF4.getReadsSet();
    fastqF1000.selectByName(qnameSet,"/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/tmp2.fastq");

*/
/*
    //做大模式的筛选
    //argv[1] is fq; argv[2] is bamFile; argv[3] is refFq
    fastqFile fqF(argv[1]);
    set<string> nameSet;
    nameSet.insert(string(argv[4]));
    auto readRel=fqF.selectByName(nameSet);
    int siz=readRel.size();
    for(int i=0;i<siz;i++){
        alignManager man(readRel[i],argv[3],argv[2],refCache,true);
    }

*/

/*
    //在所有的bam中挑选长度大于1000的
    //识别有complementary alignment的reads
    //将有双向比对的reads去掉

    bamFile bamF5("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.sam");
    bamF5.selectByReadLen(1000,"/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.len1000.sam");
    bamFile bamF6("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.len1000.sam");
    fastqFile fastqF("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/total_SA.fq");
    fastqF.selectByLen(1000,"/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/tmp.fastq");
    fastqFile fastqF1000("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/tmp.fastq");
    set<string> fastqF1000Name=fastqF1000.getQnameSet();
    bamFile bamF("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.SA.sam");
    bamF.selectByName(fastqF1000Name,"/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sam");
    bamFile bamF2("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sam");
    bamF2.selectByComAlign("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp2.sam");
    bamFile bamF3("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp2.sam");
    set<string> qnameSet=bamF3.getReadsSet();
    bamF6.selectByName(qnameSet,"/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/bam/tmp.sorted.len1000.normal.sam",false);
    return 0;
    fastqF1000.selectByName(qnameSet,"/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/tmp2.fastq");

*/

/*
    //对comRel结果进行批量作图评估
    //argv[1] is fq; argv[2] is bamFile; argv[3] is refFq
    comRel rel("/home/ttbond/test/comRel.dat");
    vector<comRel::oneRel> selectedRel=rel.selectByDis(-100,0);
    set<string> selectedRelQname;
    for(auto i:selectedRel){
        selectedRelQname.insert(i.qname);
    }
    fastqFile fqF("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/fast5/tmp2.fastq");
    vector<oneFastq> fqVec=fqF.selectByName(selectedRelQname);
    for(auto i:fqVec){
        //if(i.qname==string("0a44741a-ebdd-41d0-95f3-ac4371644166"))
        alignManager manager(i,argv[3],argv[2],refCache,true);
    }
*/

///*
    //进行reads过滤的主程序
    //argv[1] is fq; argv[2] is bamFile; argv[3] is refFq ; argv[4] outFile; argv[5] core num; argv[6] qname file
    //FILE *fp=fopen("/media/ttbond/642A8ADA2A8AA91E/fast5Test/NA12878_r9.4/vcf/NA12878_rel6.sorted.bam.inv.vcf.qname","r");
    FILE *fp=fopen(argv[6],"r");
    set<string> qnameSet;
    puts("111111");
    char cache[1000];
    while(fscanf(fp,"%s",cache)!=EOF){
        qnameSet.insert(string(cache));
    }
    puts("aaaaa");
    fastqFile fqF(argv[1]);
    auto fqVec=fqF.selectByName(qnameSet);
    puts("bbbb");
    printf("fqVecSize:%d\n",fqVec.size());
    omp_set_num_threads(atoi(argv[5]));
    #pragma omp parallel for
    for(int i=0;i<fqVec.size();i++){
        //if(i.qname==string("a0bfbe15-7ff5-4bf0-aaa5-f33a8ae6843f"))
        alignManager manager(fqVec[i],argv[3],argv[2],refCache,argv[4]);
    }

//*/
/*
    omp_set_num_threads(4);
    #pragma omp parallel for
    for (int i = 0; i < 12; i++)
        printf("i = %d, I am Thread %d\n", i, omp_get_thread_num());
*/


/*
    //使用qname文件筛选bam，产生sub bam文件
    //bamFilterFromQname
    //argv[1] totablBamFile; argv[2] outBamFile; argv[3] qnameFile
    bamFile totalBam(argv[1]);
    FILE *fp=fopen(argv[3],"r");
    char cache[1000];
    set<string> selectedQname;
    while(fscanf(fp,"%s",cache)!=EOF){
        selectedQname.insert(string(cache));
    }
    fclose(fp);
    totalBam.selectByName(selectedQname,argv[2]);
*/

/*
    //使用qname文件筛选fastq，产生sub fastq文件
    //fastqFilterFromQname
    //argv[1] totablFastqFile; argv[2] outFastqFile; argv[3] qnameFile
    fastqFile totalFq(argv[1]);
    FILE *fp=fopen(argv[3],"r");
    char cache[1000];
    set<string> selectedQname;
    while(fscanf(fp,"%s",cache)!=EOF){
        selectedQname.insert(string(cache));
    }
    fclose(fp);
    totalFq.selectByName(selectedQname,argv[2]);
*/

/*
    //argv[1] is fq; argv[2] is bamFile; argv[3] is refFq; argv[4] outDir; argv[5] core num; argv[6] qname file; argv[7] workDir
    FILE *fp=fopen(argv[6],"r");
    set<string> qnameSet;
    char cache[1000];
    while(fscanf(fp,"%s",cache)!=EOF){
        qnameSet.insert(string(cache));
    }
    fastqFile fqF(argv[1]);
    auto fqVec=fqF.selectByName(qnameSet);
    printf("fqVecSize:%d\n",fqVec.size());
    //omp_set_num_threads(atoi(argv[5]));
    //#pragma omp parallel for
    for(int i=0;i<fqVec.size();i++){
        alignManager manager(fqVec[i]);
        manager.getDotplot(argv[3],argv[2],refCache,argv[4],argv[7]);
    }
*/

/*
    //在bamfile中产生qname文件
    //bam2qname
    //argv[1] bamFile; argv[2] qnameFile;
    FILE *fp=fopen(argv[2],"w");
    bamFile bamF(argv[1]);
    puts("abcdefg");
    auto qnameSet=bamF.getQnameSet();
    for(auto i:qnameSet){
        fprintf(fp,"%s\n",i.c_str());
    }
    fclose(fp);
*/
/*
    //筛选符合model1的潜在qname
    //bamModel1Qname
    //argv[1] is bamFile; argv[2] outFile
    bamFile bamf(argv[1]);
    bamf.selectByModel1(argv[2]);
*/
/*
    //在某条染色体上随机生成在非gap区域的位点
    //argv[1] refFile; argv[2] chr; argv[3] outFile; argv[4] posNum
    //ranPosFromRef
    faFile faF(refCache);
    faF.getValidatedRegionByChr(atoi(argv[2]),argv[1]);
    vector<basicInfo> ranPos;
    faF.getRandomPos(atoi(argv[4]),ranPos);
    FILE *fp=fopen(argv[3],"w");
    for(auto it:ranPos){
        fprintf(fp,"1\t%lld\t%lld\n",it.st-250,it.ed+250);
    }
    fclose(fp);
*/

}
