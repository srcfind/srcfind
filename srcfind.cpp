/*
#include<stdio.h>
#include<string.h>
#include<vector>
#include<stdlib.h>


#include"ttbond_fa.h"
#include"vcfFile.h"
#include"vcfFile.cpp"
#include"detectRegion.h"
#include"detectRegion.cpp"

char agct[300000000];




int main()
{
 */
    /*
    FILE *fp=fopen("/home/ttbond/桌面/From JD/201809/repeat Regions/2/_leftBp.fa","r");
    fgets(agct,1000000,fp);
    agct[strlen(agct)-1]='\0';
    FILE *mirFile=fopen("chr1MirRep.scr","w");
    FILE *dirFile=fopen("chr1DirRep.scr","w");
    FILE *revFile=fopen("chr1RevRep.scr","w");
    detectRegion region(agct);
    region.getMirrorRepeatScore(mirFile);
    region.getDirectRepeatScore(dirFile);
    region.getReverseComScore(revFile);
    fclose(mirFile);
    fclose(dirFile);
    fclose(revFile);
    return 0;
    */

/*
    long long agctLen=loadAgctByChr(1,"GRCh38.d1.vd1.fa",agct);
    vcfFile myVcfFile("NA19240.vcf");
    vcfFile chr1SV(myVcfFile.selectByChr(1));
    vcfFile chr1NormalSV(chr1SV.selectBySv(DEL,true));
    std::vector<int>chr;
    std::vector<long long>stL,stR;
    //std::vector<long long>st2;
    std::vector<long long>edL,edR;
    std::vector<svType>svs;
    //std::vector<long long>ed2;
    chr1NormalSV.getDetRegBpSides('L',500,chr,stL,edL,svs,agctLen);
    chr.clear();
    svs.clear();
    chr1NormalSV.getDetRegBpSides('R',500,chr,stR,edR,svs,agctLen);
    int regionNum=chr.size();

     //
    //
    FILE *fp=fopen("chr1.rcf","w");
    FILE *mirFile=fopen("chr1MirRep.scr","w");
    FILE *dirFile=fopen("chr1DirRep.scr","w");
    FILE *revFile=fopen("chr1RevRep.scr","w");
    double maxLen=0;
    for(int i=0;i<regionNum;i++){
        //detectRegion region(agct,chr[i],st[i],ed[i],agct,chr[i],st2[i],ed2[i]);
        detectRegion region(agct,chr[i],st[i],ed[i],svs[i]);
        region.getMirrorRepeatScore(mirFile);
        region.getDirectRepeatScore(dirFile);
        region.getReverseComScore(revFile);
        //region.printDetectRel(fp);
        //getchar();
    }
    fclose(fp);
    fclose(mirFile);
    fclose(dirFile);
    fclose(revFile);
     //
    //printf("%d\n",chr1NormalSV.vcfNum);
    //ab.print();
    //printf("%d\n",oneVCF::knowFormat);
    //char *fileName="GRCh38.d1.vd1.fa";
    //loadAgctByChr(11,fileName,agct);
    return 0;

}
*/
