cmake_minimum_required(VERSION 3.12)
project(srcfind)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(LINK_DIR /home/xutun/miniconda3/lib)
set(Rcpp_DIR /home/xutun/miniconda3/envs/cppEnv/lib/R/library/Rcpp)
set(R_DIR /home/xutun/miniconda3/envs/cppEnv/lib/R)
set(RInside_DIR /home/xutun/miniconda3/envs/cppEnv/lib/R/library/RInside)
include_directories(.)
include_directories(${Rcpp_DIR}/include)
include_directories(${R_DIR}/include)
include_directories(${RInside_DIR}/include)
link_directories(${LINK_DIR})
link_directories(${Rcpp_DIR}/libs)
link_directories(${R_DIR}/lib)
link_directories(${RInside_DIR}/lib)
LINK_LIBRARIES(hts)
LINK_LIBRARIES(R)
LINK_LIBRARIES(RInside)
LINK_LIBRARIES(Rcpp)

add_executable(srcfind
        agctNode.cpp
        agctNode.h
        agctTree.cpp
        agctTree.h
        backgroundGenerator.cpp
        backgroundGenerator.h
        bamFile.cpp
        bamFile.h
        basicInfo.cpp
        basicInfo.h
        bedFile.cpp
        bedFile.h
        bpDecReg.cpp
        bpDecReg.h
        bpSideSolRel.cpp
        bpSideSolRel.h
        bpSolRel.cpp
        bpSolRel.h
        bpSolRelSet.cpp
        bpSolRelSet.h
        contig.h
        detectRegion.cpp
        detectRegion.h
        faFile.cpp
        faFile.h
        gcSta.cpp
        gcSta.h
        geneRegion.cpp
        geneRegion.h
        geneRegionFile.cpp
        geneRegionFile.h
        nanoSpeed.cpp
        nonB.cpp
        nonB.h
        nonbPotDetector.cpp
        nonbPotDetector.h
        oneVCF.cpp
        oneVCF.h
        relSet.cpp
        relSet.h
        solRel.cpp
        solRel.h
        srcfind.cpp
        ttbond_fa.cpp
        ttbond_math.h
        vcfFile.cpp
        vcfFile.h ttbond_fa.h nanoSpeedFile.cpp nanoSpeedFile.h ttbond_math.cpp nanoSpeedBF.cpp nanoSpeedBF.h fastqFile.cpp fastqFile.h fastqReads.cpp fastqReads.h oneFastq.cpp oneFastq.h alignManager.cpp alignManager.h comRel.cpp comRel.h alignInfo.cpp alignInfo.h)

target_link_libraries(srcfind hts)
