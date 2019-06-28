cmake_minimum_required(VERSION 2.8)
project(Sniffles)

find_package (ZLIB)
find_package (SeqAn REQUIRED)
find_package (htslib REQUIRED)


configure_file( Version.h.in ${CMAKE_SOURCE_DIR}/src/Version.h )

set(SEQAN_INCLUDE_DIRS "/data2/junwenwang/m204333/Project/seqan-v2.2.0/include")
set(HTSlib_INCLUDE_DIRS "/data2/junwenwang/m204333/Project/htslib/include")

include_directories (../lib/bamtools-2.3.0/src)
include_directories(../lib/tclap-1.2.1/include)
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories(${HTSlib_INCLUDE_DIRS})

add_definitions (${SEQAN_DEFINITIONS})

 IF(STATIC)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        SET(BUILD_SHARED_LIBRARIES OFF)
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
ENDIF()


# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Add executable and link against SeqAn dependencies.

add_executable(sniffles
                    tree/Breakpoint_Tree.cpp
                    Genotyper/Genotyper.cpp
					Alignment.cpp
					BamParser.cpp
					Sniffles.cpp
					Ignore_Regions.cpp
					tree/Intervall_bed.cpp
					sub/Detect_Breakpoints.cpp
					sub/Breakpoint.cpp
					tree/IntervallTree.cpp
					tree/IntervallList.cpp
					realign/SWCPU.cpp
					realign/Realign.cpp
					print/VCFPrinter.cpp
					print/BedpePrinter.cpp
					print/IPrinter.cpp
					tree/BinTree.cpp
					print/NGMPrinter.cpp
					plane-sweep/PlaneSweep_slim.cpp
					cluster/Cluster_SVs.cpp
					force_calling/Force_calling.cpp
					force_calling/VCF_parser.cpp
					)										
					
#target_link_libraries(ngm-core pthread)	
TARGET_LINK_LIBRARIES(sniffles BamTools-static)
TARGET_LINK_LIBRARIES(sniffles zlibstatic)
TARGET_LINK_LIBRARIES (sniffles ${SEQAN_LIBRARIES})
target_link_libraries(sniffles ${HTSlib_LIBRARIES})
	
add_executable(sniffles-debug 	
                    tree/Breakpoint_Tree.cpp
                    Genotyper/Genotyper.cpp
					Alignment.cpp
					BamParser.cpp
					Sniffles.cpp
					Ignore_Regions.cpp
					tree/Intervall_bed.cpp
					tree/IntervallList.cpp
					sub/Detect_Breakpoints.cpp
					sub/Breakpoint.cpp
					tree/IntervallTree.cpp
					realign/SWCPU.cpp
					realign/Realign.cpp
					print/VCFPrinter.cpp
					print/BedpePrinter.cpp
					print/IPrinter.cpp
					tree/BinTree.cpp
					print/NGMPrinter.cpp
					plane-sweep/PlaneSweep_slim.cpp
					cluster/Cluster_SVs.cpp
					force_calling/Force_calling.cpp
					force_calling/VCF_parser.cpp
					)	

SET_TARGET_PROPERTIES(sniffles-debug PROPERTIES COMPILE_FLAGS "-g3 -O0")



#target_link_libraries(sniffles-debug pthread)
TARGET_LINK_LIBRARIES(sniffles-debug BamTools-static)
TARGET_LINK_LIBRARIES(sniffles-debug zlibstatic)
TARGET_LINK_LIBRARIES (sniffles-debug ${SEQAN_LIBRARIES})
target_link_libraries(sniffles-debug ${HTSlib_LIBRARIES})



