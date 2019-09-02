/*
 * Detect_Breakpoints.h
 *
 *  Created on: Jun 19, 2015
 *      Author: fsedlaze
 */

#ifndef SUB_DETECT_BREAKPOINTS_H_
#define SUB_DETECT_BREAKPOINTS_H_

#include "../BamParser.h"
#include "../Parser.h"
#include "../Alignment.h"
#include "../plane-sweep/Plane-sweep.h"
#include "../tree/IntervallTree.h"
#include "../tree/TNode.h"
#include "../tree/IntervallContainer.h"
#include "../tree/IntervallList.h"
#include "../Paramer.h"
#include "../print/IPrinter.h"

#include <iostream>
#include <omp.h>

struct hist_str{
	long position;
	int hits;
	std::vector<std::string> names;
};

struct tra_str{
    long position;
    int hits;
    int sameStrand_hits = 0;
    int diffStrand_hits = 0;
    std::map<long, vector<long>> map_pos;
    std::vector<std::string> names;

};


struct BpRln{
//    int originalSupport;
    bool isSameStrand;
    pair<long, long> coordinate;
    pair<long, long> chr_pos;
    pair<std::string, std::string> chr;
    pair<int, int> chr_idx;
    Breakpoint * bp;

    BpRln(bool strand, pair<long, long> coor, const RefVector ref, Breakpoint * breakpoint){
        isSameStrand = strand;
        coordinate = coor;
        chr_pos.first = IPrinter::calc_pos(coordinate.first, ref, chr_idx.first);
        chr_pos.second = IPrinter::calc_pos(coordinate.second, ref, chr_idx.second);
        chr.first = ref[chr_idx.first].RefName;
        chr.second =ref[chr_idx.second].RefName;

        bp = breakpoint;
    }

    bool operator < (const BpRln& bp) const {
        if (chr.first == bp.chr.first)
            return chr_pos.first < bp.chr_pos.first;
//        else return std::stoi(chr.first.substr(3, chr.first.size())) < std::stoi(bp.chr.first.substr(3, bp.chr.first.size()));
        else return chr_idx.first < bp.chr_idx.first;
    }
};


void clarify(std::vector<Breakpoint *> & points);
void detect_breakpoints(std::string filename, IPrinter *& printer);
//void screen_for_events(Node * list,IntervallTree & bst ,TNode *&root, int cov, int lowMQ_cov,RefVector ref);
bool screen_for_events(Alignment * tmp, IntervallTree & bst, TNode *&root, RefVector ref, int cov);
void add_events(Alignment *& tmp, std::vector<str_event> events, short type, long ref_space, IntervallTree & bst, TNode *&root,long read_id,bool add);
void add_splits(Alignment *& tmp, std::vector<aln_str> events, short type, RefVector ref, IntervallTree & bst, TNode *&root,long read_id,bool add);
void estimate_parameters(std::string read_filename);
bool overlaps(aln_str prev,aln_str curr);
void detect_merged_svs(Breakpoint * point);

std::string TRANS_type(char type);



#endif /* SUB_DETECT_BREAKPOINTS_H_ */
