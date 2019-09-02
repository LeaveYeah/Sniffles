/*
 * Detect_Breapoints.cpp
 *
 *  Created on: Jun 19, 2015
 *      Author: fsedlaze
 */

#include "Detect_Breakpoints.h"
#include "../print/IPrinter.h"
#include "seqan/align.h"
#include "seqan/bam_io.h"
#include "../Alignment.h"
#include "../bioio.hpp"
#include <chrono>

struct map_example {
    bool strand;
    string chr;
    pair<int, int> diff;
    pair<long, long>ref_pos;
    pair<long, long> read_pos;
    Alignment * tmp_aln;
    string bases;
};

void store_pos(vector<hist_str> &positions, long pos, std::string read_name) {
	for (size_t i = 0; i < positions.size(); i++) {
		if (abs(positions[i].position - pos) < Parameter::Instance()->min_length) {
			positions[i].hits++;
			positions[i].names.push_back(read_name);
			return;
		}
	}
	hist_str tmp;
	tmp.position = pos;
	tmp.hits = 1;
	tmp.names.push_back(read_name);
	positions.push_back(tmp);
}

std::string reverse_complement(std::string sequence) {
	std::string tmp_seq;
	for (std::string::reverse_iterator i = sequence.rbegin(); i != sequence.rend(); i++) {
		switch ((*i)) {
		case 'A':
			tmp_seq += 'T';
			break;
		case 'C':
			tmp_seq += 'G';
			break;
		case 'G':
			tmp_seq += 'C';
			break;
		case 'T':
			tmp_seq += 'A';
			break;
		default:
			tmp_seq += (*i);
			break;
		}
	}
	return tmp_seq;
}

Breakpoint * split_points(vector<std::string> names, std::map<std::string, read_str> support) {
	std::map<std::string, read_str> new_support;
	for (size_t i = 0; i < names.size(); i++) {
		new_support[names[i]] = support[names[i]];
	}
	position_str svs;
	svs.start.min_pos = 0; //just to initialize. Should not be needed anymore at this stage of the prog.
	svs.stop.max_pos = 0;
	svs.support = new_support;
	Breakpoint * point = new Breakpoint(svs, (*new_support.begin()).second.coordinates.second - (*new_support.begin()).second.coordinates.first);
	return point;
}


void detect_merged_svs(position_str point, RefVector ref, vector<Breakpoint *> & new_points) {
	new_points.clear(); //just in case!
	vector<hist_str> pos_start;
	vector<hist_str> pos_stop;
	for (std::map<std::string, read_str>::iterator i = point.support.begin(); i != point.support.end(); ++i) {

		store_pos(pos_start, (*i).second.coordinates.first, (*i).first);
		store_pos(pos_stop, (*i).second.coordinates.second, (*i).first);
//		string chr_start, chr_stop;
//		long pos_start = IPrinter::calc_pos((*i).second.coordinates.first, ref, chr_start);
//        long pos_stop = IPrinter::calc_pos((*i).second.coordinates.second, ref, chr_stop);
//        std::cout << (*i).first << endl;
//        std::cout << chr_start << " " << pos_start << endl;
//        std::cout << chr_stop << " " << pos_stop << endl;

	}
    vector<size_t> start_idx, stop_idx;
	int start_count = 0;
	for (size_t i = 0; i < pos_start.size(); i++) {
		//std::cout<<pos_start[i].hits <<",";
		if (pos_start[i].hits > Parameter::Instance()->min_support) {
			start_count++;
            start_idx.push_back(i);
		}

	}
	int stop_count = 0;
	for (size_t i = 0; i < pos_stop.size(); i++) {
		//	std::cout << pos_stop[i].hits << ",";
		if (pos_stop[i].hits > Parameter::Instance()->min_support) {
			stop_count++;
            stop_idx.push_back(i);
		}
	}
	if (stop_count > 1 || start_count > 1) {
		std::cout << "\tprocessing merged TRA" << std::endl;
		if (start_count > 1) {
			new_points.push_back(split_points(pos_start[start_idx[0]].names, point.support));
			new_points.push_back(split_points(pos_start[start_idx[1]].names, point.support));
		} else {
			new_points.push_back(split_points(pos_stop[stop_idx[0]].names, point.support));
			new_points.push_back(split_points(pos_stop[stop_idx[1]].names, point.support));
		}
	}
}

void fix_mismatch_read_pos(vector<differences_str> &event_aln, int i, Alignment * tmp_aln){
    int read_pos = 0;
    if (i == 0) {
        if (tmp_aln->getAlignment()->CigarData[0].Type == 'S')
            read_pos += tmp_aln->getAlignment()->CigarData[0].Length;
        read_pos += event_aln[i].position - tmp_aln->getPosition();
        event_aln[i].readposition = read_pos;
    } else
        event_aln[i].readposition = event_aln[i-1].readposition + event_aln[i].position - event_aln[i-1].position -
                event_aln[i-1].type;
}
void store_tra_pos(vector<tra_str> &positions, read_str read, std::string read_name, bool start) {
    long pos, opp_pos;
    if (start) {
        pos = read.coordinates.first;
        opp_pos = read.coordinates.second;
    } else {
        opp_pos = read.coordinates.second;
        pos = read.coordinates.first;
    }


    for (size_t i = 0; i < positions.size(); i++) {
        if (abs(positions[i].position - pos) < Parameter::Instance()->min_length) {
            positions[i].hits++;
            positions[i].names.push_back(read_name);
            positions[i].map_pos[pos].push_back(opp_pos);
            if (read.read_strand.first == read.read_strand.second)
                positions[i].sameStrand_hits++;
            else positions[i].diffStrand_hits++;
            return;
        }
    }
    tra_str tmp;
    tmp.position = pos;
    tmp.hits = 1;
    tmp.names.push_back(read_name);
    tmp.map_pos[pos].push_back(opp_pos);
    if (read.read_strand.first == read.read_strand.second)
        tmp.sameStrand_hits = 1;
    else tmp.diffStrand_hits = 1;
    positions.push_back(tmp);
}

long get_max_pos(map<long, vector<long>> map_pos, int &max){
    long coordinate = 0;
    for (auto i = map_pos.begin(); i != map_pos.end(); i++){
        if ((*i).second.size() > max){
            max = (*i).second.size();
            coordinate = (*i).first;
        }
    }
    return coordinate;
}

long get_max_opp_pos(vector<long> opp_pos, int &max){
    map<long, int> map_opp_pos;
    for (auto j = opp_pos.begin(); j != opp_pos.end(); j++)
        if (map_opp_pos.find(*j) == map_opp_pos.end())
            map_opp_pos[*j] = 1;
        else map_opp_pos[*j]++;
    long coordinate = 0;

    for (auto i = map_opp_pos.begin(); i != map_opp_pos.end(); i++){
        if ((*i).second > max){
            max = (*i).second;
            coordinate = (*i).first;
        }
    }
    return coordinate;
}

bool cal_high_error_side(vector<differences_str> &event_aln, long pos, long distance, Alignment * tmp_aln){
    if (event_aln.empty())
        return false;
    for (size_t j = 0; j < event_aln.size(); j++) {
        if (event_aln[j].type == 0)
            fix_mismatch_read_pos(event_aln, j, tmp_aln);
    }
    int left_hits = 0;
    int right_hits = 0;
    differences_str evt_left = event_aln[0];
    differences_str evt_right = event_aln[event_aln.size()-1];
    for (auto i : event_aln){
        if (i.position < pos - distance) continue;
        if (i.position > pos + distance) break;
        if (i.position >= pos - distance && i.position <= pos){
            left_hits += max((int) abs(i.type), 1);
            evt_left = i;
        }
        else if (i.position >= pos && i.position <= pos + distance) {
            right_hits += max((int) abs(i.type), 1);
            if (evt_right.position > i.position) evt_right = i;
        }

    }

    if (abs(evt_left.position - pos) < abs(evt_right.position - pos)) {
        tmp_aln->bp_read_pos = evt_left.readposition;
        tmp_aln->diff = evt_left.position - pos;
    } else if (abs(evt_right.position - pos) < abs(evt_left.position - pos)) {
        tmp_aln->bp_read_pos = evt_right.readposition;
        tmp_aln->diff = evt_right.position - pos;
    }

    if (right_hits - left_hits >= distance / 5) {
        tmp_aln->high_error_region_score = -1 * right_hits;
        if (abs(evt_left.position - pos) == abs(evt_right.position - pos)) {
            tmp_aln->bp_read_pos = evt_left.readposition;
            tmp_aln->diff = evt_left.position - pos;
        }
        tmp_aln->high_error_side = true;
        return true;
    }

    if (left_hits - right_hits >= distance / 5) {
        tmp_aln->high_error_region_score = -1 * left_hits;
        if (abs(evt_left.position - pos) == abs(evt_right.position - pos)) {
            tmp_aln->bp_read_pos = evt_right.readposition;
            tmp_aln->diff = evt_right.position - pos;
        }
        tmp_aln->high_error_side = false;
        return true;
    }

    return false;
}

int map_read(Alignment  * tmp_aln, BpRln bp, int distance,
        const bioio::FastaIndex  index, std::ifstream & fasta){
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    using namespace seqan;
    typedef String<char> TSequence;
    typedef StringSet<TSequence> TStringSet;
    typedef Align<TSequence, ArrayGaps> TAlign;// container for strings
    typedef StringSet<TSequence, Dependent<>> TDepStringSet;
    typedef seqan::Alignment<TDepStringSet> TAlignStringSet; // dependent string set
    typedef Graph<TAlignStringSet> TAlignGraph;       // alignment graph// sequence type


    if (bp.isSameStrand) bp.chr_pos.second += tmp_aln->diff;
    else bp.chr_pos.second  -= tmp_aln->diff;

    if (!tmp_aln->high_error_side) tmp_aln->bp_read_pos -= distance; //lefthand side
    if (bp.isSameStrand != tmp_aln->high_error_side) bp.chr_pos.second -= distance;

    string ref_str = bioio::read_fasta_contig(fasta, index.at(bp.chr.second), bp.chr_pos.second, distance);
    transform(ref_str.begin(), ref_str.end(), ref_str.begin(), ::toupper);
    TSequence reference = ref_str;
    if (tmp_aln->bp_read_pos + distance > tmp_aln->getQueryBases().size())
        return -50;
    if (tmp_aln->bp_read_pos < 0)
        return -50;
//    std::cout << bp.chr.second << " " << bp.chr_pos.second << endl;
//    std::cout << tmp_aln->bp_read_pos << " " << distance << endl;
    string seq_str = tmp_aln->getQueryBases().substr(tmp_aln->bp_read_pos, distance);
    transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
    TSequence sequence = seq_str;
    TStringSet sequences;
    appendValue(sequences, reference);
    appendValue(sequences, sequence);

    TAlignGraph alignG(sequences);

    int score = globalAlignment(alignG, Score<int, Simple>(0, -1, -1), AlignConfig<false, false, true, true>(), LinearGaps());

    return score;

}

vector<CigarOp> map_read_example(string ref_str, string seq_str){
    using namespace seqan;
    typedef String<char> TSequence;
    typedef StringSet <TSequence> TStringSet;
    typedef Align <TSequence, ArrayGaps> TAlign;// container for strings
    typedef StringSet <TSequence, Dependent<>> TDepStringSet;
    typedef seqan::Alignment <TDepStringSet> TAlignStringSet;
    typedef Row<TAlign>::Type TRow;
    typedef Iterator<TRow>::Type TRowIterator;

    TAlign align;


    transform(ref_str.begin(), ref_str.end(), ref_str.begin(), ::toupper);
    TSequence reference = ref_str;

    transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
    TSequence sequence = seq_str;


    resize(rows(align), 2);
    assignSource(row(align, 0), reference);
    assignSource(row(align, 1), sequence);

    TRow & row1 = row(align, 0);
    TRow & row2 = row(align, 1);

    int score = globalAlignment(align, Score<int, Simple>(0, -2, -1, -2));
    cout << score << endl;
    cout << align << endl;
    TRowIterator itRef = begin(row1);
    TRowIterator itSeq = begin(row2);
    vector<CigarOp> cigar;
    if (seq_str.size() == 0)
        return cigar;

    int numMatchAndMismatches = 0;
    while (itSeq != end(row2))
    {
        // Count insertions.
        if (isGap(itRef))
        {
            int numGaps = countGaps(itRef);
            CigarOp cigarOp('I', numGaps );
            cigar.push_back(cigarOp);
            itRef += numGaps;
            itSeq += numGaps;
            continue;
        }
        // Count deletions.
        if (isGap(itSeq))
        {
            int numGaps = countGaps(itSeq);
            CigarOp cigarOp('D', numGaps );
            cigar.push_back(cigarOp);
            itRef += numGaps;
            itSeq += numGaps;
            continue;
        }

        // Count matches and  mismatches.
        while (itSeq != end(row2))
        {
            if (isGap(itSeq) || isGap(itRef))
                break;

            ++numMatchAndMismatches;
            ++itRef;
            ++itSeq;
        }
        if (numMatchAndMismatches != 0) {
            CigarOp cigarOp('M', numMatchAndMismatches );
            cigar.push_back(cigarOp);
        }
        numMatchAndMismatches = 0;
    }

    cout << endl;
    return cigar;

}

//vector<differences_str> get_aln_events(Alignment * tmp_aln){
//    vector<differences_str> event_aln;
//    vector<indel_str> dels;
//    event_aln = tmp_aln->summarizeAlignment(dels);
//    for (size_t i = 0; i < event_aln.size(); i++){
//        if (event_aln[i].readposition == -1)
//            if (i > 0)
//                event_aln[i].readposition ;
//        ev.resolved = true;
//        ev.readposition = read_pos;
//        ev.type = al->CigarData[i].Length * -1; //insertion
//    }
//}

std::string TRANS_type(char type) {
	string tmp;
	if (type & DEL) {
		tmp += "DEL";
	}
	if (type & INV) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "INV";
	}
	if (type & DUP) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "DUP";
	}
	if (type & INS) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "INS";
	}
	if (type & TRA) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "TRA";
	}
	if (type & NEST) {
		if (!tmp.empty()) {
			tmp += '/';
		}
		tmp += "NEST";
	}
	return tmp; // should not occur!
}

long get_ref_lengths(int id, RefVector ref) {
	long length = 0;

	for (size_t i = 0; i < (size_t) id && i < ref.size(); i++) {
		length += (long) ref[i].RefLength + (long) Parameter::Instance()->max_dist;
	}
	return length;
}

void detect_bp_for_realn(Breakpoint  *breakpoint, const RefVector ref, vector<BpRln> &bp_rlns) {
    auto point = breakpoint->get_coordinates();
    vector<tra_str> pos_start;

    vector<tra_str> pos_stop;


    for (auto i = point.support.begin(); i != point.support.end(); ++i) {
        store_tra_pos(pos_start, (*i).second, (*i).first, true);
        store_tra_pos(pos_stop, (*i).second, (*i).first, false);
    }

    for (size_t i = 0; i < pos_start.size(); i++) {

        if (pos_start[i].hits >= Parameter::Instance()->min_support / 3) {
//        if (pos_start[i].hits >= 1) {

            bool isSameStrand = pos_start[i].sameStrand_hits >= pos_start[i].diffStrand_hits;
            pair<long, long> coordinate;
            int max_start = 0, max_stop = 0;
            coordinate.first = get_max_pos(pos_start[i].map_pos, max_start);
            coordinate.second = get_max_opp_pos(pos_start[i].map_pos[coordinate.first], max_stop);

            BpRln start_bp(isSameStrand, coordinate, ref, breakpoint);
            std::swap(coordinate.first, coordinate.second);
            BpRln stop_bp(isSameStrand, coordinate, ref, breakpoint);

            bp_rlns.push_back(start_bp);
            bp_rlns.push_back(stop_bp);
        }
    }
}

bool should_be_stored(Breakpoint *& point) {
	point->calc_support(); // we need that before:
	//std::cout << "Stored: " << point->get_support() << " " << point->get_length() << std::endl;
	if (point->get_SVtype() & TRA) { // we cannot make assumptions abut support yet.
		point->set_valid((bool) (point->get_support() > 1)); // this is needed as we take each chr independently and just look at the primary alignment
	} else if (point->get_support() >= Parameter::Instance()->min_support) {
		point->predict_SV();
		point->set_valid((bool) (point->get_length() > Parameter::Instance()->min_length));
	}
	return point->get_valid();
}
void polish_points(std::vector<Breakpoint *> & points, RefVector ref) { //TODO might be usefull! but why does the tree not fully work??
	return;
	for (size_t i = 0; i < points.size(); i++) {
		if (points[i]->get_SVtype() & INS && (points[i]->get_length() == Parameter::Instance()->huge_ins)) {
			for (size_t j = 0; j < points.size(); j++) {
				if (i != j) {
					if (abs(points[i]->get_coordinates().start.min_pos - points[j]->get_coordinates().start.min_pos) < Parameter::Instance()->max_dist || abs(points[i]->get_coordinates().stop.max_pos - points[j]->get_coordinates().stop.max_pos) < Parameter::Instance()->max_dist) {
						std::cout << "HIT!: " << points[j]->get_coordinates().start.min_pos << " " << points[i]->get_coordinates().start.min_pos << " " << points[j]->get_coordinates().stop.max_pos << " " << points[i]->get_coordinates().stop.max_pos << " len: " << points[j]->get_length() << " " << points[i]->get_length() << std::endl;
						break;
					}
				}

			}
		}
	}
}

void add_realign_read(BpRln bp_realign,  Alignment * tmp_aln){
    read_str tmp;
    if (bp_realign.coordinate.first < bp_realign.coordinate.second) {
        tmp.coordinates.first = bp_realign.coordinate.first + tmp_aln->diff;
        if (bp_realign.isSameStrand)
            tmp.coordinates.second = bp_realign.coordinate.second + tmp_aln->diff;
        else tmp.coordinates.second = bp_realign.coordinate.second - tmp_aln->diff;
    } else {
        tmp.coordinates.second = bp_realign.coordinate.first + tmp_aln->diff;
        if (bp_realign.isSameStrand)
            tmp.coordinates.first = bp_realign.coordinate.second + tmp_aln->diff;
        else tmp.coordinates.first = bp_realign.coordinate.second - tmp_aln->diff;
    }
    tmp.SV = TRA;
    tmp.type = 1;

//                        std::cout << bp_realign.bp->get_coordinates().support.size() << endl;
    auto map = bp_realign.bp->get_coordinates().support;
    if (map.find(tmp_aln->getName()+"_rln_0") == map.end())
        bp_realign.bp->add_read(tmp, tmp_aln->getName()+"_rln_0");
    else
        bp_realign.bp->add_read(tmp, tmp_aln->getName() + "_rln_1");
}


bool detect_gap(Alignment* tmp_aln, long ref_pos, int distance, RefVector ref) {
    std::vector <aln_str> split_events = tmp_aln->getSA(ref);
    if (split_events.size() == 0) return false;
    long read_length;
    aln_str last_event = split_events[split_events.size()-1];
    if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800) )
        read_length = tmp_aln->getAlignment()->Length;
    else if (last_event.strand) {
        read_length = last_event.read_pos_stop;
        if (last_event.cigar.back().Type == 'S')
            read_length += last_event.cigar.back().Type;
    } else {
        read_length = last_event.read_pos_stop;
        if (last_event.cigar.front().Type == 'S')
            read_length += last_event.cigar.front().Type;
    }
    if (split_events[0].isMain && split_events[0].read_pos_start >= distance) {
        if (split_events[0].strand && abs(split_events[0].pos - ref_pos) < distance) {
            tmp_aln->diff = split_events[0].pos - ref_pos;
            tmp_aln->high_error_side = false;
            tmp_aln->bp_read_pos = split_events[0].read_pos_start;
            return true;
        } else if (!split_events[0].strand && abs(split_events[0].pos + split_events[0].length - ref_pos) < distance) {
            tmp_aln->diff = split_events[0].pos + split_events[0].length - ref_pos;
            tmp_aln->high_error_side = true;
            tmp_aln->bp_read_pos = read_length - 1 - split_events[0].read_pos_start;
            return true;
        }
    }

    if (last_event.isMain && tmp_aln->getAlignment()->Length - 1 - last_event.read_pos_stop >= distance) {
        if (last_event.strand && abs(last_event.pos + last_event.length - ref_pos) < distance) {
            tmp_aln->diff = last_event.pos + last_event.length - ref_pos;
            tmp_aln->high_error_side = true;
            tmp_aln->bp_read_pos = last_event.read_pos_stop;
            return true;
        } else if (!last_event.strand && abs(last_event.pos - ref_pos) < distance) {
            tmp_aln->diff = last_event.pos - ref_pos;
            tmp_aln->high_error_side = false;
            tmp_aln->bp_read_pos = read_length - 1 - last_event.read_pos_stop;
            return true;
        }
    }

    for (size_t i = 0; i < split_events.size(); i++) {
        if (split_events[i].isMain) {
            if (abs(split_events[i].pos - ref_pos) < distance) {
                tmp_aln->diff = split_events[i].pos - ref_pos;
                tmp_aln->high_error_side = false;
                if (split_events[i].strand && i > 0
                && abs(split_events[i].read_pos_start - split_events[i-1].read_pos_stop) >= distance) {
                    tmp_aln->bp_read_pos = split_events[i].read_pos_start;
                    return true;
                } else if (!split_events[i].strand && i < split_events.size() - 1
                && abs(split_events[i].read_pos_stop - split_events[i+1].read_pos_start) >= distance) {
                    tmp_aln->bp_read_pos = read_length - 1 - split_events[i].read_pos_stop;
                    return true;
                } else return false;
            } else if (abs(split_events[i].pos + split_events[i].length - ref_pos) < distance) {
                tmp_aln->diff = split_events[i].pos + split_events[i].length - ref_pos;
                tmp_aln->high_error_side = true;
                if (split_events[i].strand && i < split_events.size() - 1
                && abs(split_events[i].read_pos_stop - split_events[i+1].read_pos_start) >= distance) {
                    tmp_aln->bp_read_pos = split_events[i].read_pos_stop;
                    return true;
                } else if (!split_events[i].strand && i > 0
                && abs(split_events[i].read_pos_start - split_events[i-1].read_pos_stop) >= distance) {
                    tmp_aln->bp_read_pos =  read_length - 1 - split_events[i].read_pos_start;
                    return true;
                }
                else return false;
            }
        }
    }
    return false;
}


void realign_read(BpRln bp_realign, vector<differences_str> &event_aln, Alignment* tmp_aln,
                  const bioio::FastaIndex  index, std::ifstream & fasta, RefVector ref) {
    int distance = min(100, Parameter::Instance()->min_length);
    if (bp_realign.chr_pos.first - distance <= tmp_aln->getPosition() ||
        bp_realign.chr_pos.first + distance >= tmp_aln->getPosition() + tmp_aln->getRefLength()) {
        if (!Parameter::Instance()->realn_clipped)
            return;
        auto map = bp_realign.bp->get_coordinates().support;
        if (map.find(tmp_aln->getName()) != map.end()) return;
        bool existsGap = detect_gap(tmp_aln, bp_realign.chr_pos.first, distance,  ref);
        if (!existsGap) return;
        int alt_aln_score = map_read(tmp_aln, bp_realign,  distance, index, fasta);
        if (alt_aln_score > -0.2 * distance)
            add_realign_read(bp_realign,  tmp_aln);
    } else {

        bool exists_high_error_side = cal_high_error_side(event_aln, bp_realign.chr_pos.first, distance,  tmp_aln);
        if (!exists_high_error_side) return;
        int alt_aln_score = map_read(tmp_aln, bp_realign,  distance, index, fasta);
        if (alt_aln_score - tmp_aln->high_error_region_score > distance / 5 &&
            alt_aln_score > -0.2 * distance) {
            add_realign_read(bp_realign,  tmp_aln);
        }
    }
}

void detect_breakpoints(std::string read_filename, IPrinter *& printer) {

	estimate_parameters(read_filename);
	BamParser * mapped_file = 0;
	RefVector ref;
	if (read_filename.find("bam") != string::npos) {
		mapped_file = new BamParser(read_filename);
		ref = mapped_file->get_refInfo();
	} else {
		cerr << "File Format not recognized. File must be a sorted .bam file!" << endl;
		exit(EXIT_FAILURE);
	}
//Using PlaneSweep to comp coverage and iterate through reads:
//PlaneSweep * sweep = new PlaneSweep();
	//Using Interval tree to store and manage breakpoints:

	IntervallTree final;
	IntervallTree bst;

	TNode * root_final = NULL;
	int current_RefID = 0;


	TNode *root = NULL;
//FILE * alt_allel_reads;
	FILE * ref_allel_reads;
	if (Parameter::Instance()->genotype) {
		ref_allel_reads = fopen(Parameter::Instance()->tmp_genotyp.c_str(), "w");
		//	ref_allel_reads = fopen(Parameter::Instance()->tmp_genotyp.c_str(), "wb");
	}

	Alignment * tmp_aln = mapped_file->parseRead(Parameter::Instance()->min_mq);

//	std::deque<Alignment*> overlapAlignments;

	long ref_space = get_ref_lengths(tmp_aln->getRefID(), ref);
	long num_reads = 0;

	/*Genotyper * go;
	 if (Parameter::Instance()->genotype) {
	 go = new Genotyper();
	 }*/
	std::cout << "Start parsing... " << ref[tmp_aln->getRefID()].RefName << std::endl;

//filter and copy results:
    while (!tmp_aln->getQueryBases().empty()) {

        if ((tmp_aln->getAlignment()->IsPrimaryAlignment()) && (!(tmp_aln->getAlignment()->AlignmentFlag & 0x800) && tmp_aln->get_is_save())) {	// && (Parameter::Instance()->chr_names.empty() || Parameter::Instance()->chr_names.find(ref[tmp_aln->getRefID()].RefName) != Parameter::Instance()->chr_names.end())) {
            //change CHR:
            if (current_RefID != tmp_aln->getRefID()) {

                std::cout << "\tSwitch Chr " << ref[tmp_aln->getRefID()].RefName << std::endl;	//" " << ref[tmp_aln->getRefID()].RefLength
                std::vector<Breakpoint *> points;
                bst.get_breakpoints(root, points);
                //polish_points(points, ref);

                /*	if (Parameter::Instance()->genotype) {
                 fclose(ref_allel_reads);
                 cout<<"\t\tGenotyping"<<endl;
                 go->update_SVs(points, ref_space);
                 cout<<"\t\tGenotyping finished"<<endl;
                 ref_allel_reads = fopen(Parameter::Instance()->tmp_genotyp.c_str(), "wb");
                 }*/

                for (int i = 0; i < points.size(); i++) {
                    points[i]->calc_support();
                    if (points[i]->get_valid()) {
                        //invoke update over ref support!
                        if (points[i]->get_SVtype() & TRA) {
                            final.insert(points[i], root_final);
                        } else {
                            printer->printSV(points[i]);
                        }
                    }
                }
                bst.clear(root);
                current_RefID = tmp_aln->getRefID();
                ref_space = get_ref_lengths(tmp_aln->getRefID(), ref);
//				overlapAlignments.clear();
            }

            //SCAN read:
            std::vector<str_event> aln_event;
            std::vector<aln_str> split_events;
            if (tmp_aln->getMappingQual() > Parameter::Instance()->min_mq) {
                double score = tmp_aln->get_scrore_ratio();

                //

#pragma omp parallel // starts a new team
                {
#pragma omp sections
                    {
#pragma omp section
                        {
                            //		clock_t begin = clock();
                            if ((score == -1 || score > Parameter::Instance()->score_treshold)) {
                                aln_event = tmp_aln->get_events_Aln();
                            }
                            //		Parameter::Instance()->meassure_time(begin, " Alignment ");
                        }
#pragma omp section
                        {
                            //	clock_t begin_split = clock();
                            split_events = tmp_aln->getSA(ref);
                            //	Parameter::Instance()->meassure_time(begin_split, " Split reads ");
                        }
                    }
                }
                //tmp_aln->set_supports_SV(aln_event.empty() && split_events.empty());

                //Store reference supporting reads for genotype estimation:

                bool SV_support = (!aln_event.empty() && !split_events.empty());
                if (Parameter::Instance()->genotype && !SV_support) {
                    //	cout << "STORE" << endl;
                    //write read:
                    /*str_read tmp;
                    tmp.chr_id = tmp_aln->getRefID();	//check string in binary???
                    tmp.start = tmp_aln->getPosition();
                    tmp.length = tmp_aln->getRefLength();
                    if (tmp_aln->getStrand()) {
                        tmp.strand = 1;
                    } else {
                        tmp.strand = 2;
                    }*/
                    write_read(tmp_aln, ref_allel_reads);
                    //fwrite(&tmp, sizeof(struct str_read), 1, ref_allel_reads);
                }

                //store the potential SVs:
                if (!aln_event.empty()) {
                    add_events(tmp_aln, aln_event, 0, ref_space, bst, root, num_reads, false);
                }
                if (!split_events.empty()) {
                    add_splits(tmp_aln, split_events, 1, ref, bst, root, num_reads, false);
                    //realign the reads in overlapAlignments

                }
            }
        }

        mapped_file->parseReadFast(Parameter::Instance()->min_mq, tmp_aln);

        num_reads++;

        if (num_reads % 10000 == 0) {
            cout << "\t\t# Processed reads: " << num_reads << endl;
        }
    }
    std::cout << "Finalizing  .." << std::endl;
	std::vector<Breakpoint *> points;
	bst.get_breakpoints(root, points);

	/*	if (Parameter::Instance()->genotype) {
	 fclose(ref_allel_reads);
	 go->update_SVs(points, ref_space);
	 string del = "rm ";
	 del += Parameter::Instance()->tmp_genotyp;
	 del += "ref_allele";
	 system(del.c_str());
	 }*/

	for (int i = 0; i < points.size(); i++) {
		points[i]->calc_support();
		if (points[i]->get_valid()) {
			//invoke update over ref support!
			if (points[i]->get_SVtype() & TRA) {
				final.insert(points[i], root_final);
			} else {
				printer->printSV(points[i]);
			}
		}
	}
	bst.clear(root);
	points.clear();
	final.get_breakpoints(root_final, points);
	//std::cout<<"Detect merged tra"<<std::endl;
	size_t points_size = points.size();


    vector<BpRln> bp_rlns;
    for (size_t i = 0; i < points_size; i++) { // its not nice, but I may alter the length of the vector within the loop.
        if (points[i]->get_SVtype() & TRA) {
            detect_bp_for_realn(points[i], ref, bp_rlns);
        }
    }

    mapped_file->Rewind();
    std::sort(bp_rlns.begin(), bp_rlns.end());
    vector<BpRln> active_bp;

    auto i = bp_rlns.begin();
    int distance = min(100, Parameter::Instance()->min_length);
    mapped_file->Jump(i->chr_idx.first, max((long)0, i->chr_pos.first - distance));
    tmp_aln = mapped_file->parseRead(Parameter::Instance()->min_mq);

    const auto index = bioio::read_fasta_index(Parameter::Instance()->fasta_index_file);
    std::ifstream fasta {Parameter::Instance()->fasta_file, std::ios::binary};
    ofstream bam_out;
    bam_out.open ("/data2/junwenwang/m204333/Project/sniffles/out/analysis/example.sam");


    num_reads = 0;
    map<string, map_example> example_map;
    while (!tmp_aln->getQueryBases().empty() && (i != bp_rlns.end() || !active_bp.empty())) {
        if (tmp_aln->get_is_save()) {
//            std::cout << ref[tmp_aln->getRefID()].RefName << " " << tmp_aln->getPosition() << endl;
//            cout << "step0" << endl;
            vector<differences_str> event_aln;
            vector<indel_str> dels;
            event_aln = tmp_aln->summarizeAlignment(dels);
            while (i != bp_rlns.end()) {
                if (i->chr_idx.first != tmp_aln->getRefID())
                    break;
                if (i->chr_pos.first + distance >= tmp_aln->getPosition() &&
                    i->chr_pos.first - distance <= tmp_aln->getPosition() + tmp_aln->getRefLength())
                    active_bp.push_back(*i);
                else break;
                i++;
            }
//            cout << "step1" << endl;
            size_t num_rm = 0;
            for (auto j = active_bp.begin(); j != active_bp.end(); j++) {

                if (j->chr_idx.first != tmp_aln->getRefID() || j->chr_pos.first + distance < tmp_aln->getPosition()) {
                    num_rm++;
                } else break;
            }
//            cout << "step2" << endl;
            if (num_rm != 0) {
//                cout << "before step2.1" << endl;
                active_bp.erase(active_bp.begin(), active_bp.begin() + num_rm);
//                cout << "after step2.1" << endl;
            }

            if (!active_bp.empty()) {
                for (auto j = active_bp.begin(); j != active_bp.end(); j++) {
                    if (tmp_aln->getName() == "52382413-3302-4d5f-9a2e-1c4a244442d7") {
                        tmp_aln->bp_read_pos += 0;
                        tmp_aln->bp_read_pos += 0;
                    }

                    realign_read(*j, event_aln, tmp_aln, index, fasta, ref);
                    if (tmp_aln->bp_read_pos == 0)
                        continue;
//                    if (tmp_aln->bp_read_pos < 0 || tmp_aln->bp_read_pos < tmp_aln->getAlignment()->Length) {
//                        tmp_aln->bp_read_pos == -1;
//                         tmp_aln->bp_read_pos += 0;
//
//                     }
                    if (j->chr_pos.first == 4007175 || j->chr_pos.first == 163087117 || j->chr_pos.first == 4008128 || j->chr_pos.first == 163088072)  {
                        if (example_map.find(tmp_aln->getName()) == example_map.end()) {

                            map_example example;
                            example.strand = tmp_aln->getStrand();

                            if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800) )
                                example.bases = tmp_aln->getQueryBases();
                            example.tmp_aln = tmp_aln;
                            if (tmp_aln->high_error_side)
                                example.read_pos.first = tmp_aln->bp_read_pos;
                            else example.read_pos.first = tmp_aln->bp_read_pos + distance;
                            example.ref_pos.first = j->chr_pos.second;


                            example.diff.first = tmp_aln->diff;
                            example_map[tmp_aln->getName()] = example;
                        }  else {
                            example_map[tmp_aln->getName()].chr = "";
                            if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800))
                                example_map[tmp_aln->getName()].bases = tmp_aln->getQueryBases();

                            if (tmp_aln->high_error_side)
                                example_map[tmp_aln->getName()].read_pos.second = tmp_aln->bp_read_pos;
                            else example_map[tmp_aln->getName()].read_pos.second = tmp_aln->bp_read_pos + distance;
                            example_map[tmp_aln->getName()].ref_pos.second = j->chr_pos.second;
                            example_map[tmp_aln->getName()].diff.second = tmp_aln->diff;
                            example_map[tmp_aln->getName()].chr = j->chr.second;
                        }

                    }

//                    std::cout << "Realn: " << j->chr.first << " " << j->chr_pos.first << " " << tmp_aln->getName() << endl;
                }


//                cout << "step3" << endl;
            } else if (i != bp_rlns.end()){
                if ((tmp_aln->getPosition() + tmp_aln->getRefLength() < i->chr_pos.first - distance
                && tmp_aln->getRefID() == i->chr_idx.first) || tmp_aln->getRefID() < i->chr_idx.first) // read is behind bp
                    mapped_file->Jump(i->chr_idx.first, max((long)0, i->chr_pos.first - distance));
                else if ((tmp_aln->getPosition() > i->chr_pos.first + distance
                && tmp_aln->getRefID() == i->chr_idx.first) || tmp_aln->getRefID() > i->chr_idx.first) { // read is ahead of  bp
                    while (!(tmp_aln->getPosition() <= i->chr_pos.first + distance
                            && tmp_aln->getRefID() == i->chr_idx.first)) {
                        i++;
                        if (tmp_aln->getRefID() < i->chr_idx.first)
                            break;
                    }
                }
            }

            num_reads++;
            if (num_reads % 10000 == 0) {
                cout << ref[tmp_aln->getRefID()].RefName << " " << tmp_aln->getPosition() << endl;
            }
        }
        tmp_aln = mapped_file->parseRead(Parameter::Instance()->min_mq);

    }
    for (auto example_map_tuple: example_map) {
        string name = example_map_tuple.first;
        if (name == "14a7a43d-d5d7-4bdf-98f5-09138c52e0d5")
            continue;
        map_example example = example_map_tuple.second;
        if (example.chr != "") {
            if (example.bases.size() == 0)
                continue;
            int read_pos_first = example.read_pos.first, read_pos_second = example.read_pos.second,
            ref_pos_first = example.ref_pos.first, ref_pos_second = example.ref_pos.second,
            diff_first = example.diff.first, diff_second = example.diff.second;
            cout << ">" << name << endl;
            cout << example.bases.substr(read_pos_first, read_pos_second - read_pos_first) << endl;
            string ref_str = bioio::read_fasta_contig(fasta, index.at(example.chr), ref_pos_first + diff_first,
                                                      ref_pos_second - ref_pos_first + diff_second);
            string seq_str = example.bases.substr(read_pos_first, read_pos_second - read_pos_first);
            vector<CigarOp> cigar = map_read_example(ref_str, seq_str);
            if (cigar.empty())
                continue;

            if (cigar[0].Type == 'D') {
                ref_pos_first = ref_pos_first + cigar[0].Length;
                cigar.erase(cigar.begin());
            } else if (cigar[0].Type == 'I') {
                read_pos_first = read_pos_first + cigar[0].Length;
                cigar.erase(cigar.begin());
            }

            if (cigar[cigar.size() - 1].Type == 'D') {
                ref_pos_second = ref_pos_second - cigar[cigar.size() - 1].Length;
                cigar.erase(cigar.begin() + cigar.size() - 1);
            } else if (cigar[cigar.size() - 1].Type == 'I') {
                read_pos_second = read_pos_second - cigar[cigar.size() - 1].Length;
                cigar.erase(cigar.begin() + cigar.size() - 1);
            }
            int flag = 0;
            if (!example.strand) flag = 16;

            bam_out << name << "\t" << flag << "\t" << example.chr << "\t" << ref_pos_first + diff_first + 1
                    << "\t" << 60 << "\t";
            cout << endl;
            int cigar_len = 0;
            for (auto cigar_str: cigar) {
                bam_out << cigar_str.Length << cigar_str.Type;
                cout << cigar_str.Length << cigar_str.Type;
                if (cigar_str.Type != 'D')
                    cigar_len += cigar_str.Length;
            }

            bam_out << "\t*\t" << 0 << "\t" << 0 << "\t"
                    << example.bases.substr(read_pos_first, read_pos_second - read_pos_first)
                    << "\t" << "*" << endl;
            cout << "\t" << cigar_len << " " << read_pos_second - read_pos_first << endl;

            cout << endl;
        }
    }


    bam_out.close();



    for (size_t i = 0; i < points_size; i++) { // its not nice, but I may alter the length of the vector within the loop.
        if (points[i]->get_SVtype() & TRA) {
            vector<Breakpoint *> new_points;

            detect_merged_svs(points[i]->get_coordinates(), ref, new_points);
            if (!new_points.empty()) {                            // I only allow for 1 split!!
                points[i] = new_points[0];
                points.push_back(new_points[1]);
            }
        }
    }


    std::string chr;

	//std::cout<<"fin up"<<std::endl;
	for (size_t i = 0; i < points.size(); i++) {
		if (points[i]->get_SVtype() & TRA) {
			points[i]->calc_support();
			points[i]->predict_SV();

		}
		if (points[i]->get_support() >= Parameter::Instance()->min_support && points[i]->get_length() > Parameter::Instance()->min_length) {
			printer->printSV(points[i]);
		}

        string chr_start, chr_stop;
        long pos_start = IPrinter::calc_pos(points[i]->get_coordinates().start.most_support, ref, chr_start);
        long pos_stop = IPrinter::calc_pos(points[i]->get_coordinates().stop.most_support, ref, chr_stop);
        std::cout << chr_start << " " << pos_start << " " << chr_stop << " " << pos_stop << endl;
        std::cout << "Number of support: " << points[i]->get_support() << endl;
        for (auto j :points[i]->get_coordinates().support){
            string chr_start, chr_stop;
            long pos_start = IPrinter::calc_pos(j.second.coordinates.first, ref, chr_start);
            long pos_stop = IPrinter::calc_pos(j.second.coordinates.second, ref, chr_stop);
            std::cout << j.first << " " << chr_start << " " << pos_start << " " << chr_stop << " " << pos_stop << endl;
        }
        std::cout << endl;
	}
	//std::cout<<"Done"<<std::endl;
	if (Parameter::Instance()->genotype) {
		fclose(ref_allel_reads);
	}
}


void add_events(Alignment *& tmp, std::vector<str_event> events, short type, long ref_space, IntervallTree & bst, TNode *&root, long read_id, bool add) {

	bool flag = (strcmp(tmp->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);
	for (size_t i = 0; i < events.size(); i++) {
		//	if (events[i - 1].length > Parameter::Instance()->min_segment_size || events[i].length > Parameter::Instance()->min_segment_size) {
		position_str svs;
		read_str read;
		if (events[i].is_noise) {
			read.type = 2;
		} else {
			read.type = 0;
		}
		read.SV = events[i].type;
		read.sequence = events[i].sequence;

		if (flag) {
			std::cout << "ADD EVENT " << tmp->getName() << " " << tmp->getRefID() << " " << events[i].pos << " " << abs(events[i].length) << std::endl;
		}
		svs.start.min_pos = (long) events[i].pos + ref_space;
		svs.stop.max_pos = svs.start.min_pos + events[i].length;

		if (tmp->getStrand()) {
			read.strand.first = (tmp->getStrand());
			read.strand.second = !(tmp->getStrand());
		} else {
			read.strand.first = !(tmp->getStrand());
			read.strand.second = (tmp->getStrand());
		}
		//	start.support[0].read_start.min = events[i].read_pos;

		read.read_strand.first = tmp->getStrand();
		read.read_strand.second = tmp->getStrand();
		if (flag) {
			std::cout << tmp->getName() << " " << tmp->getRefID() << " " << svs.start.min_pos << " " << svs.stop.max_pos << " " << svs.stop.max_pos - svs.start.min_pos << std::endl;
		}

		if (svs.start.min_pos > svs.stop.max_pos) {
			//can this actually happen?
			read.coordinates.first = svs.stop.max_pos;
			read.coordinates.second = svs.start.min_pos;
		} else {
			read.coordinates.first = svs.start.min_pos;
			read.coordinates.second = svs.stop.max_pos;
		}

		svs.start.max_pos = svs.start.min_pos;
		svs.stop.min_pos = svs.stop.max_pos;

		if (svs.start.min_pos > svs.stop.max_pos) { //incase they are inverted
			svs_breakpoint_str pos = svs.start;
			svs.start = svs.stop;
			svs.stop = pos;
			pair<bool, bool> tmp = read.strand;
			read.strand.first = tmp.second;
			read.strand.second = tmp.first;
		}

		//TODO: we might not need this:
		if (svs.start.min_pos > svs.stop.max_pos) {
			read.coordinates.first = svs.stop.max_pos;
			read.coordinates.second = svs.start.min_pos;
		} else {
			read.coordinates.first = svs.start.min_pos;
			read.coordinates.second = svs.stop.max_pos;
		}

		read.id = read_id;
		svs.support[tmp->getName()] = read;
		svs.support[tmp->getName()].length = events[i].length;
		Breakpoint * point = new Breakpoint(svs, events[i].length);
		if (add) {
			bst.insert_existant(point, root);
		} else {
			bst.insert(point, root);
		}
		//std::cout<<"Print:"<<std::endl;
		//bst.print(root);
	}
//	}
}

void add_splits(Alignment *& tmp, std::vector<aln_str> events, short type, RefVector ref, IntervallTree& bst, TNode *&root, long read_id, bool add) {
	bool flag = (strcmp(tmp->getName().c_str(), Parameter::Instance()->read_name.c_str()) == 0);

	if (flag) {
		cout << "SPLIT: " << std::endl;
		for (size_t i = 0; i < events.size(); i++) {
			std::cout << events[i].pos << " stop: " << events[i].pos + events[i].length << " " << events[i].RefID << " READ: " << events[i].read_pos_start << " " << events[i].read_pos_stop;
			if (events[i].strand) {
				cout << " +" << endl;
			} else {
				cout << " -" << endl;
			}
		}
	}

	for (size_t i = 1; i < events.size(); i++) {
		//	if (events[i - 1].length > Parameter::Instance()->min_segment_size || events[i].length > Parameter::Instance()->min_segment_size) {
		position_str svs;
		//position_str stop;
		read_str read;
		read.sequence = "NA";
		//read.name = tmp->getName();
		read.type = type;
		read.SV = 0;
		read.read_strand.first = events[i - 1].strand;
		read.read_strand.second = events[i].strand;

		//stop.support.push_back(read);
		if (events[i].RefID == events[i - 1].RefID) { //IF different chr -> tra
			if (events[i - 1].strand == events[i].strand) { //IF same strand -> del/ins/dup
				if (events[i - 1].strand) {
					read.strand.first = events[i - 1].strand;
					read.strand.second = !events[i].strand;
				} else {
					read.strand.first = !events[i - 1].strand;
					read.strand.second = events[i].strand;
				}
				//	int len1 = 0;
				//int len2 = 0;
				svs.read_start = events[i - 1].read_pos_stop; // (short) events[i - 1].read_pos_start + (short) events[i - 1].length;
				svs.read_stop = events[i].read_pos_start;
				if (events[i - 1].strand) {
					svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
				} else {
					svs.start.min_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
					svs.stop.max_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
				}

				if (flag) {
					cout << "Debug: SV_Size: " << (svs.start.min_pos - svs.stop.max_pos) << " tmp: " << (svs.stop.max_pos - svs.start.min_pos) << " Ref_start: " << svs.start.min_pos - get_ref_lengths(events[i].RefID, ref) << " Ref_stop: " << svs.stop.max_pos - get_ref_lengths(events[i].RefID, ref) << " readstart: " << svs.read_start << " readstop: " << svs.read_stop << std::endl;
				}

				if ((svs.stop.max_pos - svs.start.min_pos) > Parameter::Instance()->min_length * -1 && ((svs.stop.max_pos - svs.start.min_pos) + (Parameter::Instance()->min_length) < (svs.read_stop - svs.read_start) && (svs.read_stop - svs.read_start) > (Parameter::Instance()->min_length * 2))) {
					if (!events[i].cross_N || (double) ((svs.stop.max_pos - svs.start.min_pos) + Parameter::Instance()->min_length) < ((double) (svs.read_stop - svs.read_start) * Parameter::Instance()->avg_ins)) {
						svs.stop.max_pos += (svs.read_stop - svs.read_start); //TODO check!
						if (Parameter::Instance()->print_seq) {
							svs.read_stop = events[i].read_pos_start;
							svs.read_start = events[i - 1].read_pos_stop;
							if (svs.read_stop > tmp->getAlignment()->QueryBases.size()) {
								cerr << "BUG: split read ins! " << svs.read_stop << " " << tmp->getAlignment()->QueryBases.size() << " " << tmp->getName() << endl;
							}
							if (!events[i - 1].strand) {
								std::string tmp_seq = reverse_complement(tmp->getAlignment()->QueryBases);

								read.sequence = reverse_complement(tmp_seq.substr(svs.read_start, svs.read_stop - svs.read_start));
							} else {
								read.sequence = tmp->getAlignment()->QueryBases.substr(svs.read_start, svs.read_stop - svs.read_start);
							}
							if (flag) {
								cout << "INS: " << endl;
								cout << "split read ins! " << events[i - 1].read_pos_stop << " " << events[i].read_pos_start << " " << " " << tmp->getAlignment()->QueryBases.size() << " " << tmp->getName() << endl;
								cout << "Seq+:" << read.sequence << endl;
							}
						}
						read.SV |= INS;
					} else {
						read.SV |= 'n';
					}

				} else if ((svs.start.min_pos - svs.stop.max_pos) * -1 > (svs.read_stop - svs.read_start) + (Parameter::Instance()->min_length)) {
					if (!events[i].cross_N || (double) (svs.start.min_pos - svs.stop.max_pos) * Parameter::Instance()->avg_del * -1.0 > (double) ((svs.read_stop - svs.read_start) + (Parameter::Instance()->min_length))) {
						read.SV |= DEL;
						if (flag) {
							cout << "DEL2" << endl;
						}
					} else {
						read.SV |= 'n';
					}

				} else if ((svs.start.min_pos - svs.stop.max_pos) > Parameter::Instance()->min_length && (svs.read_start - svs.read_stop) < Parameter::Instance()->min_length) { //check with respect to the coords of reads!
					if (flag) {
						cout << "DUP: " << endl;
					}
					read.SV |= DUP;
				} else {
					if (flag) {
						cout << "N" << endl;
					}
					read.SV = 'n';
				}
			} else { // if first part of read is in a different direction as the second part-> INV

				read.strand.first = events[i - 1].strand;
				read.strand.second = !events[i].strand;

				bool is_overlapping = overlaps(events[i - 1], events[i]);
				if (is_overlapping && (events[i - 1].length > Parameter::Instance()->min_segment_size || events[i].length > Parameter::Instance()->min_segment_size)) {
					if (flag) {
						std::cout << "Overlap curr: " << events[i].pos << " " << events[i].pos + events[i].length << " prev: " << events[i - 1].pos << " " << events[i - 1].pos + events[i - 1].length << " " << tmp->getName() << std::endl;
					}
					read.SV |= NEST;

					if (events[i - 1].strand) {
						svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
						svs.stop.max_pos = (events[i].pos + events[i].length) + get_ref_lengths(events[i].RefID, ref);
					} else {
						svs.start.min_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
						svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
					}

					if (svs.start.min_pos > svs.stop.max_pos) {
						long tmp = svs.start.min_pos;
						svs.start.min_pos = svs.stop.max_pos;
						svs.stop.max_pos = tmp;
					}
				} else if (!is_overlapping) {
					read.SV |= INV;
					if (events[i - 1].strand) {
						svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
						svs.stop.max_pos = (events[i].pos + events[i].length) + get_ref_lengths(events[i].RefID, ref);
					} else {
						svs.start.min_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
						svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
					}
				}
			}

		} else { //if not on the same chr-> TRA
			read.strand.first = events[i - 1].strand;
			read.strand.second = !events[i].strand;
			if (events[i - 1].strand == events[i].strand) {
				//check this with + - strands!!

				if (events[i - 1].strand) { //"++"
					svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
				} else { //"--"
					svs.start.min_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
				}
			} else {
				if (events[i - 1].strand) { //"+-"
					svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
				} else { //"-+"
					svs.start.min_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
					svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
				}
			}
			read.SV |= TRA;
		}

		if (read.SV != 'n') {
			if (flag) {
				std::cout << "SPLIT: " << TRANS_type(read.SV) << " start: " << svs.start.min_pos - get_ref_lengths(events[i].RefID, ref) << " stop: " << svs.stop.max_pos - get_ref_lengths(events[i].RefID, ref);
				if (events[i - 1].strand) {
					std::cout << " +";
				} else {
					std::cout << " -";
				}
				if (events[i].strand) {
					std::cout << " +";
				} else {
					std::cout << " -";
				}
				std::cout << " " << tmp->getName() << std::endl;
				std::cout << "READ: " << svs.read_start << " " << svs.read_stop << " " << svs.read_start - svs.read_stop << std::endl;
			}
			//std::cout<<"split"<<std::endl;

			svs.start.max_pos = svs.start.min_pos;
			svs.stop.min_pos = svs.stop.max_pos;
			if (svs.start.min_pos > svs.stop.max_pos) {
				//maybe we have to invert the directions???
				svs_breakpoint_str pos = svs.start;
				svs.start = svs.stop;
				svs.stop = pos;

				pair<bool, bool> tmp = read.strand;

				read.strand.first = tmp.second;
				read.strand.second = tmp.first;
			}

			//TODO: we might not need this:
			if (svs.start.min_pos > svs.stop.max_pos) {
				read.coordinates.first = svs.stop.max_pos;
				read.coordinates.second = svs.start.min_pos;
			} else {
				read.coordinates.first = svs.start.min_pos;
				read.coordinates.second = svs.stop.max_pos;
			}

			//pool out?
			read.id = read_id;
			svs.support[tmp->getName()] = read;
			svs.support[tmp->getName()].length = abs(read.coordinates.second - read.coordinates.first);
			Breakpoint * point = new Breakpoint(svs, abs(read.coordinates.second - read.coordinates.first));
			//std::cout<<"split ADD: " << <<" Name: "<<tmp->getName()<<" "<< svs.start.min_pos- get_ref_lengths(events[i].RefID, ref)<<"-"<<svs.stop.max_pos- get_ref_lengths(events[i].RefID, ref)<<std::endl;
			if (add) {
				bst.insert_existant(point, root);
			} else {
				bst.insert(point, root);
			}
			//	std::cout<<"Print:"<<std::endl;
			//	bst.print(root);
		}
	}
	//}
}

void estimate_parameters(std::string read_filename) {
	if (Parameter::Instance()->skip_parameter_estimation) {
		return;
	}
	cout << "Estimating parameter..." << endl;
	BamParser * mapped_file = 0;
	RefVector ref;
	if (read_filename.find("bam") != string::npos) {
		mapped_file = new BamParser(read_filename);
		ref = mapped_file->get_refInfo();
	} else {
		cerr << "File Format not recognized. File must be a sorted .bam file!" << endl;
		exit(EXIT_FAILURE);
	}

	Alignment * tmp_aln = mapped_file->parseRead(Parameter::Instance()->min_mq);
	double num = 0;
	double avg_score = 0;
	double avg_mis = 0;
	double avg_indel = 0;
	double avg_diffs_perwindow = 0;
	vector<int> mis_per_window; //histogram over #differences
	vector<int> scores;
//	std::string curr, prev = "";
	double avg_dist = 0;
	double tot_avg_ins = 0;
	double tot_avg_del = 0;
	while (!tmp_aln->getQueryBases().empty() && num < 1000) {	//1000
		//	std::cout<<"test "<<tmp_aln->getName()<<std::endl;
		if (rand() % 100 < 20 && ((tmp_aln->getAlignment()->IsPrimaryAlignment()) && (!(tmp_aln->getAlignment()->AlignmentFlag & 0x800)))) {				//}&& tmp_aln->get_is_save()))) {
			//1. check differences in window => min_treshold for scanning!
			//2. get score ration without checking before hand! (above if!)
			double dist = 0;
			double avg_del = 0;
			double avg_ins = 0;
			vector<int> tmp = tmp_aln->get_avg_diff(dist, avg_del, avg_ins);
			//	std::cout<<"Debug:\t"<<avg_del<<" "<<avg_ins<<endl;
			tot_avg_ins += avg_ins;
			tot_avg_del += avg_del;
			//
			avg_dist += dist;
			double avg_mis = 0;
			for (size_t i = 0; i < tmp.size(); i++) {
				while (tmp[i] + 1 > mis_per_window.size()) { //adjust length
					mis_per_window.push_back(0);
				}
				avg_mis += tmp[i];
				mis_per_window[tmp[i]]++;
			}
			//	std::cout <<avg_mis/tmp.size()<<"\t";
			//get score ratio
			double score = round(tmp_aln->get_scrore_ratio());
			//	std::cout<<score<<"\t"<<std::endl;;
			if (score > -1) {
				while (score + 1 > scores.size()) {
					scores.push_back(0);
				}
				scores[score]++;
			}
			num++;
		}

		mapped_file->parseReadFast(Parameter::Instance()->min_mq, tmp_aln);
	}
	if (num == 0) {
		std::cerr << "Too few reads detected in " << Parameter::Instance()->bam_files[0] << std::endl;
		exit(EXIT_FAILURE);
	}
	vector<int> nums;
	size_t pos = 0;
	Parameter::Instance()->max_dist_alns = floor(avg_dist / num) / 2;
	Parameter::Instance()->window_thresh = 50;			//25;
	if (!mis_per_window.empty()) {
		for (size_t i = 0; i < mis_per_window.size(); i++) {

			for (size_t j = 0; j < mis_per_window[i]; j++) {
				nums.push_back(i);
			}
		}
		pos = nums.size() * 0.95; //the highest 5% cutoff
		if (pos > 0 && pos <= nums.size()) {
			Parameter::Instance()->window_thresh = std::max(Parameter::Instance()->window_thresh, nums[pos]); //just in case we have too clean data! :)
		}
		nums.clear();
	}

	for (size_t i = 0; i < scores.size(); i++) {
		for (size_t j = 0; j < scores[i]; j++) {
			nums.push_back(i);
		}
	}
	pos = nums.size() * 0.05; //the lowest 5% cuttoff
	Parameter::Instance()->score_treshold = 2; //nums[pos]; //prev=2

	//cout<<"test: "<<tot_avg_ins<<" "<<num<<endl;
	//cout<<"test2: "<<tot_avg_del<<" "<<num<<endl;
	Parameter::Instance()->avg_del = tot_avg_del / num;
	Parameter::Instance()->avg_ins = tot_avg_ins / num;

	std::cout << "\tMax dist between aln events: " << Parameter::Instance()->max_dist_alns << std::endl;
	std::cout << "\tMax diff in window: " << Parameter::Instance()->window_thresh << std::endl;
	std::cout << "\tMin score ratio: " << Parameter::Instance()->score_treshold << std::endl;
	std::cout << "\tAvg DEL ratio: " << Parameter::Instance()->avg_del << std::endl;
	std::cout << "\tAvg INS ratio: " << Parameter::Instance()->avg_ins << std::endl;

}

bool overlaps(aln_str prev, aln_str curr) {

	double ratio = 0;
	double overlap = 0;
	if (prev.pos + Parameter::Instance()->min_length < curr.pos + curr.length && prev.pos + prev.length - Parameter::Instance()->min_length > curr.pos) {
		overlap = min((curr.pos + curr.length), (prev.pos + prev.length)) - max(prev.pos, curr.pos);
		ratio = overlap / (double) min(curr.length, prev.length);
	}
//	std::cout<<overlap<<" "<<ratio<<std::endl;
	return (ratio > 0.4 && overlap > 200);
}


