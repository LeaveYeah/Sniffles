//
// Created by wuyei on 9/15/2019.
//

#include "Realign.h"
#include "limits.h"
#include "seqan/align.h"
#include "seqan/bam_io.h"


using namespace std;

void ClippedRead::setChrId(int RefId){
    this->RefId = RefId;
}
void ClippedRead::setBases(string bases){
    this->bases = bases;
}
void ClippedRead::set_ref_pos(vector<aln_str> split_events){
    long max_pos=0, min_pos=LONG_MAX;
    for (aln_str split_event: split_events){
        if (split_event.RefID != this->RefId)
            continue;
        if (split_event.pos < min_pos)
            min_pos = split_event.pos;
        if (split_event.pos + split_event.length > max_pos)
            max_pos = split_event.pos + split_event.length;
    }
    this->ref_start = min_pos;
    this->ref_stop = max_pos;
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

bool cal_high_error_side(vector<differences_str> &event_aln, long pos, long distance, Alignment * tmp_aln) {
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


int map_clipped_read(rln_str event_rln, int distance, const bioio::FastaIndex  index, std::ifstream & fasta, string bases) {
    using namespace seqan;
    typedef String<char> TSequence;
    typedef StringSet<TSequence> TStringSet;
    typedef Align<TSequence, ArrayGaps> TAlign;// container for strings
    typedef StringSet<TSequence, Dependent<>> TDepStringSet;
    typedef seqan::Alignment<TDepStringSet> TAlignStringSet; // dependent string set
    typedef Graph<TAlignStringSet> TAlignGraph;       // alignment graph// sequence type

    if (event_rln.isSameStrand) event_rln.ref_pos.second += event_rln.diff;
    else event_rln.ref_pos.second -= event_rln.diff;

    if (!event_rln.side) event_rln.read_pos -= distance; //lefthand side
    if (event_rln.isSameStrand != event_rln.side) event_rln.ref_pos.second -= distance;

    string ref_str = bioio::read_fasta_contig(fasta, index.at(event_rln.ref.second), event_rln.ref_pos.second, distance);
    transform(ref_str.begin(), ref_str.end(), ref_str.begin(), ::toupper);
    TSequence reference = ref_str;
    if (event_rln.read_pos + distance > bases.size())
        return -50;

    if (event_rln.read_pos < 0)
        return -50;
//    std::cout << bp.chr.second << " " << bp.chr_pos.second << endl;
//    std::cout << tmp_aln->bp_read_pos << " " << distance << endl;
    string seq_str = bases.substr(event_rln.read_pos, distance);
    transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
    TSequence sequence = seq_str;
    TStringSet sequences;
    appendValue(sequences, reference);
    appendValue(sequences, sequence);

    TAlignGraph alignG(sequences);


    int score = globalAlignment(alignG, Score<int, Simple>(0, -1, -1), AlignConfig<false, false, true, true>(), LinearGaps());
//    std::cout << score << endl;
//    std::cout << alignG << endl;
    return score;
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

void add_realign_read(string name, rln_str event_rln){
    read_str tmp;
    if (event_rln.coordinate.first < event_rln.coordinate.second) {
        tmp.coordinates.first = event_rln.coordinate.first + event_rln.diff;
        tmp.coordinates.second = event_rln.coordinate.second + event_rln.diff;
    } else {
        tmp.coordinates.second = event_rln.coordinate.first + event_rln.diff;
        if (event_rln.isSameStrand)
            tmp.coordinates.first = event_rln.coordinate.second + event_rln.diff;
        else tmp.coordinates.first = event_rln.coordinate.second - event_rln.diff;
    }
    tmp.SV = TRA;
    tmp.type = 1;

//                        std::cout << bp_realign.bp->get_coordinates().support.size() << endl;
    auto map = event_rln.bp->get_coordinates().support;
    if (map.find(name+"_rln_0") == map.end())
        event_rln.bp->add_read(tmp, name+"_rln_0");
    else
        event_rln.bp->add_read(tmp, name + "_rln_1");
}


bool detect_gap(Alignment* tmp_aln, long ref_pos, int distance, RefVector ref) {
    if (tmp_aln->getName() == "cb169fa3-e5f9-4007-b234-089110b6f71e")
        tmp_aln->bp_read_pos += 0;
//    cout << ref_pos << endl;
//    cout << "step2" << endl;
    std::vector <aln_str> split_events = tmp_aln->getSA(ref);
    if (split_events.size() == 0) return false;
    long read_length;
    aln_str last_event = split_events[split_events.size()-1];
//    cout << "step2.1" << endl;
    if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800) )
        read_length = tmp_aln->getAlignment()->Length;
    else if (last_event.strand) {
        read_length = last_event.read_pos_stop;
        if (last_event.cigar.size() > 0 && last_event.cigar.back().Type == 'S')
            read_length += last_event.cigar.back().Length;
    } else {
        read_length = last_event.read_pos_stop;
        if (last_event.cigar.size() > 0 && last_event.cigar.front().Type == 'S')
            read_length += last_event.cigar.front().Length;
    }
//    cout << "step2.2" << endl;
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
//    cout << "step2.2" << endl;
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

void detect_clipped_reads_rln(BpRln bp_realign,  Alignment* tmp_aln,
                              RefVector ref, std::map<std::string, ClippedRead> &mapClippedRead){

    int distance = min(100, Parameter::Instance()->min_length);
    std::vector <aln_str> split_events = tmp_aln->getSA(ref);
    if (tmp_aln->getName() == "cb169fa3-e5f9-4007-b234-089110b6f71ef")
        tmp_aln->bp_read_pos += 0;
    if (bp_realign.chr_pos.first - distance <= tmp_aln->getPosition() ||
        bp_realign.chr_pos.first + distance >= tmp_aln->getPosition() + tmp_aln->getRefLength()) {
        auto map_support = bp_realign.bp->get_coordinates().support;
        if (map_support.find(tmp_aln->getName()) != map_support.end()) return;
        bool existsGap = detect_gap(tmp_aln, bp_realign.chr_pos.first, distance,  ref);
        if (!existsGap) return;
        rln_str event_rln;
        event_rln.side = tmp_aln->high_error_side;
        event_rln.ref = bp_realign.chr;
        event_rln.ref_pos = bp_realign.chr_pos;
        event_rln.coordinate = bp_realign.coordinate;
        event_rln.bp = bp_realign.bp;
        event_rln.isSameStrand = bp_realign.isSameStrand;
//        if (tmp_aln->high_error_side)
        event_rln.read_pos = tmp_aln->bp_read_pos;
//        else event_rln.read_pos = tmp_aln->bp_read_pos - distance;
        event_rln.diff = tmp_aln->diff;

        if (mapClippedRead.find(tmp_aln->getName()) == mapClippedRead.end()){
            ClippedRead clipped_read;
            clipped_read.setChrId(tmp_aln->getRefID());
            clipped_read.set_ref_pos(split_events);
            if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800) )
                clipped_read.setBases(tmp_aln->getQueryBases());
            clipped_read.events_rln.push_back(event_rln);
            mapClippedRead[tmp_aln->getName()] = clipped_read;
        } else {
            if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800) )
                mapClippedRead[tmp_aln->getName()].setBases(tmp_aln->getQueryBases());
            mapClippedRead[tmp_aln->getName()].events_rln.push_back(event_rln);
        }
    }
}

void realign_across_read(BpRln bp_realign, vector<differences_str> &event_aln, Alignment* tmp_aln,
                         const bioio::FastaIndex  index, std::ifstream & fasta, RefVector ref) {
    int distance = min(100, Parameter::Instance()->min_length);
    if (bp_realign.chr_pos.first - distance > tmp_aln->getPosition() &&
        bp_realign.chr_pos.first + distance < tmp_aln->getPosition() + tmp_aln->getRefLength()) {

        bool exists_high_error_side = cal_high_error_side(event_aln, bp_realign.chr_pos.first, distance,  tmp_aln);
        if (!exists_high_error_side) return;
        int alt_aln_score = map_read(tmp_aln, bp_realign,  distance, index, fasta);
        if (alt_aln_score - tmp_aln->high_error_region_score > distance / 5 &&
            alt_aln_score > -0.2 * distance) {
            add_realign_read(bp_realign,  tmp_aln);
        }
    }
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

