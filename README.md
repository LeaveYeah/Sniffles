# Sniffles
Sniffles is a structural variation caller using third generation sequencing (PacBio or Oxford Nanopore). It detects all types of SVs (10bp+) using evidence from split-read alignments, high-mismatch regions, and coverage analysis. Please note the current version of Sniffles requires sorted output from BWA-MEM (use -M and -x parameter) or NGMLR with the optional SAM attributes enabled! If you experience problems or have suggestions please contact: fritz.sedlazeck@gmail.com


Please see our github wiki for more information (https://github.com/fritzsedlazeck/Sniffles/wiki)


# How to build Sniffles
<pre>wget https://github.com/fritzsedlazeck/Sniffles/archive/master.tar.gz -O Sniffles.tar.gz
tar xzvf Sniffles.tar.gz
cd Sniffles-master/
mkdir -p build/
cd build/
cmake ..
make

cd ../bin/sniffles*
./sniffles</pre>

Note Mac users often have to provide parameters to the cmake command:
<pre>cmake -D CMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.7 -D CMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.7 .. 
</pre>


**************************************
## NGMLR
Sniffles performs best with the mappings of NGMLR our novel long read mapping method. 
Please see:
https://github.com/philres/ngmlr

****************************************
## Citation:
Please see and cite our paper:
https://www.nature.com/articles/s41592-018-0001-7
  
**************************************
## Poster & Talks:

[Accurate and fast detection of complex and nested structural variations using long read technologies](http://schatzlab.cshl.edu/presentations/2016/2016.10.28.BIODATA.PacBioSV.pdf)
Biological Data Science, Cold Spring Harbor Laboratory, Cold Spring Harbor, NY, 26 - 29.10.2016

[NGMLR: Highly accurate read mapping of third generation sequencing reads for improved structural variation analysis](http://www.cibiv.at/~philipp_/files/gi2016_poster_phr.pdf) 
Genome Informatics 2016, Wellcome Genome Campus Conference Centre, Hinxton, Cambridge, UK, 19.09.-2.09.2016

**************************************
## Datasets used in the mansucript:
We provide the NGMLR aligned reads and the Sniffles calls for the data sets used:  

Arabidopsis trio: 
+ [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Arabidopsis_trio](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Arabidopsis_trio) . 

Genome in the Bottle trio: 
+ Mappings: [ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/) . 

+ SV calls: [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/GiaB/](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/GiaB/)

NA12878: 
+ [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/NA12878/](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/NA12878/) .  

SKBR3: 
+ [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Skbr3/](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Skbr3/) .

CMake Options:
-DCMAKE_MODULE_PATH="/data2/junwenwang/m204333/Project/seqan-v2.2.0/util/cmake" -DSEQAN_INCLUDE_PATH="/data2/junwenwang/m204333/Project/seqan-v2.2.0/include"

GCC Path:
/projects/bsi/gentools/src/gcc/v4.9.4/bin/gcc

Scratch:
//through reads


void add_realign_read(BreakPointRealign bp_realign){
    read_str tmp;
    if (bp_realign.coordinate.first < bp_realign.coordinate.second) {
        tmp.coordinates.first = bp_realign.coordinate.first + diff;
        if (bp_realign.isSameStrand)
            tmp.coordinates.second = bp_realign.coordinate.second + diff;
        else tmp.coordinates.second = bp_realign.coordinate.second - diff;
    } else {
        tmp.coordinates.second = bp_realign.coordinate.first + diff;
        if (bp_realign.isSameStrand)
            tmp.coordinates.first = bp_realign.coordinate.second + diff;
        else tmp.coordinates.first = bp_realign.coordinate.second - diff;
    }
    tmp.SV = TRA;
    tmp.type = 1;

//                        std::cout << bp_realign.bp->get_coordinates().support.size() << endl;
    auto map = bp_realign.bp->get_coordinates().support;
    if (map.find(tmp_aln->getName()) == map.end())
        bp_realign.bp->add_read(tmp, tmp_aln->getName());
    else
        bp_realign.bp->add_read(tmp, tmp_aln->getName() + "_extra");}
}
bool detect_gaps(Alignment* tmp_aln, BreakPointRealign bp_realign, int &diff) {
    std::vector <aln_str> split_events = tmp_aln->getSA(ref);
    for (auto i: split_events) {
        if (i->isMain) {
            if (abs(i->pos - bp_realign.chr_pos.first) < distance) {
                diff = i->pos - bp_realign.chr_pos.first;
                tmp_aln->high_error_side = false;
                if (i->strand) tmp_aln->bp_read_pos = i->read_pos_start;
                else tmp_aln->bp_read_pos = i->read_pos_stop;
                return true;
            } else if (abs(i->pos + i->length - bp_realign.chr_pos.first) < distance) {
                diff = i->pos + i->length - bp_realign.chr_pos.first;
                tmp_aln->high_error_side = true;
                if (i->strand) tmp_aln->bp_read_pos = i->read_pos_stop;
                else tmp_aln->bp_read_pos = i->read_pos_start;
                return true;
            }
        }
    }
    return false;
}


void realign_read(BreakPointRealign bp_realign, vector<aln_str> &event_aln, Alignment* tmp_aln,
        const bioio::FastaIndex  index, std::ifstream & fasta) {
    int distance = min(100, Parameter::Instance()->min_length);
    int diff;
    if (bp_realign.chr_pos.first - distance <= tmp_aln->getPosition() ||
        bp_realign.chr_pos.first + distance >= tmp_aln->getPosition() + tmp_aln->getRefLength()) {
        auto map = bp_realign.bp->get_coordinates().support;
        if (map.find(tmp_aln->getName()) != map.end()) return;

    } else {

        bool exists_high_error_side = cal_high_error_side(event_aln, bp_realign.chr_pos.first, distance, diff, tmp_aln);
        if (!exists_high_error_side) return;
        int alt_aln_score = map_read(tmp_aln, bp_realign, diff, distance, index, fasta);
        if (alt_aln_score - tmp_aln->high_error_region_score > distance / 5 &&
            alt_aln_score > -0.2 * distance) {
            add_realign_read(bp_realign);
        }
    }
}


//
int diff;
bool exists_high_error_side = cal_high_error_side(event_aln, j->chr_pos.first, distance, diff,
                                                  tmp_aln);
//                    cout << "step2.2" << endl;
if (!exists_high_error_side) continue;
int alt_aln_score = map_read(tmp_aln, *j, diff, distance, index, fasta);
if (alt_aln_score - tmp_aln->high_error_region_score > distance / 5 &&
alt_aln_score > -0.2 * distance) {
read_str tmp;
if (j->coordinate.first < j->coordinate.second) {
tmp.coordinates.first = j->coordinate.first + diff;
if (j->isSameStrand)
tmp.coordinates.second = j->coordinate.second + diff;
else tmp.coordinates.second = j->coordinate.second - diff;
} else {
tmp.coordinates.second = j->coordinate.first + diff;
if (j->isSameStrand)
tmp.coordinates.first = j->coordinate.second + diff;
else tmp.coordinates.first = j->coordinate.second - diff;
}
tmp.SV = TRA;
tmp.type = 1;

//                        std::cout << j->bp->get_coordinates().support.size() << endl;
auto map = j->bp->get_coordinates().support;
if (map.find(tmp_aln->getName()) ==
map.end())
j->bp->add_read(tmp, tmp_aln->getName());
else
j->bp->add_read(tmp, tmp_aln->getName()+"_extra");

//                        std::cout << j->bp->get_coordinates().start.max_pos << " " << j->bp->get_coordinates().stop.max_pos << endl;
//                        std::cout << j->bp->get_coordinates().support.size() << endl;
//                        for (auto i: j->bp->get_coordinates().support){
//                            std::cout << i.first << " " << i.second.coordinates.first << " "
//                             << i.second.coordinates.second << endl;
//                        }
}

if (j->chr_idx.first != tmp_aln->getRefID() || j->chr_pos.first - distance < tmp_aln->getPosition() ||
j->chr_pos.first + distance < tmp_aln->getPosition()) {
