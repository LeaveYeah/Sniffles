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
int start = IPrinter::calc_pos(SV->get_coordinates().start.most_support, ref, chr);

void add_map(long pos, std::map<long, int> map){
    if (map.find(pos) == map.end() map[pos] = 1;
    else map[pos]++;
}
void store_tra_pos(vector<tra_str> &positions, read_str read, std::string read_name) {
    for (size_t i = 0; i < positions.size(); i++) {
        if (abs(positions[i].position - read.coordinates.first) < Parameter::Instance()->min_length) {
            positions[i].hits++;
            positions[i].names.push_back(read_name);
            add_map(read.coordinates.first, positions[i].starts)
            add_map(read.coordinates.second, positions[i].stops)
            if (read.read_strand.first == read.read_strand.second)
                positions[i].sameStrand_hits++;
            else positions[i].diffStrand_hits++;
            return;
        }
    }
    tra_str tmp;
    tmp.position = read.coordinates.first;
    tmp.hits = 1;
    tmp.names.push_back(read_name);
    add_map(read.coordinates.first, tmp.starts)
    add_map(read.coordinates.second, tmp.stops)
    if (read.read_strand.first == read.read_strand.second)
        tmp.sameStrand_hits = 1;
    else tmp.diffStrand_hits = 1;
    positions.push_back(tmp);
}

void detect_merged_svs(position_str point, RefVector ref, vector<Breakpoint *> & new_points) {
    new_points.clear(); //just in case!
    vector<tra_str> pos_start;
    for (std::map<std::string, read_str>::iterator i = point.support.begin(); i != point.support.end(); ++i) {
        store_tra_pos(pos_start, (*i).second, (*i).first);
    }

    int start_count = 0;
    for (size_t i = 0; i < pos_start.size(); i++) {
        //std::cout<<pos_start[i].hits <<",";
        if ((pos_start[i].hits > Parameter::Instance()->min_support / 2)
                && (pos_start[i].hits <= Parameter::Instance()->min_support)) {
            //do realignment, determine the start and the end the realignment

        }
    }
    int stop_count = 0;
    for (size_t i = 0; i < pos_stop.size(); i++) {
        //	std::cout << pos_stop[i].hits << ",";
        if (pos_stop[i].hits > Parameter::Instance()->min_support) {
            stop_count++;
        }
    }
    if (stop_count > 1 || start_count > 1) {
        std::cout << "\tprocessing merged TRA" << std::endl;
        if (start_count > 1) {
            new_points.push_back(split_points(pos_start[0].names, point.support));
            new_points.push_back(split_points(pos_start[1].names, point.support));
        } else {
            new_points.push_back(split_points(pos_stop[0].names, point.support));
            new_points.push_back(split_points(pos_stop[1].names, point.support));
        }
    }
}


