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


```
void map_read_example(string ref_str, String seq_str){
    using namespace seqan;
    typedef String<char> TSequence;
    typedef StringSet <TSequence> TStringSet;
    typedef Align <TSequence, ArrayGaps> TAlign;// container for strings
    typedef StringSet <TSequence, Dependent<>> TDepStringSet;
    typedef seqan::Alignment <TDepStringSet> TAlignStringSet; // dependent string set
    typedef Iterator<TRow>::Type TRowIterator;

    TAlign align;
    TRowIterator it = begin(row1);
    TRowIterator itEnd = end(row1);
    for (; it != itEnd; ++it)
    {
        if (isGap(it))
            std::cout << gapValue<char>();
        else
            std::cout << *it;
    }
    transform(ref_str.begin(), ref_str.end(), ref_str.begin(), ::toupper);
    TSequence reference = ref_str;

    transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
    TSequence sequence = seq_str;

    resize(rows(align), 2);
    assignSource(row(align, 0), reference);
    assignSource(row(align, 1), sequence);

    TRow & row1 = row(align, 0);
    TRow & row2 = row(align, 1);

}

int map_read(Alignment  * tmp_aln, BreakPointRealign bp, int diff, int distance,
             const bioio::FastaIndex  index, std::ifstream & fasta) {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    using namespace seqan;
    typedef String<char> TSequence;
    typedef StringSet <TSequence> TStringSet;
    typedef Align <TSequence, ArrayGaps> TAlign;// container for strings
    typedef StringSet <TSequence, Dependent<>> TDepStringSet;
    typedef seqan::Alignment <TDepStringSet> TAlignStringSet; // dependent string set
    typedef Graph <TAlignStringSet> TAlignGraph;       // alignment graph// sequence type


    if (bp.isSameStrand) bp.chr_pos.second += diff;
    else bp.chr_pos.second -= diff;

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

    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);

    TAlignGraph alignG(sequences);

    int score = globalAlignment(alignG, Score<int, Simple>(0, -1, -1), AlignConfig<false, false, true, true>(),
                                LinearGaps());

    return score;
}
```