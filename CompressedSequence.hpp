#ifndef _BLASR_COMPRESSED_SEQUENCE_HPP_
#define _BLASR_COMPRESSED_SEQUENCE_HPP_

typedef unsigned char CompressedNucleotide;

template<typename T_Sequence>
class CompressedSequence : public FASTASequence {
private:
    int hasIndex;
    int hasTitle;
    ReverseCompressIndex index;

    static const unsigned char MaskCount = 0xf;
    static const unsigned char MaskNuc   = 0xf0;
    static const unsigned char ShiftCount = 4;

public:
    //
    // This is just a placeholder for now.  
    // No extra data here, just the ability to decompress.  Right now 
    // the utilities for the compressed dna sequences
    // are in CompressedSeqUtils.h, which could move here later.
    //
    QualityValue *qual;

    void CopyConfiguration(CompressedSequence<T_Sequence> &rhs); 

    void ShallowCopy(CompressedSequence<T_Sequence> &rhs); 

    void MakeRC(CompressedSequence &rc); 

    Nucleotide operator[](DNALength i); 

    Nucleotide GetNuc(DNALength i); 

    unsigned char GetCount(DNALength i); 

    char *GetName();

    void Copy(FASTASequence &rhs); 

    float GetAverageQuality(); 

    void SortHomopolymerQualities(); 

    CompressedSequence(); 

    void SetHasTitle(); 

    void SetHasIndex(); 

    void Write(std::string outFileName); 

    void Read(std::string inFileName); 

    int BuildFourBitReverseIndex(int binSize); 

    int BuildReverseIndex(int maxRun, int binSize); 

    long Lookup4BitCompressedSequencePos(int cpPos); 

    int LookupSequencePos(int cpPos); 

    char GetCount(unsigned char ch); 

    DNALength FourBitCompressHomopolymers(); 

    static int Only4BitACTG(CompressedNucleotide *seq, int seqLength); 

    int Only4BitACTG(); 

    void RemoveCompressionCounts(); 

    DNALength FourBitDecompressHomopolymers(int start, int end,
            T_Sequence &decompSeq); 

    DNALength CondenseHomopolymers();
};



#endif // _BLASR_COMPRESSED_SEQUENCE_HPP_
