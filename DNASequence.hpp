#ifndef  _BLASR_DNA_SEQUENCE_HPP_
#define  _BLASR_DNA_SEQUENCE_HPP_

#include <stdint.h>

typedef uint32_t DNALength;
typedef unsigned char Nucleotide;

class DNASequence {
public:
    DNALength length;
    Nucleotide *seq;
    int bitsPerNuc;
    bool deleteOnExit;

    DNASequence();

    //--- functions ---//
    
    DNALength size();

    void TakeOwnership(DNASequence &rhs);

    void Append(const DNASequence &rhs, DNALength appendPos=0);

    DNASequence &Copy(const DNASequence &rhs, DNALength rhsPos=0, DNALength rhsLength=0);

    void ShallowCopy(const DNASequence &rhs);

    int GetStorageSize();

    DNASequence &operator=(const DNASequence &rhs);

    void Print(std::ostream &out, int lineLength = 50);

    void PrintSeq(std::ostream &out, int lineLength = 50);

    void Allocate(DNALength plength);

    void ReferenceSubstring(const DNASequence &rhs, int pos=0, int substrLength=0); 

    DNALength MakeRCCoordinate(DNALength forPos );

    void CopyAsRC(DNASequence &rc, DNALength pos=0, DNALength rcLength =0);

    void MakeRC(DNASequence &rc, DNALength pos=0, DNALength rcLength=0);

    void ToTwoBit();

    void ToThreeBit();

    void ToFourBit();

    void ConvertThreeBitToAscii();

    void ToAscii();

    void Assign(DNASequence &ref, DNALength start=0, DNALength plength=0);

    void ToLower();

    void ToUpper(); 

    void Concatenate(const Nucleotide *moreSeq, DNALength moreSeqLength); 

    std::string GetTitle() const; 

    void Concatenate(const Nucleotide* moreSeq);

    void Concatenate(DNASequence &seq); 

    int Compare(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length); 

    int LessThanEqual(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length); 

    int Equals(DNASequence &rhs, DNALength rhsPos, DNALength length, DNALength pos=0 ); 

    int LessThan(DNALength pos,  DNASequence &rhs, DNALength rhsPos, DNALength length); 

    void CleanupASCII(); 

    Nucleotide operator[](int i) {
        return seq[i];
    }

    Nucleotide GetNuc(DNALength i); 

    DNALength GetRepeatContent(); 

    void CleanupOnFree();

    void FreeIfControlled(); 

    virtual void Free(); 

    void Resize(DNALength newLength);

    DNALength GetSeqStorage();
};

template<typename T>
DNALength ResizeSequence(T &dnaseq, DNALength newLength); 

#endif // _BLASR_DNA_SEQUENCE_HPP_
