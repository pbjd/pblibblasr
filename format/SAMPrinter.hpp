#ifndef _BLASR_FORMAT_SAMPRINTER_HPP_
#define _BLASR_FORMAT_SAMPRINTER_HPP_ 

#include <sstream>
#include <stdint.h>
#include "SMRTSequence.hpp"
#include "datastructures/alignment/AlignmentCandidate.hpp"
#include "datastructures/alignment/AlignmentContext.hpp"


#define MULTI_SEGMENTS 0x1
#define ALL_SEGMENTS_ALIGNED 0x2
#define SEGMENT_UNMAPPED 0x4
#define NEXT_SEGMENT_UNMAPPED 0x8
#define SEQ_REVERSED 0x10
#define SEQ_NEXT_REVERSED 0x20
#define FIRST_SEGMENT 0x40
#define LAST_SEGMENT 0x80
#define SECONDARY_ALIGNMENT 0x100
#define NO_PASS_QUALITY 0x200
#define PCR_OR_OPTICAL_DUPLICATE 0x400

namespace SAMOutput {

enum Clipping {hard, soft, none};

void BuildFlag(T_AlignmentCandidate &alignment, AlignmentContext &context, uint16_t &flag); 

//
// Trimming is used for both hard non-clipping
// so it is called trim instead of clip.
//
void CreateDNAString(DNASequence &seq, DNASequence &clippedSeq, 
    int trimFront=0, int trimEnd=0); 

void AddGaps(T_AlignmentCandidate &alignment, int gapIndex,
        vector<int> &opSize, vector<char> &opChar); 

void CreateNoClippingCigarOps(T_AlignmentCandidate &alignment, 
        vector<int> &opSize, vector<char> &opChar); 
//
// 
// The aligned sequence is either the sequence from the first
// aligned base to the last (hard and no clipping), or first high
// quality base to the last high quality base (soft clipping).
//
template<typename T_Sequence>
void SetAlignedSequence(T_AlignmentCandidate &alignment, T_Sequence &read,
    T_Sequence &alignedSeq, Clipping clipping = none); 

template<typename T_Sequence>
void SetSoftClip(T_AlignmentCandidate &alignment, T_Sequence &read, 
    DNALength &softClipPrefix, DNALength &softClipSuffix); 

template<typename T_Sequence>
void SetHardClip(T_AlignmentCandidate &alignment, T_Sequence &read, 
    int &prefixClip, int &suffixClip); 

void CigarOpsToString(vector<int> &opSize, vector<char> &opChar, 
        string &cigarString); 

//
// Straight forward: create the cigar string allowing some clipping
// The read is provided to give length and hq information.
//
template<typename T_Sequence>
void CreateCIGARString(T_AlignmentCandidate &alignment, T_Sequence &read,
        string &cigarString, Clipping clipping=none); 

template<typename T_Sequence>
void PrintAlignment(T_AlignmentCandidate &alignment, T_Sequence &read,
        ostream &samFile, AlignmentContext &context, Clipping clipping = none,
        int subreadIndex = 0, int nSubreads = 0); 
}
#endif // _BLASR_FORMAT_SAMPRINTER_HPP_
