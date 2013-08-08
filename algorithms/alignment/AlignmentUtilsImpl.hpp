#include <cassert>
#include <vector>
#include <stdint.h>
#include <iostream>
#include <ostream>
#include <cstdlib>
#include <fstream>
#include "Types.h"
#include "defs.h"
#include "utils.hpp"
#include "NucConversion.hpp"
#include "DNASequence.hpp"
#include "datastructures/matrix/FlatMatrix.hpp"
#include "datastructures/alignment/Alignment.hpp"
#include "algorithms/alignment/AlignmentUtils.hpp"

using namespace blasr;

template<typename T_QuerySequence, typename T_TargetSequence, typename T_ScoreFn>
int ComputeAlignmentScore(blasr::Alignment &alignment,
        T_QuerySequence &query,
        T_TargetSequence &text,
        T_ScoreFn &scoreFn) {
    VectorIndex b, q, t, l, bi, gi;
    int alignmentScore = 0;

    for (b = 0; b < alignment.blocks.size(); b++ ) {
        for ((q = alignment.qPos + alignment.blocks[b].qPos,
                    t = alignment.tPos + alignment.blocks[b].tPos,
                    l = 0) ; 
                l < alignment.blocks[b].length; q++, t++, l++) {
            alignmentScore += scoreFn.Match(text, t, query, q);
        }
        if (alignment.gaps.size() == alignment.blocks.size() + 1) {
            for (gi = 0; gi < alignment.gaps[b+1].size(); gi++) {
                if (alignment.gaps[b+1][gi].seq == Gap::Target) {
                    alignmentScore += alignment.gaps[b+1][gi].length * scoreFn.ins;
                }
                else {
                    alignmentScore += alignment.gaps[b+1][gi].length * scoreFn.del;
                }
            }
        }
    }
    return alignmentScore;

}


template<typename T_QuerySequence, typename T_TargetSequence>
int ComputeAlignmentScore(blasr::Alignment &alignment,
        T_QuerySequence &query,
        T_TargetSequence &text,
        int matchScores[5][5],
        int ins,
        int del) {
    VectorIndex b;
    int alignmentScore = 0;
    VectorIndex q, t;
    DNALength l;

    //
    // Compute the average match score for the case that
    // there are unassigned alignments.
    //
    int totalMismatchPenalty = 0, averageMismatchPenalty;
    for (q = 0; q < 4; q++ ) {
        for (t = 0; t < 4; t++ ){ 
            if (q != t) totalMismatchPenalty += matchScores[q][t];
        }
    }

    averageMismatchPenalty = (int) (totalMismatchPenalty/ 12.0);

    for (b = 0; b < alignment.blocks.size(); b++ ) {
        for ((q = alignment.qPos + alignment.blocks[b].qPos,
                    t = alignment.tPos + alignment.blocks[b].tPos,
                    l = 0) ; 
                l < alignment.blocks[b].length; q++, t++, l++) {
            alignmentScore += matchScores[ThreeBit[query.seq[q]]][ThreeBit[text.seq[t]]];
        }
        if (alignment.blocks.size() > 0 and b < alignment.blocks.size() - 1) {
            // 
            // Add the gap penalties.
            //
            int queryGap = (alignment.blocks[b+1].qPos - 
                    (alignment.blocks[b].qPos +
                     alignment.blocks[b].length));
            int textGap  = (alignment.blocks[b+1].tPos - 
                    (alignment.blocks[b].tPos + 
                     alignment.blocks[b].length));

            int commonGap = queryGap;
            if (textGap < queryGap) {
                commonGap = textGap;
            }

            queryGap -= commonGap;
            textGap  -= commonGap;
            alignmentScore += averageMismatchPenalty * commonGap;

            alignmentScore += (queryGap * ins);
            alignmentScore += (textGap * del);
        }
    }
    return alignmentScore;
}

template<typename T_Alignment>
inline void PrintCompareSequencesAlignmentStats(T_Alignment &alignment, std::ostream &out) {
    int lastBlock;

    lastBlock = alignment.blocks.size() -1;

    int qLength, tLength;
    if (lastBlock >= 0) {
        qLength = (alignment.blocks[lastBlock].qPos 
                + alignment.blocks[lastBlock].length) ;
        tLength =  (alignment.blocks[lastBlock].tPos 
                + alignment.blocks[lastBlock].length);
    }
    else {
        qLength = tLength = 0;
    }

    //
    // First print the query 
    //
    int alignmentQStart;
    if (lastBlock >= 0) {
        alignmentQStart = alignment.qPos + alignment.qAlignedSeqPos;
    }
    else {
        alignmentQStart = alignment.qAlignedSeqPos;
    }
    out << alignment.qName
        << " " << alignment.qLength
        << " " << alignmentQStart
        << " " << alignmentQStart + qLength;

    if (alignment.qStrand == 0) {
        out << " + ";
    }
    else {
        out << " - ";
    }

    int alignmentTStart;
    if (lastBlock >= 0) {
        alignmentTStart = alignment.tPos + alignment.tAlignedSeqPos;
    }
    else {
        alignmentTStart = 0;
    }
    out << " " << alignment.tName 
        << " " << alignment.tLength;
    if (alignment.tStrand == 0) {
        out << " " << alignmentTStart
            << " " << alignmentTStart + tLength;
        out << " + ";
    }
    else {
        out	<< " " << alignment.tLength - (alignmentTStart + tLength)
            << " " << alignment.tLength - (alignmentTStart);
        out << " - ";
    }

    out << alignment.score 
        << " " << alignment.nMatch 
        << " " << alignment.nMismatch
        << " " << alignment.nIns
        << " " << alignment.nDel
        << " " << (int) alignment.mapQV 
        << " ";
}

template<typename T_Alignment>
inline int ReadCompareSequencesAlignmentStats(std::istream &in, T_Alignment &alignment) {
    int qEnd, tEnd;
    int qLength, tLength;
    char qStrand, tStrand;
    if (!(in >> alignment.qName )) return 0;
    if (!(in >> qLength)) return 0;
    if (!(in >> alignment.qPos)) return 0;
    if (!(in >> qEnd)) return 0;
    if (!(in >> qStrand)) return 0;
    if (!(in >> alignment.tName)) return 0;
    if (!(in >> tLength)) return 0;
    if (!(in >> alignment.tPos)) return 0;
    if (!(in >> tEnd)) return 0;
    if (!(in >> tStrand)) return 0;
    if (!(in >> alignment.score)) return 0;
    if (!(in >> alignment.nMatch)) return 0;
    if (!(in >> alignment.nMismatch)) return 0;
    if (!(in >> alignment.nIns)) return 0;
    if (!(in >> alignment.nDel)) return 0;
    return 1;
}


template<typename T_Alignment>
inline int ReadCompSeqAlignment(std::istream &in, T_Alignment &alignment) {
    if (!ReadCompareSequencesAlignmentStats(in, alignment)) return 0;
    std::string alignStr;
    if (!(in >> alignment.qString)) return 0;
    if (!(in >> alignStr)) return 0;
    if (!(in >> alignment.tString)) return 0;
    std::string eol;
    std::getline(in, eol);
    return 1;
}


template<typename T_QuerySequence, typename T_TargetSequence>
void  AppendGapCharacters(Gap &gap, 
        T_QuerySequence &query, T_TargetSequence &text, 
        DNALength &q, DNALength &t,
        char mismatchChar, char gapChar,
        std::string &textStr, std::string &alignStr, std::string &queryStr) {
    int gp;
    for (gp = 0; gp < gap.length; gp++) {
        if (gap.seq == Gap::Query) {
            textStr.push_back(text[t]);
            alignStr.push_back(mismatchChar);
            queryStr.push_back(gapChar);
            t++;
        }
        else if (gap.seq == Gap::Target) {
            textStr.push_back(gapChar);
            alignStr.push_back(mismatchChar);
            queryStr.push_back(query[q]);
            q++;
        }
    }
}

template<typename T_Alignment, typename T_QuerySequence, typename T_TargetSequence>
void CreateAlignmentStrings(T_Alignment &alignment, 
        T_QuerySequence &query, T_TargetSequence &text, 
        std::string &textStr, std::string &alignStr, std::string &queryStr, 
        DNALength queryLength, DNALength textLength) {
    DNALength q = alignment.qPos;
    DNALength t = alignment.tPos;
    DNALength qPos, tPos;
    DNALength  g;
    char mismatchChar = '*';
    char matchChar = '|';
    char gapChar = '-';
    char gapSeparationChar = ' ';
    if (alignment.blocks.size() == 0) {
        textStr = "";
        alignStr = "";
        queryStr = "";
        return;
    }

    if (alignment.gaps.size() == 0) {
        //
        // If there is no gap list, add the gaps as an offset here.
        //
        if (alignment.blocks[0].qPos > 0 or
                alignment.blocks[0].tPos > 0) {
            // commonGapLen should be the shorter gap.
            qPos = alignment.blocks[0].qPos;
            tPos = alignment.blocks[0].tPos;
            DNALength commonGapLen = qPos;
            if (commonGapLen > tPos) {
                commonGapLen = tPos;
            }
            for (g = 0; g < commonGapLen; g++ ) {
                textStr.push_back(text[t]);
                alignStr.push_back(mismatchChar);
                queryStr.push_back(query[q]);
                t++;
                q++;
            }
            tPos -= commonGapLen;
            qPos -= commonGapLen;

            //
            // one of tPos or qPos is now 0.
            // The other represents extra sequence
            // that should be output before starting the alignment
            DNALength p;
            for (p = 0; p < tPos; p++) {
                textStr.push_back(text[t]);
                alignStr.push_back(gapSeparationChar);
                queryStr.push_back(gapChar);
                t++;
            }
            for (p = 0; p < qPos; p++) {
                textStr.push_back(gapChar);
                alignStr.push_back(gapSeparationChar);
                queryStr.push_back(query[q]);
                q++;
            }
        }
    }

    //
    // Add gap characters if they are before the beginning of the alignment. 
    // This shouldn't happen, but for some local alignments, it can.
    //
    DNALength b, bl, gi;
    if (alignment.gaps.size() > 0) {
        // The first gap is before the first block of characters.
        DNALength gi;
        for (gi = 0; gi < alignment.gaps[0].size(); gi++) {
            AppendGapCharacters(alignment.gaps[0][gi], query, text, q, t, mismatchChar, gapChar, textStr, alignStr, queryStr);
        }
    }

    for (b = 0; b < alignment.size() ; b++) {
        for (bl = 0; bl < alignment.blocks[b].length; bl++ ) {
            queryStr.push_back(query[q]);
            textStr.push_back(text[t]);
            assert(queryLength == 0 or q < queryLength);
            assert(textLength == 0 or t < textLength);
            if (TwoBit[query[q]] != 
                    TwoBit[text[t]])
                alignStr.push_back(mismatchChar);
            else
                alignStr.push_back(matchChar);
            q++;
            t++;
        }
        //
        //  There are no gaps to count after the last block, so 
        //  don't add the gapped characters for this.
        //
        if (alignment.blocks.size() == 0)
            continue;
        if (b == alignment.blocks.size() - 1) {
            continue;
        }
        if (alignment.gaps.size() > 0) {
            for (gi = 0; gi < alignment.gaps[b+1].size(); gi++) {
                AppendGapCharacters(alignment.gaps[b+1][gi], query, text, q, t, gapSeparationChar, gapChar, textStr, alignStr, queryStr);
            }
        }
        else {


            DNALength queryGapLen = (alignment.blocks[b+1].qPos - 
                    alignment.blocks[b].qPos - alignment.blocks[b].length);
            DNALength textGapLen  = (alignment.blocks[b+1].tPos - 
                    alignment.blocks[b].tPos - alignment.blocks[b].length);

            if (queryGapLen > 0 or textGapLen > 0) {
                // commonGapLen should be the shorter gap.
                DNALength commonGapLen = queryGapLen; 
                if (queryGapLen > textGapLen) {
                    commonGapLen = textGapLen;
                }
                textGapLen -= commonGapLen;
                queryGapLen -= commonGapLen;

                for (g = 0; g < queryGapLen; g++, q++ ){
                    textStr.push_back(gapChar);
                    alignStr.push_back(gapSeparationChar);
                    queryStr.push_back(query[q]);
                }
                for (g = 0; g < textGapLen; g++, t++ ){
                    textStr.push_back(text[t]);
                    alignStr.push_back(gapSeparationChar);
                    queryStr.push_back(gapChar);

                }

                for (g = 0; g < commonGapLen; g++ ) {
                    textStr.push_back(text[t]);
                    alignStr.push_back(gapSeparationChar);
                    queryStr.push_back(query[q]);
                    t++;
                    q++;
                }
            }
        }
    }
}

template<typename T_Alignment>
void ComputeAlignmentStats(T_Alignment &alignment, Nucleotide* qSeq, Nucleotide *tSeq, int matchMatrix[5][5], int ins, int del) {
    int qp = 0, tp = 0;
    int nMatch = 0, nMismatch = 0, nIns =0, nDel = 0;
    float pctSimilarity;
    std::string textStr, alignStr, queryStr;
    CreateAlignmentStrings(alignment, qSeq, tSeq, textStr, alignStr, queryStr);
    int i;
    int alignLength = textStr.size();

    for (i = 0; i < alignLength; i++ ) {
        if ((textStr[i] != '-') and (queryStr[i] != '-')) {
            int ti = (int)textStr[i];
            int qi = (int)queryStr[i];
            if (ThreeBit[ti] == ThreeBit[qi]) {
                nMatch++;
            }
            else {
                nMismatch++;
            }
            tp++;
            qp++;
        }
        else if (textStr[i] == '-' and queryStr[i] != '-') {
            nIns++;
            qp++;
        }
        else if (queryStr[i] == '-' and textStr[i] != '-') {
            nDel++;
            tp++;
        }
    }
    if (tp + qp > 0) {
        if (textStr.size() + queryStr.size() > 0) {
            pctSimilarity = (nMatch*2.0) / (textStr.size() + queryStr.size()) * 100;
        }
        else {
            pctSimilarity = 0;
        }
    }
    else {
        pctSimilarity = 0;
    }

    alignment.score = ComputeAlignmentScore(queryStr, textStr, matchMatrix, ins, del);
    alignment.nMatch = nMatch;
    alignment.nMismatch = nMismatch;
    alignment.nDel = nDel;
    alignment.nIns = nIns;
    alignment.pctSimilarity = pctSimilarity;
}

template<typename T_Alignment>
int ComputeDrift(T_Alignment &alignment) {
    VectorIndex b;
    int qGap = 0, tGap = 0, commonGap = 0;
    int drift = 0;
    int maxDrift = 0;
    int driftBetweenBlocks;
    if (alignment.blocks.size() == 0)
        return 0;
    for (b = 0; b < alignment.blocks.size() - 1; b++) { 
        driftBetweenBlocks = ComputeDrift(alignment.blocks[b], alignment.blocks[b+1]);
        drift += driftBetweenBlocks;
        if (abs(drift) > maxDrift) 
            maxDrift = abs(drift);
    }
    return maxDrift;
}

template<typename T_Alignment>
void RemoveAlignmentPrefixGaps(T_Alignment &alignment) {
    if (alignment.gaps.size()  == 0) {
        return;
    }

    unsigned int g;
    int tStart = 0, qStart = 0;
    for (g = 0; g < alignment.gaps[0].size(); g++) {
        if (alignment.gaps[0][g].seq == Gap::Target) {
            qStart += alignment.gaps[0][g].length;
        }
        else if (alignment.gaps[0][g].seq == Gap::Query) {
            tStart += alignment.gaps[0][g].length;
        }
    }
    int b;
    int nBlocks;
    for (b = 0, nBlocks = alignment.blocks.size(); b < nBlocks; b++) {
        alignment.blocks[b].qPos -= qStart;
        alignment.blocks[b].tPos -= tStart;
    }
    alignment.gaps[0].clear();
    alignment.tPos += tStart;
    alignment.qPos += qStart;
}
