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

int ComputeAlignmentScore(
        std::string &queryStr, std::string &textStr,
        int matchScores[5][5], int ins, int del) { 
    if (queryStr.size() != textStr.size()) {
        std::cout << "Computing alignment score using invalid alignment string." << std::endl;
        std::cout << "Bailing out."<<std::endl;
        exit(1);
    }
    VectorIndex i;
    int score = 0;
    int alignStrLen = queryStr.size();
    for(i=0; i < alignStrLen; i++) {
        if (queryStr[i] != '-' and
                textStr[i] != '-') {
            score += matchScores[ThreeBit[(int)queryStr[i]]][ThreeBit[(int)textStr[i]]];
        }
        else {
            if (queryStr[i] == '-' and textStr[i] != '-') {
                score += del;
            }
            else if (queryStr[i] != '-' and textStr[i] == '-') {
                score += ins;
            }
            else {
                score += matchScores[4][4];
            }
        }
    }
    return score;
}

inline int GetNumberWidth(unsigned int value) {
    // 0 has a width of 1.
    int width = 1;
    while (value / 10 > 0)  {
        value = value / 10;
        ++width;
    }
    return width;
}

/*
 * This should be changed to read any type of alignment, since templates are being
 * used.
 */
inline void ReadCompSeqAlignments(std::string &compSeqAlignmentFileName, std::vector<CompSeqAlignment> &alignments) {
    std::ifstream in;
    CrucialOpen(compSeqAlignmentFileName, in);
    CompSeqAlignment alignment;
    while (ReadCompSeqAlignment(in, alignment)) {
        alignments.push_back(alignment);
    }
}

inline void PrintAlignmentStats(Alignment &alignment, std::ostream &out) {
    out << "    nMatch: " << alignment.nMatch << std::endl;
    out << " nMisMatch: " << alignment.nMismatch << std::endl;
    out << "      nIns: " << alignment.nIns << std::endl;
    out << "      nDel: " << alignment.nDel << std::endl;
    out << "      %sim: " << alignment.pctSimilarity << std::endl;
    out << "     Score: " << alignment.score << std::endl; 
}

int ComputeDrift(Block &cur,  Block &next ) {

    int tGap  = (next.tPos - cur.TEnd());
    int qGap = (next.qPos - cur.QEnd());

    int commonGap = 0;

    if (tGap > 0 and qGap > 0) {
        commonGap = abs(tGap - qGap);
    }
    tGap -= commonGap;
    qGap -= commonGap;

    return tGap - qGap;
}
