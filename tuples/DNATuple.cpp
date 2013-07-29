#include <cassert>
#include <vector>
#include <stdint.h>
#include <ostream>
#include <string>
#include "Types.h"
#include "DNASequence.hpp"
#include "SeqUtils.hpp"
#include "NucConversion.hpp"
#include "tuples/TupleMetrics.hpp"
#include "tuples/TupleList.hpp"
#include "tuples/TupleOperations.h"
#include "tuples/DNATuple.hpp"

int DNATuple::FromStringLR(Nucleotide *strPtr, TupleMetrics &tm) {
    DNASequence tmpSeq;
    tmpSeq.seq = strPtr;
    tmpSeq.length = tm.tupleSize;
    if (!OnlyACTG(tmpSeq))
        return 0;

    if (tm.tupleSize == 0)
        return 1;

    tuple = 0;
    Nucleotide *p;
    Nucleotide *endPtr = &strPtr[tm.tupleSize - 1];
    for (p = strPtr; p != endPtr; p++) {
        // If it is not possible to convert this string, return null.
        if (ThreeBit[*p] > 3) {
            return 0;
        }
        tuple += TwoBit[*p];
        tuple <<=2;
    }
    //
    // The tuple size is guaranteed to be at least 
    // 1, so it's safe to add the last value.
    // This cannot be in the previous loop since
    // the shift shouldn't happen.
    tuple += TwoBit[*p];
    return 1;
}

int DNATuple::FromStringRL(Nucleotide *strPtr, TupleMetrics &tm) {

    //
    // Tuples are created with the right-most character
    // in the most significant bit in the sequence.
    //
    DNASequence tmpSeq;
    tmpSeq.seq = strPtr;
    tmpSeq.length = tm.tupleSize;
    if (!OnlyACTG(tmpSeq))
        return 0;

    if (tm.tupleSize == 0)
        return 1;

    tuple = 0;
    Nucleotide *p;
    for (p = strPtr + tm.tupleSize - 1; p > strPtr; p--) {
        tuple += TwoBit[*p];
        tuple <<=2;
    }
    //
    // The tuple size is guaranteed to be at least 
    // 1, so it's safe to add the last value.
    // This cannot be in the previous loop since
    // the shift shouldn't happen.
    tuple += TwoBit[*p];
    return 1;
}

int DNATuple::ShiftAddRL(Nucleotide nuc, TupleMetrics &tm) {
    if (ThreeBit[nuc] > 3) {
        return 0;
    }
    else {
        tuple >>= 2;
        tuple += (TwoBit[nuc] << ((tm.tupleSize-1)*2));
        return 1;
    }
}

int DNATuple::MakeRC(DNATuple &dest, TupleMetrics &tm) {
    int i;
    ULong tempTuple = tuple;
    dest.tuple = 0;
    ULong mask = 0x3;
    if (tm.tupleSize == 0)
        return 0;
    for (i = 0; i < tm.tupleSize - 1; i++ ) {
        dest.tuple += (~tempTuple & mask);
        tempTuple >>= 2;
        dest.tuple <<= 2;
    }
    dest.tuple += (~tempTuple & mask);
    return 1;
}

std::string DNATuple::ToString(TupleMetrics &tm) {
    int i;
    std::string s;
    ULong tempTuple = tuple;
    for (i = 0;i < tm.tupleSize;i++) {
        s.insert(s.begin(), TwoBitToAscii[tempTuple & 3]);
        tempTuple = tempTuple >> 2;
    }
    return s;
}

bool 
CompareByTuple::operator()(const DNATuple &lhs, const DNATuple &rhs) 
const {
    return lhs.tuple < rhs.tuple;
}

int PositionDNATuple::operator==(const DNATuple &pTuple) const {
    return tuple == pTuple.tuple;
}

PositionDNATuple::PositionDNATuple() : DNATuple() {
    pos = -1;
}

PositionDNATuple::PositionDNATuple(PositionDNATuple &tupleP, DNALength posP) {
    tuple = tupleP.tuple;
    pos   = posP;
}


int 
OrderPositionDNATuplesByPosition::operator()(
    const PositionDNATuple &lhs, const PositionDNATuple &rhs) 
const {
    return lhs.pos < rhs.pos;
}

int
OrderPositionDNATuplesByTuple::operator()(
    const PositionDNATuple &lhs, const PositionDNATuple &rhs) 
const {
    return lhs.tuple < rhs.tuple;
}
