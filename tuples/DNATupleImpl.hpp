#include <cassert>
#include "NucConversion.hpp"

template<typename Sequence> 
int SearchSequenceForTuple(Sequence &seq, TupleMetrics &tm, DNATuple &queryTuple) {
    DNALength p;
    PositionDNATuple tempTuple, upperTuple;

    p = 0;
    DNALength cur = 0;
    DNALength curValidEnd = 0;

    //
    // Construct the mask-off bit pair for the shifted tuple.
    //
    PositionDNATuple maskLeftTuple;
    maskLeftTuple.tuple = 3;
    maskLeftTuple.tuple = maskLeftTuple.tuple << 2*tm.tupleSize;
    maskLeftTuple.tuple = ~maskLeftTuple.tuple;
    PositionDNATuple testTuple;
    while (curValidEnd < seq.length) {
        //
        // Search for the next available window that can be translated into a tuple.
        //
        cur = curValidEnd;
        while(curValidEnd < seq.length and IsACTG[seq.seq[curValidEnd]]) {
            curValidEnd++;
        }
        if (curValidEnd - cur >= tm.tupleSize) {
            //
            // Found a span that does not have N's in it, 
            //
            assert (tempTuple.FromStringRL(&(seq.seq[cur]), tm) == 1);
            p = cur;
            if (tempTuple.tuple == queryTuple.tuple) {
                return 1;
            }
            for (p++; p < curValidEnd - tm.tupleSize + 1; p++) {
                tempTuple.tuple >>=2;
                //				tempTuple.tuple &= maskLeftTuple.tuple;
                upperTuple.tuple = TwoBit[seq.seq[p+tm.tupleSize-1]];
                upperTuple.tuple = upperTuple.tuple << (2 * (tm.tupleSize-1));
                tempTuple.tuple += upperTuple.tuple;
                if (tempTuple.tuple == queryTuple.tuple) {
                    return 1;
                }
            }
        }
        else {
            ++curValidEnd;
        }
    }
}


template<typename Sequence>
int SequenceToTupleList(Sequence &seq, TupleMetrics &tm, TupleList<DNATuple> &tupleList) {
    DNALength p;
    PositionDNATuple tempTuple, upperTuple;

    p = 0;
    DNALength cur = 0;
    DNALength curValidEnd = 0;

    //
    // Construct the mask-off bit pair for the shifted tuple.
    //
    PositionDNATuple maskLeftTuple;
    maskLeftTuple.tuple = 3;
    maskLeftTuple.tuple = maskLeftTuple.tuple << 2*tm.tupleSize;
    maskLeftTuple.tuple = ~maskLeftTuple.tuple;
    PositionDNATuple testTuple;
    while (curValidEnd < seq.length) {
        //
        // Search for the next available window that can be translated into a tuple.
        //
        cur = curValidEnd;
        while(curValidEnd < seq.length and IsACTG[seq.seq[curValidEnd]]) {
            curValidEnd++;
        }
        if (curValidEnd - cur >= tm.tupleSize) {
            //
            // Found a span that does not have N's in it, 
            //
            assert (tempTuple.FromStringRL(&(seq.seq[cur]), tm) == 1);
            p = cur;
            tupleList.Append(tempTuple);
            for (p++; p < curValidEnd - tm.tupleSize + 1; p++) {
                tempTuple.tuple >>=2;
                //				tempTuple.tuple &= maskLeftTuple.tuple;
                upperTuple.tuple = TwoBit[seq.seq[p+tm.tupleSize-1]];
                upperTuple.tuple = upperTuple.tuple << (2 * (tm.tupleSize-1));
                tempTuple.tuple += upperTuple.tuple;
                //testTuple.FromStringRL(&seq.seq[p], tm);
                //assert(testTuple.tuple == tempTuple.tuple);
                tupleList.Append(tempTuple);
            }
        }
        else {
            ++curValidEnd;
        }
    }
    return tupleList.size();
}


template<typename Sequence>
int SequenceToTupleList(Sequence &seq, TupleMetrics &tm, TupleList<PositionDNATuple> &tupleList) {
    DNALength p;
    PositionDNATuple tempTuple, upperTuple;

    p = 0;
    DNALength cur = 0;
    DNALength curValidEnd = 0;

    //
    // Construct the mask-off bit pair for the shifted tuple.
    //
    PositionDNATuple maskLeftTuple;
    maskLeftTuple.tuple = 3;
    maskLeftTuple.tuple = maskLeftTuple.tuple << 2*tm.tupleSize;
    maskLeftTuple.tuple = ~maskLeftTuple.tuple;
    PositionDNATuple testTuple;
    while (curValidEnd < seq.length) {
        //
        // Search for the next available window that can be translated into a tuple.
        //
        cur = curValidEnd;
        while(curValidEnd < seq.length and IsACTG[seq.seq[curValidEnd]]) {
            curValidEnd++;
        }
        if (curValidEnd - cur >= tm.tupleSize) {
            //
            // Found a span that does not have N's in it, 
            //
            assert (tempTuple.FromStringRL(&(seq.seq[cur]), tm) == 1);
            p = cur;
            tempTuple.pos = p;
            tupleList.Append(tempTuple);
            for (p++; p < curValidEnd - tm.tupleSize + 1; p++) {
                tempTuple.tuple >>=2;
                //				tempTuple.tuple &= maskLeftTuple.tuple;
                upperTuple.tuple = TwoBit[seq.seq[p+tm.tupleSize-1]];
                upperTuple.tuple = upperTuple.tuple << (2 * (tm.tupleSize-1));
                tempTuple.tuple += upperTuple.tuple;
                tempTuple.pos = p;
                //testTuple.FromStringRL(&seq.seq[p], tm);
                //assert(testTuple.tuple == tempTuple.tuple);
                tupleList.Append(tempTuple);
            }
        }
        else {
            ++curValidEnd;
        }
    }
    return tupleList.size();
}
