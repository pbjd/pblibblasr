#include <cassert>
#include <vector>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdint.h>
#include "Types.h"
#include "PlatformId.h"
#include "NucConversion.hpp"
#include "DNASequence.hpp"
#include "Enumerations.h"
#include "ReverseCompressIndex.hpp"
#include "FASTASequence.hpp"
#include "FASTQSequence.hpp"
#include "CompressedSequence.hpp"

using namespace std;

//
// Initialize a read with quality probabilities from one with quality values.
//
int FASTQSequence::charToQuality = 33;

QVScale FASTQSequence::GetQVScale() {
    return qvScale;
}

void FASTQSequence::SetQVScale(QVScale qvScaleP) {
    qvScale                   = qvScaleP;
    qual.qvScale              = qvScale;
    deletionQV.qvScale        = qvScale;
    preBaseDeletionQV.qvScale = qvScale;
    insertionQV.qvScale       = qvScale;
    substitutionQV.qvScale    = qvScale;
    mergeQV.qvScale           = qvScale;
}

int FASTQSequence::GetStorageSize() {
    int total = 0;
    int nQV = 0;
    int nTag =0;
    if (!qual.Empty()) {
        nQV++;
    }
    if (!deletionQV.Empty()) { 
        nQV++;
    }
    if (!preBaseDeletionQV.Empty()) {
        nQV+=4;
    }
    if (!insertionQV.Empty()) {
        nQV++;
    }
    if (!substitutionQV.Empty()) {
        nQV++;
    }
    if (!mergeQV.Empty()) {
        nQV++;
    }
    if (deletionTag != NULL) {
        nTag++;
    }
    if (substitutionTag !=NULL) {
        nTag++;
    }
    total = nQV*sizeof(QualityValue)*length + nTag*sizeof(Nucleotide)*length;
    return total + FASTASequence::GetStorageSize();
}

FASTQSequence::FASTQSequence() : FASTASequence() {
    deletionTag       = NULL;
    substitutionTag   = NULL;

    //
    // For now assume a prior distribution to be the variation of the human genome.
    //
    deletionQVPrior = 0.001; 
    insertionQVPrior = 0.001; 
    substitutionQVPrior = 0.001;
    preBaseDeletionQVPrior = 0.001;

    subreadStart = subreadEnd = 0;
    qvScale = PHRED;
}

QualityValue FASTQSequence::GetDeletionQV(DNALength pos) {
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    if (deletionQV.Empty()) {
        return deletionQVPrior;
    }
    else {
        return deletionQV[pos];
    }
}

QualityValue FASTQSequence::GetMergeQV(DNALength pos) {
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    if (mergeQV.Empty()) {
        return 0;
    }
    else {
        return mergeQV[pos];
    }
}

Nucleotide FASTQSequence::GetSubstitutionTag(DNALength pos) {
    if (substitutionTag == NULL) {
        return 'N';
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return substitutionTag[pos];
}

Nucleotide FASTQSequence::GetDeletionTag(DNALength pos) {
    if (deletionTag == NULL) {
        return 'N';
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return deletionTag[pos];
}

QualityValue FASTQSequence::GetInsertionQV(DNALength pos) {
    if (insertionQV.Empty()) {
        return insertionQVPrior;
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return insertionQV[pos];
}

QualityValue FASTQSequence::GetSubstitutionQV(DNALength pos) {
    if (substitutionQV.Empty()) {
        return substitutionQVPrior;
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return substitutionQV[pos];
}

QualityValue FASTQSequence::GetPreBaseDeletionQV(DNALength pos, Nucleotide nuc) {
    if (preBaseDeletionQV.Empty()) {
        return preBaseDeletionQVPrior;
    }
    assert(pos < ((unsigned int)-1));
    assert(pos < length);
    return preBaseDeletionQV[pos*4 + TwoBit[nuc]];
}

void FASTQSequence::ShallowCopy(const FASTQSequence &rhs) {
    qual.ShallowCopy(rhs.qual);
    FASTASequence::ShallowCopy(rhs);
}

void FASTQSequence::ReferenceSubstring(const FASTQSequence &rhs) {
    ReferenceSubstring(rhs, 0, rhs.length);
}

void FASTQSequence::ReferenceSubstring(const FASTQSequence &rhs, DNALength pos) {
    ReferenceSubstring(rhs, pos, rhs.length - pos);
}

void FASTQSequence::ReferenceSubstring(const FASTQSequence &rhs, DNALength pos, DNALength substrLength) {
    SetQVScale(rhs.qvScale);
    if (substrLength == 0) {
        substrLength = rhs.length - pos;
    }
    FASTASequence::ReferenceSubstring(rhs,pos,substrLength);
    if (rhs.qual.Empty() == false) {
        qual.ShallowCopy(rhs.qual, pos);
    }
    if (rhs.deletionQV.Empty() == false) {
        deletionQV.ShallowCopy(rhs.deletionQV, pos);
    }
    if (rhs.mergeQV.Empty() == false) {
        mergeQV.ShallowCopy(rhs.mergeQV, pos);
    }
    if (rhs.insertionQV.Empty() == false) {
        insertionQV.ShallowCopy(rhs.insertionQV, pos);
    }
    if (rhs.preBaseDeletionQV.Empty() == false ){
        preBaseDeletionQV.ShallowCopy(rhs.preBaseDeletionQV, pos);
    }
    if (rhs.deletionTag != NULL) {
        deletionTag = &rhs.deletionTag[pos];
    }
    if (rhs.substitutionTag != NULL) {
        substitutionTag = &rhs.substitutionTag[pos];
    }
    if (rhs.substitutionQV.Empty() == false) {
        substitutionQV.ShallowCopy(rhs.substitutionQV, pos);
    }
    deletionQVPrior = rhs.deletionQVPrior;
    insertionQVPrior = rhs.insertionQVPrior;
    substitutionQVPrior = rhs.substitutionQVPrior;
    preBaseDeletionQVPrior = rhs.preBaseDeletionQVPrior;
}

void FASTQSequence::ClearAndNull(QualityValue *value) {
    if (value != NULL) {
        delete[] value;
    }
    value = NULL;
}
void FASTQSequence::CopyQualityValues(const FASTQSequence &rhs) {
    SetQVScale(rhs.qvScale);
    qual.Copy(rhs.qual, rhs.length);
    deletionQV.Copy(rhs.deletionQV, rhs.length);
    insertionQV.Copy(rhs.insertionQV, rhs.length);
    substitutionQV.Copy(rhs.substitutionQV, rhs.length);
    mergeQV.Copy(rhs.mergeQV, rhs.length);
    //
    // Handle the tags separtely (and verbosely)
    //
    if (rhs.deletionTag) {
        AllocateDeletionTagSpace(rhs.length);
        memcpy(deletionTag, rhs.deletionTag, sizeof(Nucleotide)*rhs.length);
    }
    else { 
        ClearAndNull(deletionTag);
    }

    if (rhs.substitutionTag) {
        AllocateSubstitutionTagSpace(rhs.length);
        memcpy(substitutionTag, rhs.substitutionTag, sizeof(Nucleotide)*rhs.length);
    }
    else {
        ClearAndNull(substitutionTag);
    }
}

void FASTQSequence::AllocateQualitySpace(DNALength qualLength) {
    qual.Allocate(qualLength);
}

void FASTQSequence::AllocateDeletionQVSpace(DNALength qualLength) {
    deletionQV.Allocate(qualLength);
}

void FASTQSequence::AllocateMergeQVSpace(DNALength len) {
    mergeQV.Allocate(len);
}

void FASTQSequence::AllocateDeletionTagSpace(DNALength qualLength) {
    if (deletionTag != NULL) delete[] deletionTag;
    deletionTag = new Nucleotide[qualLength];
}

void FASTQSequence::AllocatePreBaseDeletionQVSpace(DNALength qualLength) {
    preBaseDeletionQV.Allocate(qualLength);
}

void FASTQSequence::AllocateInsertionQVSpace(DNALength qualLength) {
    insertionQV.Allocate(qualLength);
}

void FASTQSequence::AllocateSubstitutionQVSpace(DNALength qualLength ){ 
    substitutionQV.Allocate(qualLength);
}

void FASTQSequence::AllocateSubstitutionTagSpace(DNALength qualLength ){ 
    if (substitutionTag != NULL) delete[] substitutionTag;
    substitutionTag = new Nucleotide[qualLength];
}

void FASTQSequence::AllocateRichQualityValues(DNALength qualLength) {
    AllocateDeletionQVSpace(qualLength);
    AllocateDeletionTagSpace(qualLength);
    AllocatePreBaseDeletionQVSpace(qualLength);
    AllocateInsertionQVSpace(qualLength);
    AllocateSubstitutionQVSpace(qualLength);
    AllocateSubstitutionTagSpace(qualLength);
    AllocateMergeQVSpace(qualLength);
}

void FASTQSequence::Copy(const FASTQSequence &rhs) {
    FASTASequence::Copy(rhs);
    CopyQualityValues(rhs);
}

FASTQSequence& FASTQSequence::operator=(const FASTQSequence &rhs) {
    this->Copy(rhs);
    return *this;
}

FASTQSequence::FASTQSequence(const FASTQSequence &rhs) {
    substitutionTag = NULL;
    deletionTag = NULL;
    this->Copy(rhs);
}

void FASTQSequence::Assign(FASTQSequence &rhs) {
    // copy the nucleotide part
    FASTASequence::Assign(rhs);
    // copy the qual part
    CopyQualityValues(rhs);
    SetQVScale(rhs.qvScale);
}

void FASTQSequence::PrintFastq(ostream &out, int lineLength) {
    PrintSeq(out, lineLength, '@');
    if (lineLength == 0) { 
        out << endl;
    }
    PrintFastqQuality(out, lineLength);
    if (lineLength == 0) {
        out << endl;
    }
}

void FASTQSequence::PrintFastqQuality(ostream &out, int lineLength) {
    out << "+" << endl;
    PrintAsciiQual(out, lineLength);
}

void FASTQSequence::PrintAsciiQual(ostream &out, int lineLength) {
    int i;
    if (lineLength == 0) {
        for (i = 0; i < length; i++) {
            out << (char)(qual[i]+charToQuality);
        }
    }
    else {
        for (i = 0; i < length; i++) {
            out << (char)(qual[i]+charToQuality);
            assert(lineLength != 0);
            if (i > 0 and (i+1) % lineLength==0) {
                out << endl;
            }
        }
        if (i == 0 or i % lineLength != 0) {
            out << endl;
        }
    }
}


void FASTQSequence::PrintQual(ostream &out, int lineLength) {
    out << ">" << this->title << endl;
    DNALength i;
    for (i = 0; i < length; i++ ){
        out << (int) qual[i];
        if (i > 0 and (i+1) % lineLength == 0)
            out << endl;
        else 
            out << " ";
    }
    if (i == 0 or i % lineLength != 0) {
        out << endl;
    }
}

void FASTQSequence::PrintQualSeq(ostream &out, int lineLength) {
    FASTASequence::PrintSeq(out, lineLength);
    int i;
    lineLength /= 4;
    PrintQual(out, lineLength);
}

void FASTQSequence::MakeRC(FASTQSequence &rc) {
    rc.SetQVScale(qvScale);
    FASTASequence::MakeRC(rc);
    if (qual.Empty() == true) {
        // there is no quality values, don't make an rc.
        return;
    }

    if (rc.qual.Empty() == false) {
        rc.qual.Free();
    }
    rc.AllocateQualitySpace(length);
    int i;
    for (i = 0; i < length; i++ ){
        rc.qual.data[length - i - 1] = qual[i];
    }

    if (deletionQV.Empty() == false) {
        //
        // The read contains rich quality values. Reverse them here.
        //
        rc.AllocateRichQualityValues(length);
        DNALength pos;

        for (pos = 0; pos < length; pos++) {
            if (insertionQV.Empty() == false) {
                rc.insertionQV[length - pos - 1] = insertionQV[pos];
            }
            if (substitutionQV.Empty() == false) {
                rc.substitutionQV[length - pos - 1]           = substitutionQV[pos];
            }
            if (deletionQV.Empty() == false) {
                rc.deletionQV[length - pos - 1]               = deletionQV[pos];
            }

            if (mergeQV.Empty() == false) {
                rc.mergeQV[length - pos - 1]               = mergeQV[pos];
            }


            if (substitutionTag != NULL) {
                rc.substitutionTag[length - pos - 1] = ReverseComplementNuc[substitutionTag[pos]];
            }
            if (deletionTag != NULL) {
                rc.deletionTag[length - pos - 1]     = ReverseComplementNuc[deletionTag[pos]];
            }
        }
    }
    deletionQVPrior = rc.deletionQVPrior;
    insertionQVPrior = rc.insertionQVPrior;
    substitutionQVPrior = rc.substitutionQVPrior;
    preBaseDeletionQVPrior = rc.preBaseDeletionQVPrior;
}

void FASTQSequence::Free() {
    FASTASequence::Free();
    if (deleteOnExit == true) {
        qual.Free();
        deletionQV.Free();
        preBaseDeletionQV.Free();
        insertionQV.Free();
        substitutionQV.Free();
        mergeQV.Free();
        if (deletionTag != NULL) {
            delete[] deletionTag;
            deletionTag = NULL;
        }
        if (substitutionTag != NULL) {
            delete[] substitutionTag;
            substitutionTag = NULL;
        }
    }
}

void FASTQSequence::LowerCaseMask(int qThreshold) {
    int i;
    if (qual.Empty() == true) return;

    for (i = 0; i < length; i++ ){
        if (qual[i] < qThreshold) {
            seq[i] = tolower(seq[i]);
        }
    }
}

float FASTQSequence::GetAverageQuality() {
    DNALength p;
    float totalQ;
    if (qual.Empty() == true) { return 0.0; }
    assert(qual.Empty() == false);
    assert(length > 0);
    for (p = 0, totalQ = 0.0; p < length; p++) {
        totalQ += qual[p];
    }
    return totalQ / length;
}
