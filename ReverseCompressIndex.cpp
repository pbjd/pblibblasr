#include <iostream>
#include <fstream>
#include "ReverseCompressIndex.hpp"

ReverseCompressIndex::ReverseCompressIndex() {
    index = NULL;
    indexLength = binSize = maxRun = 0;
}

void ReverseCompressIndex::Write(std::ofstream &out) {
    out.write((char*) &indexLength, sizeof(int));
    out.write((char*) &binSize, sizeof(int));
    out.write((char*) &maxRun, sizeof(int));
    out.write((char*) index, sizeof(int) * indexLength);
}

void ReverseCompressIndex::Read(std::ifstream &in) {
    in.read((char*) &indexLength, sizeof(int));
    in.read((char*) &binSize, sizeof(int));
    in.read((char*) &maxRun, sizeof(int));
    index = new int[indexLength];
    in.read((char*) index, sizeof(int) *indexLength);
} 

void ReverseCompressIndex::ShallowCopy(ReverseCompressIndex &rhs) {
    index = rhs.index;
    indexLength = rhs.indexLength;
    binSize = rhs.binSize;
    maxRun  = rhs.maxRun;
}
