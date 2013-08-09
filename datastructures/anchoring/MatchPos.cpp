#include "MatchPos.hpp"

MatchPos::MatchPos(DNALength pt, DNALength pq, DNALength pl, int pm) :
    t(pt),
    q(pt),
    l(pl),
    m(pm),
    w(0)
 { }

MatchPos::MatchPos(const MatchPos &rhs) {
    (*this) = rhs;
}

MatchWeight MatchPos::Size() {
    return l;
}

MatchPos::MatchPos() {
    t = q = -1;
    l = 0;
    w = 0;
    m = 0;
}

MatchPos& MatchPos::operator=(const MatchPos &rhs) {
    t = rhs.t; q = rhs.q; w = rhs.w;
    l = rhs.l;
    m = rhs.m;
    return *this;
}

DNALength MatchPos::GetLength() const {
    return l;
}

int MatchPos::GetMultiplicity() const {
    return m;
}

MatchWeight MatchPos::GetWeight() const {
    if (m > 0) {
        return (1.0*l)/m;
    } else {
        return 0;
    }
}

DNALength MatchPos::GetX() const {
    return q;
}

DNALength MatchPos::GetY() const {
    return t;
}

UInt MatchPos::GetT() {
    return t;
}

UInt MatchPos::GetQ() {
    return (UInt) q;
}

UInt MatchPos::GetW() {
    return w;
}

std::ostream& operator<<(std::ostream & out, MatchPos &p) {
    out << p.q << "\t" << p.t <<"\t"<< p.l << "\t"<< p.m;
    return out;
}


ChainedMatchPos::ChainedMatchPos(DNALength pt, DNALength pq, DNALength pl, int pm) : MatchPos(pt, pq, pl, pm) {
    score = 0; 
    chainPrev = NULL;
}

ChainedMatchPos::ChainedMatchPos() : MatchPos() {
    score = 0;
    chainPrev = NULL;
}

ChainedMatchPos::ChainedMatchPos(const ChainedMatchPos &rhs) {
    (*this) = rhs;
}

int ChainedMatchPos::GetScore() {
    return score;
}

int ChainedMatchPos::SetScore(int _score) {
    return (score = _score);
}

ChainedMatchPos* ChainedMatchPos::SetChainPrev(ChainedMatchPos *_chainPrev) {
    return (chainPrev = _chainPrev);
}

ChainedMatchPos* ChainedMatchPos::GetChainPrev() {
    return chainPrev;
}

ChainedMatchPos &ChainedMatchPos::operator=(const ChainedMatchPos &rhs) {
    ((MatchPos&)(*this)) = ((MatchPos&)rhs);
    return *this;
}

std::ostream& operator<<(std::ostream & out, ChainedMatchPos &p) {
    out << p.q << "\t" << p.t <<"\t"<< p.l << "\t"<< p.m;
    return out;
}
