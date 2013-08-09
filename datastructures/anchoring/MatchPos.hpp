#ifndef _BLASR_MATCH_POS_HPP_
#define _BLASR_MATCH_POS_HPP_

#include <vector>
#include <algorithm>
#include <ostream>
#include "Types.h"
#include "DNASequence.hpp"

class MatchPos {
public:
    DNALength t, q;
    MatchWeight w;
    DNALength l;
    int m; // multiplicity

    MatchPos(DNALength pt, DNALength pq, DNALength pl, int pm = 0); 

    MatchPos (const MatchPos &rhs); 

    MatchWeight Size(); 

    MatchPos(); 

    MatchPos& operator=(const MatchPos &rhs); 

    DNALength GetLength() const; 
    
    int GetMultiplicity() const; 
    
    MatchWeight GetWeight() const; 

    DNALength GetX() const; 
    
    DNALength GetY() const; 
    
    UInt GetT(); 

    UInt GetQ(); 

    UInt GetW(); 

    friend std::ostream& operator<<(std::ostream & out, MatchPos &p); 
};


class ChainedMatchPos : public MatchPos {
private:
    int score;
    ChainedMatchPos *chainPrev;

public:
    ChainedMatchPos(DNALength pt, DNALength pq, DNALength pl, int pm);

    ChainedMatchPos();

    ChainedMatchPos(const ChainedMatchPos &rhs); 

    int GetScore();

    int SetScore(int _score); 

    ChainedMatchPos* SetChainPrev(ChainedMatchPos *_chainPrev); 

    ChainedMatchPos* GetChainPrev();

    ChainedMatchPos &operator=(const ChainedMatchPos &rhs); 

    friend std::ostream& operator<<(std::ostream & out, ChainedMatchPos &p); 

};

template<typename T_MatchPos>
class CompareMatchPos {
    public:
        int operator()(const T_MatchPos &lhs, const T_MatchPos &rhs) const {
            if (lhs.t < rhs.t) 
                return 1;
            else if (lhs.t > rhs.t)
                return 0;
            else { 
                return lhs.q < rhs.q;
            }
        }
};

typedef std::vector<MatchPos> MatchPosList;

template<typename T_MatchPos>
class CompareMatchPosByWeight {
    public:
        int operator()(const T_MatchPos &a, const T_MatchPos &b) const {
            return a.l < b.l;
        }
};

template<typename T_MatchPos>
class CompareMatchPosIndexByWeight {
    public:
        std::vector<T_MatchPos> *list;
        int operator()(const int i, const int j) const {
            return ((*list)[i].w > (*list)[j].w);
        }
};


template<typename T_MatchPos>
class CompareMatchPosIndexByTextPos {
    public:
        std::vector<T_MatchPos> *list;
        int operator()(const int i, const int j) const {
            return (*list)[i].t < (*list)[j].t;
        }
};


template<typename T_MatchPos>
void SortMatchPosList(std::vector<T_MatchPos> &mpl) {
    std::sort(mpl.begin(), mpl.end(), CompareMatchPos<T_MatchPos>());
}

template<typename T_MatchPos>
void SortMatchPosListByWeight(std::vector<T_MatchPos> &mpl) {
    std::sort(mpl.begin(), mpl.end(), CompareMatchPosByWeight<T_MatchPos>());
}

template<typename T_MatchPos>
void SortMatchPosIndexListByWeight(std::vector<T_MatchPos> &mpl, std::vector<int> &indices) {
    CompareMatchPosIndexByWeight<T_MatchPos> cmp;
    cmp.list = &mpl;
    std::sort(indices.begin(), indices.end(), cmp);
}

template<typename T_MatchPos>
void SortMatchPosIndexListByTextPos(std::vector<T_MatchPos> &mpl, std::vector<int> &indices) {
    CompareMatchPosIndexByTextPos<T_MatchPos> cmp;
    cmp.list = &mpl;
    std::sort(indices.begin(), indices.end(), cmp);
}

#endif  // _BLASR_MATCH_POS_HPP_
