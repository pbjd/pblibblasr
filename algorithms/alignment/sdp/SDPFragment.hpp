#ifndef _BLASR_SDP_FRAGMENT_HPP_
#define _BLASR_SDP_FRAGMENT_HPP_

class Fragment {
 public:
	unsigned int x;
	unsigned int y;
	unsigned int weight;
	unsigned int length;
	int index;
	int chainPrev;
	int cost;
	unsigned int chainLength;

	Fragment(unsigned int px, unsigned int py, int pweight=0); 

	Fragment();

	unsigned int GetX() const; 

	unsigned int GetY() const; 

    int LessThan(const Fragment &f) const {
        if (x < f.x)
            return 1;
        else if (x == f.x) 
            return y < f.y;
        else 
            return 0;
    }

    int operator<(const Fragment &f) const {
        // 
        // Sort fragments by diagonal:
        //
        int diag, fDiag;
        diag = (y - x);
        fDiag = f.y - f.x;
        if (diag < fDiag)
            return 1;
        else if (diag == fDiag)
            return (x < f.x);
        else
            return 0;
    }

    Fragment& operator=(const Fragment &rhs) {
        x           = rhs.x;
        y           = rhs.y;
        index       = rhs.index;
        cost        = rhs.cost;
        weight      = rhs.weight;
        chainLength = rhs.chainLength;
        chainPrev   = rhs.chainPrev;
        return *this;
    }
		
    int operator==(const Fragment &f) const {
        return (x == f.x and y == f.y);
    }

	int operator>(const Fragment &f) const; 

	int GetLength(); 

	void SetLength(int _length); 
};


#endif // _BLASR_SDP_FRAGMENT_HPP_
