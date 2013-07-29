#ifndef BASE_TUPLE_H_
#define BASE_TUPLE_H_


class BaseTuple {
public:
    ULong  tuple;

    ULong  HashPowerOfFour(int nBases, TupleMetrics &tm);

    int operator<(const BaseTuple &rhs) const; 

    int operator==(const BaseTuple &rhs) const; 

    int operator!= (const BaseTuple &rhs) const;

    BaseTuple ShiftLeft(TupleMetrics &tm, ULong shift=1L); 

    BaseTuple ShiftRight(ULong shift=1L); 

    BaseTuple Append(ULong val, TupleMetrics &tm, ULong nBits); 

    long ToLongIndex(); 
};


#endif
