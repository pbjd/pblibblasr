#include "FragmentSort.hpp"

template<typename T_Fragment>
int 
LexicographicFragmentSort<T_Fragment>::operator()
    (const T_Fragment &a, const T_Fragment &b) const {
    return a.LessThan(b);
}
