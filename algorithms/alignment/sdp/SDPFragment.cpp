#include "SDPFragment.hpp"

Fragment::Fragment(unsigned int px, unsigned int py, int pweight) {
    x = px;
    y = py;
    weight = pweight;
    length = index = 0;
    chainPrev = cost = chainLength = 0;
}
//
// Provide default constructor that will
// give bad results if members are not properly initialized
// later on.
//
Fragment::Fragment() {
    x = -1;
    y = -1;
    weight = length = index = 0;
    chainPrev = cost = chainLength = 0;
}

int Fragment::operator>(const Fragment &f) const {
    return (!(*this < f) &&  !(*this == f));
}
int Fragment::GetLength() {
    return length;
}
void Fragment::SetLength(int _length) {
    length = _length;
}
