#ifndef UTILS_STRING_UTILS_H_
#define UTILS_STRING_UTILS_H_

#include <cassert>
#include <string>
#include <vector>
#include <sstream>
#include "MD5Utils.hpp"
#include "Types.h"

using namespace std;

int ExactPatternMatch(string orig, string pattern); 

void MakeMD5(const char *data, unsigned int dataLength, string &md5Str, int nChars = 0); 

void MakeMD5(string &data, string &md5Str, int nChars=0); 

int IsWhitespace(char c); 

int IsSpace(char c); 

int ToWords(string &orig, vector<string> &words); 

int Tokenize(string orig, string pattern, vector<string> &tokens); 

void ParseSeparatedList(const string &csl, vector<string> &values, char delim=','); 

int AssignUntilFirstSpace(char *orig, int origLength, string &result); 

template<typename T_Value>
void ParseSeparatedList(const string &csl, vector<T_Value> &values, char delim=',') {
    stringstream cslStrm(csl);
    T_Value val;
    string valString;
    string next;
    do {
        getline(cslStrm, valString, delim);
        if (cslStrm and valString.size() > 0) {
            stringstream valStrm(valString);
            if (! (valStrm  >> val) ) {
                cout << "Error, value " << valString << " is malformatted." << endl;
            }
            else {
                values.push_back(val);
            }
        }
    }
    while (cslStrm);
}

#endif
