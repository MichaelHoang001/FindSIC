/* 
 * File:   misc.cpp
 * Author: Andrew Scott
 * 
 * Created on June 9, 2011, 2:12 PM
 */

#include "misc.h"

bool fileexists(const char *filename)
{
    std::ifstream ifile(filename);
    return ifile;
}

unsigned numchars(unsigned n)
{
    unsigned c = 0;
    if (n==0) {
        c = 1;
    }
    else {
        while (n) {
            n /= 10;
            c++;
        }
    }
    return c;
}


