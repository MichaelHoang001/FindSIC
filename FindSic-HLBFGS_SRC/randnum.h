/* 
 * File:   randnum.h
 * Author: Andrew Scott
 *
 * Created on June 9, 2011, 8:19 PM
 */

#ifndef RANDNUM_H
#define	RANDNUM_H

#include <cstdlib>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>


// Base random number generator class
template<class T> class randBase {
public:
    randBase(unsigned d = 1);
    virtual ~randBase();
    unsigned getDim() const;
    virtual void setDim(unsigned d)=0;
    virtual void getNextD(std::vector<T> &x)=0;
    virtual T getNext()=0;
protected:
    unsigned dim;
};
template<class T> randBase<T>::randBase(unsigned d) : dim(d) {}
template<class T> randBase<T>::~randBase() {}
template<class T> unsigned randBase<T>::getDim() const {
    return dim;
}

// Standard random number generator
template<class T> class randDefault: public randBase<T> {
public:
    randDefault(unsigned d = 1);
    void setDim(unsigned d);
    void getNextD(std::vector<T> &x);
    T getNext();
}; 
template<class T> randDefault<T>::randDefault(unsigned d) : randBase<T>(d) {}
template<class T> void randDefault<T>::setDim(unsigned d) {
    randBase<T>::dim = d;
}
template<class T> void randDefault<T>::getNextD(std::vector<T> &x) {
    typename std::vector<T>::iterator px;
    for (px = x.begin(); px != x.end(); ++px) {
        *px = static_cast<T>(std::rand())/static_cast<T>(RAND_MAX);
    }
}
template<class T> T randDefault<T>::getNext() {
    return static_cast<T>(std::rand())/static_cast<T>(RAND_MAX);
}

// Sobol quasi-random number generator
template<class T> class randSobol: public randBase<T> {
public:
    randSobol(unsigned d = 1, std::string dir_file_ = "");
    virtual ~randSobol();
    void setDim(unsigned d);
    void getNextD(std::vector<T> &x);
    T getNext();
private:
    std::string dir_file;
    std::ifstream infile;
    unsigned L;
    T LF;
    unsigned N;
    unsigned C;
    unsigned **V;
    unsigned *X, *pX;
    typename std::vector<T>::iterator px;
    void initialise();
    void clear();
}; 
template<class T> void randSobol<T>::initialise() {
  std::ifstream infile(dir_file.c_str());
  if (!infile) {
    std::cout << "Input file containing direction numbers cannot be found!" << std::endl;
    exit(1);
  }
  infile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  L = 32;
  LF = std::pow(static_cast<T>(2),static_cast<T>(L));
  // Compute direction numbers V[1] to V[L], scaled by LF=pow(2,L)
  V = new unsigned * [L+1];
  for (unsigned i = 0; i <= L; i++) {
      V[i] = new unsigned [randBase<T>::dim];
  }
  for (unsigned i = 1; i <= L; i++) { 
      V[i][0] = 1 << (L-i); // all m's = 1
  }
  for (unsigned j = 1; j < randBase<T>::dim; j++) {
    // Read in parameters from file 
    unsigned d, s, a;
    infile >> d >> s >> a;
    unsigned *m = new unsigned [s+1];
    for (unsigned i = 1; i<=s;i++) infile >> m[i];
    if (L <= s) {
      for (unsigned i=1;i<=L;i++) V[i][j] = m[i] << (L-i); 
    }
    else {
      for (unsigned i=1;i<=s;i++) V[i][j] = m[i] << (L-i); 
      for (unsigned i=s+1;i<=L;i++) {
	V[i][j] = V[i-s][j] ^ (V[i-s][j] >> s); 
	for (unsigned k=1;k<=s-1;k++) 
	  V[i][j] ^= (((a >> (s-1-k)) & 1) * V[i-k][j]); 
      }
    }
    delete [] m;    
  }
  X = new unsigned [randBase<T>::dim];
  pX = &X[0];
  *pX = 0;
  for (unsigned i = 1; i < randBase<T>::dim; i++) {
      ++pX;
      *pX = 0;
  }  
  N = 1;
  C = 1;  
}
template<class T> void randSobol<T>::clear() {
    for (unsigned i = 0; i <= L; i++) {
        delete [] V[i];
    }
    delete [] V;
    delete [] X;
}
template<class T> randSobol<T>::randSobol(unsigned d, std::string dir_file_) : randBase<T>(d), dir_file(dir_file_) {
    initialise();
}
template<class T> randSobol<T>::~randSobol() {
    clear();
}
template<class T> void randSobol<T>::setDim(unsigned d) {
    randBase<T>::dim = d;
    clear();
    initialise();
}
template<class T> void randSobol<T>::getNextD(std::vector<T> &x) {    
    pX = &X[0];
    px = x.begin();
    *pX ^= V[C][0];
    *px = static_cast<T>(*pX) / LF;
    for (unsigned j = 1; j < randBase<T>::dim; ++j) {
        ++pX;++px;
        *pX ^= V[C][j];
        *px = static_cast<T>(*pX) / LF;
    }
    C = 1;
    unsigned v = N;
    while (v & 1) {
        v >>= 1;
        C++;
    }
    ++N;
}
template<class T> T randSobol<T>::getNext() {
    return static_cast<T>(0);
}


// Mersenne Twister random number generator
template<class T> class randMT: public randBase<T> {
public:
    randMT(unsigned d = 1);
    virtual ~randMT();
    void setDim(unsigned d);
    void getNextD(std::vector<T> &x);
    T getNext();
private:
    static const int N                    = 624;
    static const int M                    = 397;
    // constant vector a
    static const unsigned long MATRIX_A   = 0x9908b0dfUL;
    // most significant w-r bits
    static const unsigned long UPPER_MASK = 0x80000000UL;
    // least significant r bits
    static const unsigned long LOWER_MASK = 0x7fffffffUL;

    unsigned long* mt_;                  // the state vector
    int mti_;                            // mti == N+1 means mt not initialized

    unsigned long* init_key_;            // Storage for the seed vector
    int key_length_;                     // Seed vector length
    unsigned long s_;                    // Seed integer
    bool seeded_by_array_;               // Seeded by an array
    bool seeded_by_int_;                 // Seeded by an integer
    
    void init_genrand(unsigned long s);
    void init_by_array(unsigned long* init_key, int key_length);

    unsigned long genrand_int32(); // Generates a random number on [0,0xffffffff]
    T genrand_real1(); // Generates a random real number on [0,1]
    T genrand_real2(); // Generates a random real number on [0,1)
    T genrand_real3(); // Generates a random real number on (0,1)
    T genrand_res53(); // Generates a random real number on [0,1) with 53-bit precision

    void initialise();
    void clear();
    
}; 
template<class T> void randMT<T>::initialise() {
    unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 };
    unsigned long length = 4;
    init_by_array(init, length);  
}
template<class T> void randMT<T>::clear() {
    delete[] mt_;
    delete[] init_key_;
}
template<class T> randMT<T>::randMT(unsigned d) : randBase<T>(d), mt_(new unsigned long[N]), mti_(N+1), init_key_(NULL), key_length_(0), s_(0), seeded_by_array_(false), seeded_by_int_(false) {
    initialise();
}
template<class T> randMT<T>::~randMT() {
    clear();
}
template<class T> void randMT<T>::setDim(unsigned d) {
    randBase<T>::dim = d;
}
/**
 * Initializes the Mersenne Twister with a seed.
 *
 * \param s seed
 */
template<class T> void randMT<T>::init_genrand(unsigned long s)
{
    mt_[0]= s & 0xffffffffUL;
    for (mti_=1; mti_<N; mti_++) {
        mt_[mti_] = 
	    (1812433253UL * (mt_[mti_-1] ^ (mt_[mti_-1] >> 30)) + mti_); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt_[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt_[mti_] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
    // Store the seed
    s_ = s;
    seeded_by_array_ = false;
    seeded_by_int_ = true;
}

/**
 * Seed the Mersenne Twister using an array.
 *
 * \param init_key an array for initializing keys
 * \param key_length the length of \a init_key
 */
template<class T> void randMT<T>::init_by_array(unsigned long* init_key, int key_length)
{
    // Store the key array
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt_[i] = (mt_[i] ^ ((mt_[i-1] ^ (mt_[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt_[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt_[0] = mt_[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt_[i] = (mt_[i] ^ ((mt_[i-1] ^ (mt_[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt_[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt_[0] = mt_[N-1]; i=1; }
    }

    mt_[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 

    // Store the seed
    if (init_key_ != NULL) {
        delete[] init_key_;
    }
    init_key_ = new unsigned long[key_length];
    for (int k = 0; k < key_length; k++) {
        init_key_[k] = init_key[k];
    }
    key_length_ = key_length;
    seeded_by_int_ = false;
    seeded_by_array_ = true;
}
template<class T> unsigned long randMT<T>::genrand_int32() {
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti_ >= N) { /* generate N words at one time */
        int kk;

        if (mti_ == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt_[kk]&UPPER_MASK)|(mt_[kk+1]&LOWER_MASK);
            mt_[kk] = mt_[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt_[kk]&UPPER_MASK)|(mt_[kk+1]&LOWER_MASK);
            mt_[kk] = mt_[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt_[N-1]&UPPER_MASK)|(mt_[0]&LOWER_MASK);
        mt_[N-1] = mt_[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti_ = 0;
    }
  
    y = mt_[mti_++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}
template<class T> T randMT<T>::genrand_real1() {
    return static_cast<T>(genrand_int32())*(static_cast<T>(1.0)/static_cast<T>(4294967295.0)); // divided by 2^32-1 
}
template<class T> T randMT<T>::genrand_real2() {
    return static_cast<T>(genrand_int32())*(static_cast<T>(1.0)/static_cast<T>(4294967296.0)); // divided by 2^32 
}
template<class T> T randMT<T>::genrand_real3() {
    return (static_cast<T>(genrand_int32()) + 0.5)*(static_cast<T>(1.0)/static_cast<T>(4294967296.0)); // divided by 2^32 
}
template<class T> T randMT<T>::genrand_res53() { 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return (static_cast<T>(a)*static_cast<T>(67108864.0)+static_cast<T>(b))*(static_cast<T>(1.0)/static_cast<T>(9007199254740992.0)); 
} 
template<class T> void randMT<T>::getNextD(std::vector<T> &x) {
    typename std::vector<T>::iterator px;
    for (px = x.begin(); px != x.end(); ++px) {
        *px = genrand_real1();
    }
}
template<class T> T randMT<T>::getNext() {
    return genrand_real1();
}




#endif	/* RANDNUM_H */

