/* 
 * File:   misc.h
 * Author: Andrew Scott
 *
 * Created on June 9, 2011, 2:12 PM
 */

#ifndef MISC_H
#define	MISC_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cmath>

template<class T> void printx(std::vector<T> &x) {
    
    typename std::vector<T>::iterator px; 
        
    for (px = x.begin(); px != x.end(); ++px) {
        std::cout << *px << std::endl;
    }
    
}

template<class T> void printu(std::vector<T> &x) {
    
    const unsigned d = x.size()/2;
    typename std::vector<T>::iterator px; 
    
    px = x.begin();    
    for (unsigned i = d; i > 0; --i) {
        std::cout << *px << " + " << *(px+d) << "i" << std::endl;
        ++px;
    }
    
}

template<class T> T normx(std::vector<T> &x) {
    typename std::vector<T>::iterator px;
    T r = 0.0;
    for (px = x.begin(); px != x.end(); ++px) {
        r += (*px)*(*px);
    }
    return std::sqrt(r);
}

template<class T> T normxmy(std::vector<T> &x, std::vector<T> &y) {
    typename std::vector<T>::iterator px, py;
    T r = 0.0;
    for (px = x.begin(), py = y.begin(); px != x.end(); ++px, ++py) {
        r += (*px-*py)*(*px-*py);
    }
    return std::sqrt(r);
}

template<class T> void normalx(std::vector<T> &x) {
    
    typename std::vector<T>::iterator px;
    T r = 0.0;
        
    for (px = x.begin(); px != x.end(); ++px) {
        r += (*px)*(*px);
    }
    r = std::sqrt(r);
    for (px = x.begin(); px != x.end(); ++px) {
        *px /= r;
    }  
}

template<class T> void projectZx(std::vector<T> &x) {
    
    std::complex<T> Ipi;
    Ipi = std::complex<T>(0,3.14159265358979323846264338327950288419716939937510);
    
    unsigned d = x.size()/2;
    
    std::complex<T> * Zc = new std::complex<T> [d*d];
    std::complex<T> *PZc_jk;
    
    int j, k;
    
    for (k=0, PZc_jk=&Zc[0]; k<d; ++k) {
        for (j=0; j<d; j++, ++PZc_jk) {
            *PZc_jk = (std::exp(-Ipi*(T((2*k+(d+1)*j)*j)/T(d)+T(d-1)/T(12)))+std::exp(Ipi*(T((2*j+(d+1)*k)*k)/T(d)+T(d-1)/T(12))))/std::sqrt(T(d));
            if (j==k) {
                *PZc_jk += 1; 
            }
            *PZc_jk /= 3;
        }
    }
    
    std::complex<T> * u = new std::complex<T> [d];
    std::complex<T> *Pu;
    
    typename std::vector<T>::iterator Ix = x.begin(), Iux, Iuy;
    
    for (k=d, Pu=&u[0], PZc_jk=&Zc[0]; k>0; --k, ++Pu) {
      *Pu = 0;
      for (j=d, Iux=Ix, Iuy=Ix+d; j>0; --j, ++Iux, ++Iuy, ++PZc_jk) {
          *Pu += (*PZc_jk)*(std::complex<T>(*Iux, *Iuy));
      }
    }
 
    for (k=d, Pu=&u[0], Iux=Ix, Iuy=Ix+d; k>0; --k, ++Pu, ++Iux, ++Iuy) {
        *Iux = (*Pu).real();
        *Iuy = (*Pu).imag();
    } 
    
    delete [] u;
    delete [] Zc;
  
}


template<class T> void rnumstox(std::vector<T> &x, std::vector<T> &r) {
    
    const T twicepi = 2*3.14159265358979323846264338327950288419716939937510;
    const unsigned d = x.size()/2;
    unsigned i,j;
    typename std::vector<T>::iterator pr, px, pe;
    px = x.begin();
    pe = x.begin()+d;
    pr = r.begin();
         
    *pe = *pr;
    for (i = 2; i < d; ++i) {
        ++pe;
        ++pr;
        *pe = *pr;
    }
    
    pe = x.begin()+d;
    *px = *pe;
    ++px;
    *px = 1-*pe;
    for (i = 2; i < d; ++i) {
        ++px;
        ++pe;
        *px = 1-pow(*pe,1.0/i);
    }
   
    for (j = d-1; j > 1; --j) {
        px = x.begin();
        for (i = 0; i < j; ++i) {
            *px *= pow(*pe,1.0/j);
            ++px;
        }
        --pe;
    }
   
    px = x.begin();
    pe = x.begin() + d;
    for (i = 0; i < d; ++i) {
        ++pr;
        *px = sqrt(*px);
        *pe = (*px)*sin(twicepi*(*pr));
        *px *= cos(twicepi*(*pr));
        ++px;
        ++pe;
    }
}

bool fileexists(const char *filename);

unsigned numchars(unsigned n);



#endif	/* MISC_H */

/*
template<class T> void funcSIC<T>::getFuncGradCC(T* Px, T* Pfunc, T* Pgrad) {

  T *Pux_j, *Puy_j, *Pgradx, *Pgrady; 
  std::complex<T> *Pu_j, *Pu_jpd, *Pu_jpdpd, *Pu_jpj, *Pu_l, *Pu_jpl, *Pu_jpjpl, *Pu_k, *Pu_jpk, *Pu_kpl, *Pu_jpkpl, *Ps_jpdj, *Ps_jpdk, *Ps_kpdj, *Pu_dmj, *Pu_dpdmjmj, *Pu_dmjpl, *Pu_dpdmjmjpl, *Pu_dmk, *Pu_dpdmjmk, *Pu_dmjpk, *Pu_dmkpj, *Pu_dmkpl, *Pu_dpdmjmkpl, *Pu_dmjpkpl, *Pu_dmkpjpl, *Pds__jpdj, *Pdsc__jpdj, *Pds_l_jpdj, *Pdsc_l_jpdj, *Pds__jpdk, *Pds__kpdj, *Pdsc__jpdk, *Pdsc__kpdj, *Pds_l_jpdk, *Pds_l_kpdj, *Pdsc_l_jpdk, *Pdsc_l_kpdj, *Pdssdz_j, *Ps_k, *Pds_j_k, *Pdsc_j_k;
  unsigned twod = 2*d, dp1 = d+1, dm1 = d-1, dd = d*d;
  
  su = 0;
  for (j=0, Pu_j = &u[0], Pux_j=Px, Puy_j=Px+d; j<d; ++j, ++Pu_j, ++Pux_j, ++Puy_j) {
   *Pu_j = std::complex<T>(*Pux_j, *Puy_j);
   su += (*Pux_j)*(*Pux_j) + (*Puy_j)*(*Puy_j);
  }
  su = std::sqrt(su);
  for (j=0, Pu_j=&u[0], Pu_jpd=&u[0]+d, Pu_jpdpd=&u[0]+twod; j<d; ++j, ++Pu_j, ++Pu_jpd, ++Pu_jpdpd) {
   *Pu_j /= std::complex<T>(su);
   *Pu_jpd = *Pu_j;
   *Pu_jpdpd = *Pu_j;
  }  

  *Pfunc = 0;
  for (j=0, Pu_j=&u[0], Pu_jpj=&u[0], Ps_jpdj=&s[0]; j<d; ++j, ++Pu_j, ++Pu_jpj, ++Pu_jpj, Ps_jpdj+=dp1) {
   Sl = 0;
   for (l=0, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_jpjpl=Pu_jpj; l<d; ++l, ++Pu_l, ++Pu_jpl, ++Pu_jpjpl) {
    Sl += std::conj((*Pu_jpl)*(*Pu_jpl))*(*Pu_l)*(*Pu_jpjpl);
   }
   *Ps_jpdj = Sl;
   *Pfunc += Sl.real()*Sl.real() + Sl.imag()*Sl.imag();
   for (k=j+1, Pu_k=Pu_j+1, Pu_jpk=Pu_jpj+1, Ps_jpdk=Ps_jpdj+d, Ps_kpdj=Ps_jpdj+1; k<d; ++k, ++Pu_k, ++Pu_jpk, Ps_jpdk+=d, ++Ps_kpdj) {
    Sl = 0;
    for (l=0, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_kpl=Pu_k, Pu_jpkpl=Pu_jpk; l<d; ++l, ++Pu_l, ++Pu_jpl, ++Pu_kpl, ++Pu_jpkpl) {
     Sl += std::conj((*Pu_jpl)*(*Pu_kpl))*(*Pu_l)*(*Pu_jpkpl);
    }
    *Ps_jpdk = Sl;
    *Ps_kpdj = Sl;
    *Pfunc += 2.0*(Sl.real()*Sl.real()+Sl.imag()*Sl.imag());
   }
  }
  *Pfunc *= d+1;
  *Pfunc -= 2;

  for (j=0, Pu_j=&u[0], Pu_jpj=&u[0], Pu_dmj=&u[0]+d, Pu_dpdmjmj=&u[0]+twod, Pds__jpdj=&ds2[0], Pdsc__jpdj=&dsc2[0]; j<d; ++j, ++Pu_j, ++Pu_jpj, ++Pu_jpj, --Pu_dmj, --Pu_dpdmjmj, --Pu_dpdmjmj, Pds__jpdj+=dp1, Pdsc__jpdj+=dp1) {
   for (l=0, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_jpjpl=Pu_jpj, Pu_dmjpl=Pu_dmj, Pu_dpdmjmjpl=Pu_dpdmjmj, Pds_l_jpdj=Pds__jpdj, Pdsc_l_jpdj=Pdsc__jpdj; l<d; ++l, ++Pu_l, ++Pu_jpl, ++Pu_jpjpl, ++Pu_dmjpl, ++Pu_dpdmjmjpl, Pds_l_jpdj+=dd, Pdsc_l_jpdj+=dd) {
    *Pds_l_jpdj = std::conj((*Pu_jpl)*(*Pu_jpl))*(*Pu_jpjpl) + std::conj((*Pu_dmjpl)*(*Pu_dmjpl))*(*Pu_dpdmjmjpl);
    *Pdsc_l_jpdj = std::complex<T>(2.0)*std::conj((*Pu_dmjpl)*(*Pu_jpl))*(*Pu_l);
   }
   for (k=j+1, Pu_k=Pu_j+1, Pu_jpk=Pu_jpj+1, Pu_dmk=Pu_dmj-1, Pu_dpdmjmk=Pu_dpdmjmj-1, Pu_dmjpk=&u[0]+dp1, Pu_dmkpj=&u[0]+dm1, Pds__jpdk=Pds__jpdj+d, Pds__kpdj=Pds__jpdj+1, Pdsc__jpdk=Pdsc__jpdj+d, Pdsc__kpdj=Pdsc__jpdj+1; k<d; ++k, ++Pu_k, ++Pu_jpk, --Pu_dmk, --Pu_dpdmjmk, ++Pu_dmjpk, --Pu_dmkpj, Pds__jpdk+=d, ++Pds__kpdj, Pdsc__jpdk+=d, ++Pdsc__kpdj) {
    for (l=0, Pu_jpl=Pu_j, Pu_kpl=Pu_k, Pu_jpkpl=Pu_jpk, Pu_dmjpl=Pu_dmj, Pu_dmkpl=Pu_dmk, Pu_dpdmjmkpl=Pu_dpdmjmk, Pu_dmjpkpl=Pu_dmjpk, Pu_dmkpjpl=Pu_dmkpj, Pds_l_jpdk=Pds__jpdk, Pds_l_kpdj=Pds__kpdj, Pdsc_l_jpdk=Pdsc__jpdk, Pdsc_l_kpdj=Pdsc__kpdj; l<d; ++l, ++Pu_jpl, ++Pu_kpl, ++Pu_jpkpl, ++Pu_dmjpl, ++Pu_dmkpl, ++Pu_dpdmjmkpl, ++Pu_dmjpkpl, ++Pu_dmkpjpl, Pds_l_jpdk+=dd, Pds_l_kpdj+=dd, Pdsc_l_jpdk+=dd, Pdsc_l_kpdj+=dd) {
     *Pds_l_jpdk = std::conj((*Pu_jpl)*(*Pu_kpl))*(*Pu_jpkpl) + std::conj((*Pu_dmjpl)*(*Pu_dmkpl))*(*Pu_dpdmjmkpl);
     *Pdsc_l_jpdk = std::conj((*Pu_dmjpl)*(*Pu_kpl))*(*Pu_dmjpkpl) + std::conj((*Pu_dmkpl)*(*Pu_jpl))*(*Pu_dmkpjpl);
     *Pds_l_kpdj = *Pds_l_jpdk;
     *Pdsc_l_kpdj = *Pdsc_l_jpdk;
    }
   }
  }

  for (j=0, Pdssdz_j=&dssdz[0], Pds_j_k=&ds2[0], Pdsc_j_k=&dsc2[0]; j<d; ++j, ++Pdssdz_j) {
   *Pdssdz_j = 0;
   for (k=0, Ps_k=&s[0]; k<dd; ++k, ++Ps_k, ++Pds_j_k, ++Pdsc_j_k) {
    *Pdssdz_j += (*Pdsc_j_k)*(*Ps_k) + (*Pds_j_k)*std::conj(*Ps_k);
   }
  }

  a = 0;
  for (j=0, Pdssdz_j=&dssdz[0], Pux_j=Px, Puy_j=Px+d; j<d; ++j, ++Pdssdz_j, ++Pux_j, ++Puy_j) {
   a += (*Pdssdz_j).real()*(*Pux_j) - (*Pdssdz_j).imag()*(*Puy_j);
  }

  a /= su*su;
  b = 2.0*dp1/su;

  for (j=0, Pdssdz_j=&dssdz[0], Pux_j=Px, Puy_j=Px+d, Pgradx=Pgrad, Pgrady=Pgrad+d; j<d; ++j, ++Pdssdz_j, ++Pux_j, ++Puy_j, ++Pgradx, ++Pgrady) {
   *Pgradx =  b*( (*Pdssdz_j).real() - a*(*Pux_j) );
   *Pgrady = -b*( (*Pdssdz_j).imag() + a*(*Puy_j) );
  }

}
*/