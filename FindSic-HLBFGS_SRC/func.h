/* 
 * File:   func.h
 * Author: Andrew Scott
 *
 * Created on June 10, 2011, 1:44 PM
 */

#ifndef FUNC_H
#define	FUNC_H

//#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <complex>


template<class T> class funcBase {
public:
    funcBase(unsigned n = 1);
    virtual ~funcBase();
    unsigned getNumVars() const;
    virtual void setNumVars(unsigned n)=0;
    virtual void getFuncGrad(std::vector<T> &x, T &func, std::vector<T> &grad)=0;
    virtual void getFunc(std::vector<T> &x, T &func)=0;
protected:
    unsigned numvars;
};
template<class T> funcBase<T>::funcBase(unsigned n) : numvars(n) {}
template<class T> funcBase<T>::~funcBase() {}
template<class T> unsigned funcBase<T>::getNumVars() const {
    return numvars;
}

template<class T> class funcSIC: public funcBase<T> {
public:
    funcSIC(unsigned n = 4);
    virtual ~funcSIC();
    void setNumVars(unsigned n);
    void getFuncGrad(std::vector<T> &x, T &func, std::vector<T> &grad);
    void getFunc(std::vector<T> &x, T &func);
    void getFuncGradC(T* xptr, T* funcptr, T* gradptr);
private:
    unsigned d;
    std::complex<T> *u, *s, **ds, **dsc, *dssdz; 
    unsigned j, k, l;
    typename std::vector<T>::iterator Px, Py;
    typename std::vector<T>::iterator PDE, Pr1, Pr2;
    T *CPx, *CPy, *CPr1, *CPr2, *CPDE;
    T *PE, su, a, b;
    std::complex<T> *Pu, *Pc1, *Pc2, *Pc3, *Pc4, Sl, *Pl, *Pj, *Pk, *Pjpk, *Pjpl, *Pkpl, *Pjpkpl, *P2j, *Pdmj, *P2dm2j, *Pdmjpl, *P2dm2jpl, *Pdmk, *P2dmjmk, *Pdmjpk, *Pdmkpj, *Pdmkpl, *P2dmjmkpl, *Pdmjpkpl, *Pdmkpjpl;
    void initialise();
    void clear();
};

template<class T> void funcSIC<T>::initialise() {
    u = new std::complex<T>[3*d];
    s = new std::complex<T>[d*d];
    ds = new std::complex<T>*[d];
    for (j=0; j<d; j++) {
        ds[j] = new std::complex<T>[d*d];
    }
    dsc = new std::complex<T>*[d];
    for (j=0; j<d; j++) {
        dsc[j] = new std::complex<T>[d*d];
    }
    dssdz = new std::complex<T>[d];
}

template<class T> void funcSIC<T>::clear() {    
    delete [] u;
    delete [] s;
    for (j=0; j<d; j++) {
        delete [] ds[j];
    }
    delete [] ds;
    for (j=0; j<d; j++) {
        delete [] dsc[j];
    }
    delete [] dsc;
    delete [] dssdz;
    u=0;
    s=0;
    ds=0;
    dsc=0;
    dssdz=0;
}

template<class T> funcSIC<T>::funcSIC(unsigned n) : funcBase<T>(n) {
    d = n/2;
    initialise();
}

template<class T> funcSIC<T>::~funcSIC() { 
    clear(); 
}

template<class T> void funcSIC<T>::setNumVars(unsigned n) {
    if (n != funcBase<T>::numvars) {
        funcBase<T>::numvars = n;
        clear();
        d = n/2;
        initialise();
    }
}

template<class T> void funcSIC<T>::getFuncGrad(std::vector<T> &x, T &func, std::vector<T> &grad) {
    
  Px=x.begin();
  PE=&func;
  PDE=grad.begin();

  Pu=&u[0];
  Py=Px+d;
  su=0;
  for (j=0; j<d; j++) {
   *Pu=std::complex<T>(*Px,*Py);
   su+=*Px*(*Px)+*Py*(*Py);
   Pu++;
   Px++;
   Py++;
  }
  su=std::sqrt(su);
  Pc1=Pu-d;
  for (j=0; j<d; j++) {
   *Pc1/=std::complex<T>(su);
   *Pu=*Pc1;
   Pu++;
   Pc1++;
  }
  Pc1=Pu-d;
  for (j=0; j<d; j++) {
   *Pu=*Pc1;
   Pu++;
   Pc1++;
  }
  Pu-=3*d;

  *PE=0;
  Pj=Pu;
  P2j=Pu;
  for (j=0; j<d; j++) {
   Sl=0;
   Pl=Pu;
   Pjpl=Pj;
   Pjpkpl=P2j;
   for (l=0; l<d; l++) {
    Sl+=std::conj(*Pjpl*(*Pjpl))*(*Pl)*(*Pjpkpl);
    ++Pl;
    ++Pjpl;
    ++Pjpkpl;
   }
   s[j+d*j]=Sl;
   *PE+=Sl.real()*Sl.real()+Sl.imag()*Sl.imag();
   Pk=Pj+1;
   Pjpk=P2j+1;
   for (k=(j+1); k<d; k++) {
    Sl=0;
    Pl=Pu;
    Pjpl=Pj;
    Pkpl=Pk;
    Pjpkpl=Pjpk;
    for (l=0; l<d; l++) {
     Sl+=std::conj(*Pjpl*(*Pkpl))*(*Pl)*(*Pjpkpl);
     ++Pl;
     ++Pjpl;
     ++Pkpl;
     ++Pjpkpl;
    }
    s[j+d*k]=Sl;
    s[k+d*j]=Sl;
    *PE+=2*(Sl.real()*Sl.real()+Sl.imag()*Sl.imag());
   ++Pk;
   ++Pjpk;
   }
  ++Pj;
  ++P2j;++P2j;
  }
  *PE*=d+1;
  *PE-=2;

  Pj=Pu;
  P2j=Pu;
  Pdmj=Pu+d;
  P2dm2j=Pu+2*d;
  for (j=0; j<d; j++) {
   Pl=Pu;
   Pjpl=Pj;
   Pjpkpl=P2j;
   Pdmjpl=Pdmj;
   P2dm2jpl=P2dm2j;
   for (l=0; l<d; l++) {
    ds[l][j+j*d]=std::conj(*Pjpl*(*Pjpl))*(*Pjpkpl)+std::conj(*Pdmjpl*(*Pdmjpl))*(*P2dm2jpl);
    dsc[l][j+j*d]=std::complex<T>(2)*std::conj((*Pdmjpl)*(*Pjpl))*(*Pl);
    ++Pl;
    ++Pjpl;
    ++Pjpkpl;
    ++Pdmjpl;
    ++P2dm2jpl;
   }
   Pk=Pj+1;
   Pjpk=P2j+1;
   Pdmk=Pdmj-1;
   P2dmjmk=P2dm2j-1;
   Pdmjpk=Pu+d+1;
   Pdmkpj=Pu+d-1;
   for (k=(j+1); k<d; k++) {
    Pjpl=Pj;
    Pkpl=Pk;
    Pjpkpl=Pjpk;
    Pdmjpl=Pdmj;
    Pdmkpl=Pdmk;
    P2dmjmkpl=P2dmjmk;
    Pdmjpkpl=Pdmjpk;
    Pdmkpjpl=Pdmkpj;
    for (l=0; l<d; l++) {
     ds[l][j+k*d]=std::conj(*Pjpl*(*Pkpl))*(*Pjpkpl)+std::conj(*Pdmjpl*(*Pdmkpl))*(*P2dmjmkpl);
     dsc[l][j+k*d]=std::conj(*Pdmjpl*(*Pkpl))*(*Pdmjpkpl)+std::conj(*Pdmkpl*(*Pjpl))*(*Pdmkpjpl);
     ds[l][k+j*d]=ds[l][j+d*k];
     dsc[l][k+j*d]=dsc[l][j+d*k];
     ++Pjpl;
     ++Pkpl;
     ++Pjpkpl;
     ++Pdmjpl;
     ++Pdmkpl;
     ++P2dmjmkpl;
     ++Pdmjpkpl;
     ++Pdmkpjpl;
    }
    ++Pk;
    ++Pjpk;
    --Pdmk;
    --P2dmjmk;
    ++Pdmjpk;
    --Pdmkpj;
   }
   ++Pj;
   ++P2j;++P2j;
   --Pdmj;
   --P2dm2j;--P2dm2j;
  }

  Pc1=&dssdz[0];
  for (j=0; j<d; j++) {
   Pc3=&ds[j][0];
   Pc4=&dsc[j][0];
   *Pc1=0;
   Pc2=&s[0];
   for (k=0; k<(d*d); k++) {
    *Pc1+=*Pc4*(*Pc2)+*Pc3*std::conj(*Pc2);
    ++Pc2;
    ++Pc3;
    ++Pc4;
   }
   ++Pc1;
  }

  a=0;
  Pc1=&dssdz[0];
  Py=Px;
  Px-=d;
  for (j=0; j<d; j++) {
   a+=(*Pc1).real()*(*Px);
   a-=(*Pc1).imag()*(*Py);
   ++Px;
   ++Py;
   ++Pc1;
  }
  Py=Px;
  Px-=d;

  a/=su*su;
  b=2*(d+1)/su;

  Pc1=&dssdz[0];
  Pr1=PDE;
  Pr2=PDE+d;
  for (j=0; j<d; j++) {
   *Pr1=b*((*Pc1).real()-a*(*Px));
   *Pr2=b*(-(*Pc1).imag()-a*(*Py));
   ++Px;
   ++Py;
   ++Pc1;
   ++Pr1;
   ++Pr2;
  }

}

template<class T> void funcSIC<T>::getFunc(std::vector<T> &x, T &func) {
    
  Px=x.begin();
  PE=&func;

  Pu=&u[0];
  Py=Px+d;
  su=0;
  for (j=0; j<d; j++) {
   *Pu=std::complex<T>(*Px,*Py);
   su+=*Px*(*Px)+*Py*(*Py);
   Pu++;
   Px++;
   Py++;
  }
  su=std::sqrt(su);
  Pc1=Pu-d;
  for (j=0; j<d; j++) {
   *Pc1/=std::complex<T>(su);
   *Pu=*Pc1;
   Pu++;
   Pc1++;
  }
  Pc1=Pu-d;
  for (j=0; j<d; j++) {
   *Pu=*Pc1;
   Pu++;
   Pc1++;
  }
  Pu-=3*d;

  *PE=0;
  Pj=Pu;
  P2j=Pu;
  for (j=0; j<d; j++) {
   Sl=0;
   Pl=Pu;
   Pjpl=Pj;
   Pjpkpl=P2j;
   for (l=0; l<d; l++) {
    Sl+=std::conj(*Pjpl*(*Pjpl))*(*Pl)*(*Pjpkpl);
    ++Pl;
    ++Pjpl;
    ++Pjpkpl;
   }
   s[j+d*j]=Sl;
   *PE+=Sl.real()*Sl.real()+Sl.imag()*Sl.imag();
   Pk=Pj+1;
   Pjpk=P2j+1;
   for (k=(j+1); k<d; k++) {
    Sl=0;
    Pl=Pu;
    Pjpl=Pj;
    Pkpl=Pk;
    Pjpkpl=Pjpk;
    for (l=0; l<d; l++) {
     Sl+=std::conj(*Pjpl*(*Pkpl))*(*Pl)*(*Pjpkpl);
     ++Pl;
     ++Pjpl;
     ++Pkpl;
     ++Pjpkpl;
    }
    s[j+d*k]=Sl;
    s[k+d*j]=Sl;
    *PE+=2*(Sl.real()*Sl.real()+Sl.imag()*Sl.imag());
   ++Pk;
   ++Pjpk;
   }
  ++Pj;
  ++P2j;++P2j;
  }
  *PE*=d+1;
  *PE-=2;

}

template<class T> void funcSIC<T>::getFuncGradC(T* xptr, T* funcptr, T* gradptr) {

  CPx=xptr;
  PE=funcptr;
  CPDE=gradptr;
  
  Pu=&u[0];
  CPy=CPx+d;
  su=0;
  for (j=0; j<d; j++) {
   *Pu=std::complex<T>(*CPx,*CPy);
   su+=*CPx*(*CPx)+*CPy*(*CPy);
   Pu++;
   CPx++;
   CPy++;
  }
  su=std::sqrt(su);
  Pc1=Pu-d;
  for (j=0; j<d; j++) {
   *Pc1/=std::complex<T>(su);
   *Pu=*Pc1;
   Pu++;
   Pc1++;
  }
  Pc1=Pu-d;
  for (j=0; j<d; j++) {
   *Pu=*Pc1;
   Pu++;
   Pc1++;
  }
  Pu-=3*d;

  *PE=0;
  Pj=Pu;
  P2j=Pu;
  for (j=0; j<d; j++) {
   Sl=0;
   Pl=Pu;
   Pjpl=Pj;
   Pjpkpl=P2j;
   for (l=0; l<d; l++) {
    Sl+=std::conj(*Pjpl*(*Pjpl))*(*Pl)*(*Pjpkpl);
    ++Pl;
    ++Pjpl;
    ++Pjpkpl;
   }
   s[j+d*j]=Sl;
   *PE+=Sl.real()*Sl.real()+Sl.imag()*Sl.imag();
   Pk=Pj+1;
   Pjpk=P2j+1;
   for (k=(j+1); k<d; k++) {
    Sl=0;
    Pl=Pu;
    Pjpl=Pj;
    Pkpl=Pk;
    Pjpkpl=Pjpk;
    for (l=0; l<d; l++) {
     Sl+=std::conj(*Pjpl*(*Pkpl))*(*Pl)*(*Pjpkpl);
     ++Pl;
     ++Pjpl;
     ++Pkpl;
     ++Pjpkpl;
    }
    s[j+d*k]=Sl;
    s[k+d*j]=Sl;
    *PE+=2*(Sl.real()*Sl.real()+Sl.imag()*Sl.imag());
   ++Pk;
   ++Pjpk;
   }
  ++Pj;
  ++P2j;++P2j;
  }
  *PE*=d+1;
  *PE-=2;

  Pj=Pu;
  P2j=Pu;
  Pdmj=Pu+d;
  P2dm2j=Pu+2*d;
  for (j=0; j<d; j++) {
   Pl=Pu;
   Pjpl=Pj;
   Pjpkpl=P2j;
   Pdmjpl=Pdmj;
   P2dm2jpl=P2dm2j;
   for (l=0; l<d; l++) {
    ds[l][j+j*d]=std::conj(*Pjpl*(*Pjpl))*(*Pjpkpl)+std::conj(*Pdmjpl*(*Pdmjpl))*(*P2dm2jpl);
    dsc[l][j+j*d]=std::complex<T>(2)*std::conj((*Pdmjpl)*(*Pjpl))*(*Pl);
    ++Pl;
    ++Pjpl;
    ++Pjpkpl;
    ++Pdmjpl;
    ++P2dm2jpl;
   }
   Pk=Pj+1;
   Pjpk=P2j+1;
   Pdmk=Pdmj-1;
   P2dmjmk=P2dm2j-1;
   Pdmjpk=Pu+d+1;
   Pdmkpj=Pu+d-1;
   for (k=(j+1); k<d; k++) {
    Pjpl=Pj;
    Pkpl=Pk;
    Pjpkpl=Pjpk;
    Pdmjpl=Pdmj;
    Pdmkpl=Pdmk;
    P2dmjmkpl=P2dmjmk;
    Pdmjpkpl=Pdmjpk;
    Pdmkpjpl=Pdmkpj;
    for (l=0; l<d; l++) {
     ds[l][j+k*d]=std::conj(*Pjpl*(*Pkpl))*(*Pjpkpl)+std::conj(*Pdmjpl*(*Pdmkpl))*(*P2dmjmkpl);
     dsc[l][j+k*d]=std::conj(*Pdmjpl*(*Pkpl))*(*Pdmjpkpl)+std::conj(*Pdmkpl*(*Pjpl))*(*Pdmkpjpl);
     ds[l][k+j*d]=ds[l][j+d*k];
     dsc[l][k+j*d]=dsc[l][j+d*k];
     ++Pjpl;
     ++Pkpl;
     ++Pjpkpl;
     ++Pdmjpl;
     ++Pdmkpl;
     ++P2dmjmkpl;
     ++Pdmjpkpl;
     ++Pdmkpjpl;
    }
    ++Pk;
    ++Pjpk;
    --Pdmk;
    --P2dmjmk;
    ++Pdmjpk;
    --Pdmkpj;
   }
   ++Pj;
   ++P2j;++P2j;
   --Pdmj;
   --P2dm2j;--P2dm2j;
  }

  Pc1=&dssdz[0];
  for (j=0; j<d; j++) {
   Pc3=&ds[j][0];
   Pc4=&dsc[j][0];
   *Pc1=0;
   Pc2=&s[0];
   for (k=0; k<(d*d); k++) {
    *Pc1+=*Pc4*(*Pc2)+*Pc3*std::conj(*Pc2);
    ++Pc2;
    ++Pc3;
    ++Pc4;
   }
   ++Pc1;
  }

  a=0;
  Pc1=&dssdz[0];
  CPy=CPx;
  CPx-=d;
  for (j=0; j<d; j++) {
   a+=(*Pc1).real()*(*CPx);
   a-=(*Pc1).imag()*(*CPy);
   ++CPx;
   ++CPy;
   ++Pc1;
  }
  CPy=CPx;
  CPx-=d;

  a/=su*su;
  b=2*(d+1)/su;

  Pc1=&dssdz[0];
  CPr1=CPDE;
  CPr2=CPDE+d;
  for (j=0; j<d; j++) {
   *CPr1=b*((*Pc1).real()-a*(*CPx));
   *CPr2=b*(-(*Pc1).imag()-a*(*CPy));
   ++CPx;
   ++CPy;
   ++Pc1;
   ++CPr1;
   ++CPr2;
  }

}

// Faster version

template<class T> class funcSICf: public funcBase<T> {
public:
    funcSICf(unsigned n = 4);
    virtual ~funcSICf();
    void setNumVars(unsigned n);
    void getFuncGrad(std::vector<T> &x, T &func, std::vector<T> &grad);
    void getFunc(std::vector<T> &x, T &func);
    void getFuncGradC(T* Px, T* Pfunc, T* Pgrad);
    void getFuncGradC2(T* Px, T* Pfunc, T* Pgrad);
private:
    unsigned d;
    std::complex<T> *u, *s, *ds, *dsc, *dssdz; 
    unsigned j, k, l;
    T su, a, b;
    std::complex<T> Sl;
    unsigned twod, dp1, dm1, dd;
    T *Pux_j, *Puy_j, *Pgradx, *Pgrady; 
    std::complex<T> *Pu_j, *Pu_jpd, *Pu_jpdpd, *Pu_jpj, *Pu_l, *Pu_jpl, *Pu_jpjpl, *Pu_k, *Pu_jpk, *Pu_kpl, *Pu_jpkpl, *Ps_jpdj, *Ps_jpdk, *Ps_kpdj, *Pu_dmj, *Pu_dpdmjmj, *Pu_dmjpl, *Pu_dpdmjmjpl, *Pu_dmk, *Pu_dpdmjmk, *Pu_dmjpk, *Pu_dmkpj, *Pu_dmkpl, *Pu_dpdmjmkpl, *Pu_dmjpkpl, *Pu_dmkpjpl, *Pds__jpdj, *Pdsc__jpdj, *Pds_l_jpdj, *Pdsc_l_jpdj, *Pds__jpdk, *Pds__kpdj, *Pdsc__jpdk, *Pdsc__kpdj, *Pds_l_jpdk, *Pds_l_kpdj, *Pdsc_l_jpdk, *Pdsc_l_kpdj, *Pdssdz_j, *Ps_k, *Pds_j_k, *Pdsc_j_k;
    void initialise();
    void clear();
};

template<class T> void funcSICf<T>::initialise() {
    u = new std::complex<T>[3*d];
    s = new std::complex<T>[d*d];
    ds = new std::complex<T>[d*d*d];
    dsc = new std::complex<T>[d*d*d];
    dssdz = new std::complex<T>[d];
}

template<class T> void funcSICf<T>::clear() {    
    delete [] u;
    delete [] s;
    delete [] ds;
    delete [] dsc;
    delete [] dssdz;
    u=0;
    s=0;
    ds=0;
    dsc=0;
    dssdz=0;
}

template<class T> funcSICf<T>::funcSICf(unsigned n) : funcBase<T>(n) {
    d = n/2;
    twod = 2*d;
    dp1 = d+1;
    dm1 = d-1;
    dd = d*d;
    initialise();
}

template<class T> funcSICf<T>::~funcSICf() { 
    clear(); 
}

template<class T> void funcSICf<T>::setNumVars(unsigned n) {
    if (n != funcBase<T>::numvars) {
        funcBase<T>::numvars = n;
        clear();
        d = n/2;
        twod = 2*d;
        dp1 = d+1;
        dm1 = d-1;
        dd = d*d;
        initialise();
    }
}

template<class T> void funcSICf<T>::getFuncGrad(std::vector<T> &x, T &func, std::vector<T> &grad) {}

template<class T> void funcSICf<T>::getFunc(std::vector<T> &x, T &func) {
      
  typename std::vector<T>::iterator Ix=x.begin(), Iux_j, Iuy_j;  
  
  su = 0;
  for (j=0, Pu_j = &u[0], Iux_j=Ix, Iuy_j=Ix+d; j<d; ++j, ++Pu_j, ++Iux_j, ++Iuy_j) {
   *Pu_j = std::complex<T>(*Iux_j, *Iuy_j);
   su += (*Iux_j)*(*Iux_j) + (*Iuy_j)*(*Iuy_j);
  }
  su = std::sqrt(su);
  for (j=0, Pu_j=&u[0], Pu_jpd=&u[0]+d, Pu_jpdpd=&u[0]+twod; j<d; ++j, ++Pu_j, ++Pu_jpd, ++Pu_jpdpd) {
   *Pu_j /= std::complex<T>(su);
   *Pu_jpd = *Pu_j;
   *Pu_jpdpd = *Pu_j;
  }  

  func = 0;
  for (j=0, Pu_j=&u[0], Pu_jpj=&u[0]; j<d; ++j, ++Pu_j, ++Pu_jpj, ++Pu_jpj) {
   Sl = 0;
   for (l=0, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_jpjpl=Pu_jpj; l<d; ++l, ++Pu_l, ++Pu_jpl, ++Pu_jpjpl) {
    Sl += std::conj((*Pu_jpl)*(*Pu_jpl))*(*Pu_l)*(*Pu_jpjpl);
   }
   func += Sl.real()*Sl.real() + Sl.imag()*Sl.imag();
   for (k=j+1, Pu_k=Pu_j+1, Pu_jpk=Pu_jpj+1; k<d; ++k, ++Pu_k, ++Pu_jpk) {
    Sl = 0;
    for (l=0, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_kpl=Pu_k, Pu_jpkpl=Pu_jpk; l<d; ++l, ++Pu_l, ++Pu_jpl, ++Pu_kpl, ++Pu_jpkpl) {
     Sl += std::conj((*Pu_jpl)*(*Pu_kpl))*(*Pu_l)*(*Pu_jpkpl);
    }
    func += 2.0*(Sl.real()*Sl.real()+Sl.imag()*Sl.imag());
   }
  }
  func *= d+1;
  func -= 2;

}

template<class T> void funcSICf<T>::getFuncGradC(T* Px, T* Pfunc, T* Pgrad) {
  
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

  for (j=0, Pu_j=&u[0], Pu_jpj=&u[0], Pu_dmj=&u[0]+d, Pu_dpdmjmj=&u[0]+twod, Pds__jpdj=&ds[0], Pdsc__jpdj=&dsc[0]; j<d; ++j, ++Pu_j, ++Pu_jpj, ++Pu_jpj, --Pu_dmj, --Pu_dpdmjmj, --Pu_dpdmjmj, Pds__jpdj+=dp1, Pdsc__jpdj+=dp1) {
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

  for (j=0, Pdssdz_j=&dssdz[0], Pds_j_k=&ds[0], Pdsc_j_k=&dsc[0]; j<d; ++j, ++Pdssdz_j) {
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

// Faster version without "new"

template<class T> class funcSICff: public funcBase<T> {
public:
    funcSICff(unsigned n = 4);
    virtual ~funcSICff();
    void setNumVars(unsigned n);
    void getFuncGrad(std::vector<T> &x, T &func, std::vector<T> &grad);
    void getFunc(std::vector<T> &x, T &func);
    void getFuncGradC(T* Px, T* Pfunc, T* Pgrad);
    void getFuncGradC2(T* Px, T* Pfunc, T* Pgrad);
private:
    int d;
    std::complex<T> u[480], s[25600], dssdz[160], ds[4096000], dsc[4096000]; 
    int j, k, l;
    T su, a, b;
    std::complex<T> Sl;
    int twod, dp1, dm1, dd;
    T *Pux_j, *Puy_j, *Pgradx, *Pgrady; 
    std::complex<T> *Pu_j, *Pu_jpd, *Pu_jpdpd, *Pu_jpj, *Pu_l, *Pu_jpl, *Pu_jpjpl, *Pu_k, *Pu_jpk, *Pu_kpl, *Pu_jpkpl, *Ps_jpdj, *Ps_jpdk, *Ps_kpdj, *Pu_dmj, *Pu_dpdmjmj, *Pu_dmjpl, *Pu_dpdmjmjpl, *Pu_dmk, *Pu_dpdmjmk, *Pu_dmjpk, *Pu_dmkpj, *Pu_dmkpl, *Pu_dpdmjmkpl, *Pu_dmjpkpl, *Pu_dmkpjpl, *Pds__jpdj, *Pdsc__jpdj, *Pds_l_jpdj, *Pdsc_l_jpdj, *Pds__jpdk, *Pds__kpdj, *Pdsc__jpdk, *Pdsc__kpdj, *Pds_l_jpdk, *Pds_l_kpdj, *Pdsc_l_jpdk, *Pdsc_l_kpdj, *Pdssdz_j, *Ps_k, *Pds_j_k, *Pdsc_j_k;
    void initialise();
    void clear();
};

template<class T> void funcSICff<T>::initialise() {}

template<class T> void funcSICff<T>::clear() {}

template<class T> funcSICff<T>::funcSICff(unsigned n) : funcBase<T>(n) {
    d = n/2;
    twod = 2*d;
    dp1 = d+1;
    dm1 = d-1;
    dd = d*d;
    initialise();
}

template<class T> funcSICff<T>::~funcSICff() { 
    clear(); 
}

template<class T> void funcSICff<T>::setNumVars(unsigned n) {
    if (n != funcBase<T>::numvars) {
        funcBase<T>::numvars = n;
        clear();
        d = n/2;
        twod = 2*d;
        dp1 = d+1;
        dm1 = d-1;
        dd = d*d;
        initialise();
    }
}

template<class T> void funcSICff<T>::getFuncGrad(std::vector<T> &x, T &func, std::vector<T> &grad) {}

template<class T> void funcSICff<T>::getFunc(std::vector<T> &x, T &func) {
  
  typename std::vector<T>::iterator Ix=x.begin(), Iux_j, Iuy_j;  
  
  su = 0;
  for (j=0, Pu_j = &u[0], Iux_j=Ix, Iuy_j=Ix+d; j<d; ++j, ++Pu_j, ++Iux_j, ++Iuy_j) {
   *Pu_j = std::complex<T>(*Iux_j, *Iuy_j);
   su += (*Iux_j)*(*Iux_j) + (*Iuy_j)*(*Iuy_j);
  }
  su = std::sqrt(su);
  for (j=0, Pu_j=&u[0], Pu_jpd=&u[0]+d, Pu_jpdpd=&u[0]+twod; j<d; ++j, ++Pu_j, ++Pu_jpd, ++Pu_jpdpd) {
   *Pu_j /= std::complex<T>(su);
   *Pu_jpd = *Pu_j;
   *Pu_jpdpd = *Pu_j;
  }  

  func = 0;
  for (j=0, Pu_j=&u[0], Pu_jpj=&u[0]; j<d; ++j, ++Pu_j, ++Pu_jpj, ++Pu_jpj) {
   Sl = 0;
   for (l=0, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_jpjpl=Pu_jpj; l<d; ++l, ++Pu_l, ++Pu_jpl, ++Pu_jpjpl) {
    Sl += std::conj((*Pu_jpl)*(*Pu_jpl))*(*Pu_l)*(*Pu_jpjpl);
   }
   func += Sl.real()*Sl.real() + Sl.imag()*Sl.imag();
   for (k=j+1, Pu_k=Pu_j+1, Pu_jpk=Pu_jpj+1; k<d; ++k, ++Pu_k, ++Pu_jpk) {
    Sl = 0;
    for (l=0, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_kpl=Pu_k, Pu_jpkpl=Pu_jpk; l<d; ++l, ++Pu_l, ++Pu_jpl, ++Pu_kpl, ++Pu_jpkpl) {
     Sl += std::conj((*Pu_jpl)*(*Pu_kpl))*(*Pu_l)*(*Pu_jpkpl);
    }
    func += 2.0*(Sl.real()*Sl.real()+Sl.imag()*Sl.imag());
   }
  }
  func *= d+1;
  func -= 2;

}

template<class T> void funcSICff<T>::getFuncGradC(T* Px, T* Pfunc, T* Pgrad) {
  
  su = 0;
  for (j=d, Pu_j = &u[0], Pux_j=Px, Puy_j=Px+d; j>0; --j, ++Pu_j, ++Pux_j, ++Puy_j) {
   *Pu_j = std::complex<T>(*Pux_j, *Puy_j);
   su += (*Pux_j)*(*Pux_j) + (*Puy_j)*(*Puy_j);
  }
  su = std::sqrt(su);
  for (j=d, Pu_j=&u[0], Pu_jpd=&u[0]+d, Pu_jpdpd=&u[0]+twod; j>0; --j, ++Pu_j, ++Pu_jpd, ++Pu_jpdpd) {
   *Pu_j /= std::complex<T>(su);
   *Pu_jpd = *Pu_j;
   *Pu_jpdpd = *Pu_j;
  }  

  *Pfunc = 0;
  for (j=d, Pu_j=&u[0], Pu_jpj=&u[0], Ps_jpdj=&s[0]; j>0; --j, ++Pu_j, ++Pu_jpj, ++Pu_jpj, Ps_jpdj+=dp1) {
   Sl = 0;
   for (l=d, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_jpjpl=Pu_jpj; l>0; --l, ++Pu_l, ++Pu_jpl, ++Pu_jpjpl) {
    Sl += std::conj((*Pu_jpl)*(*Pu_jpl))*(*Pu_l)*(*Pu_jpjpl);
   }
   *Ps_jpdj = Sl;
   *Pfunc += Sl.real()*Sl.real() + Sl.imag()*Sl.imag();
   for (k=j-1, Pu_k=Pu_j+1, Pu_jpk=Pu_jpj+1, Ps_jpdk=Ps_jpdj+d, Ps_kpdj=Ps_jpdj+1; k>0; --k, ++Pu_k, ++Pu_jpk, Ps_jpdk+=d, ++Ps_kpdj) {
    Sl = 0;
    for (l=d, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_kpl=Pu_k, Pu_jpkpl=Pu_jpk; l>0; --l, ++Pu_l, ++Pu_jpl, ++Pu_kpl, ++Pu_jpkpl) {
     Sl += std::conj((*Pu_jpl)*(*Pu_kpl))*(*Pu_l)*(*Pu_jpkpl);
    }
    *Ps_jpdk = Sl;
    *Ps_kpdj = Sl;
    *Pfunc += 2.0*(Sl.real()*Sl.real()+Sl.imag()*Sl.imag());
   }
  }
  *Pfunc *= d+1;
  *Pfunc -= 2;

  for (j=d, Pu_j=&u[0], Pu_jpj=&u[0], Pu_dmj=&u[0]+d, Pu_dpdmjmj=&u[0]+twod, Pds__jpdj=&ds[0], Pdsc__jpdj=&dsc[0]; j>0; --j, ++Pu_j, ++Pu_jpj, ++Pu_jpj, --Pu_dmj, --Pu_dpdmjmj, --Pu_dpdmjmj, Pds__jpdj+=dp1, Pdsc__jpdj+=dp1) {
   for (l=d, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_jpjpl=Pu_jpj, Pu_dmjpl=Pu_dmj, Pu_dpdmjmjpl=Pu_dpdmjmj, Pds_l_jpdj=Pds__jpdj, Pdsc_l_jpdj=Pdsc__jpdj; l>0; --l, ++Pu_l, ++Pu_jpl, ++Pu_jpjpl, ++Pu_dmjpl, ++Pu_dpdmjmjpl, Pds_l_jpdj+=dd, Pdsc_l_jpdj+=dd) {
    *Pds_l_jpdj = std::conj((*Pu_jpl)*(*Pu_jpl))*(*Pu_jpjpl) + std::conj((*Pu_dmjpl)*(*Pu_dmjpl))*(*Pu_dpdmjmjpl);
    *Pdsc_l_jpdj = std::complex<T>(2.0)*std::conj((*Pu_dmjpl)*(*Pu_jpl))*(*Pu_l);
   }
   for (k=j-1, Pu_k=Pu_j+1, Pu_jpk=Pu_jpj+1, Pu_dmk=Pu_dmj-1, Pu_dpdmjmk=Pu_dpdmjmj-1, Pu_dmjpk=&u[0]+dp1, Pu_dmkpj=&u[0]+dm1, Pds__jpdk=Pds__jpdj+d, Pds__kpdj=Pds__jpdj+1, Pdsc__jpdk=Pdsc__jpdj+d, Pdsc__kpdj=Pdsc__jpdj+1; k>0; --k, ++Pu_k, ++Pu_jpk, --Pu_dmk, --Pu_dpdmjmk, ++Pu_dmjpk, --Pu_dmkpj, Pds__jpdk+=d, ++Pds__kpdj, Pdsc__jpdk+=d, ++Pdsc__kpdj) {
    for (l=d, Pu_jpl=Pu_j, Pu_kpl=Pu_k, Pu_jpkpl=Pu_jpk, Pu_dmjpl=Pu_dmj, Pu_dmkpl=Pu_dmk, Pu_dpdmjmkpl=Pu_dpdmjmk, Pu_dmjpkpl=Pu_dmjpk, Pu_dmkpjpl=Pu_dmkpj, Pds_l_jpdk=Pds__jpdk, Pds_l_kpdj=Pds__kpdj, Pdsc_l_jpdk=Pdsc__jpdk, Pdsc_l_kpdj=Pdsc__kpdj; l>0; --l, ++Pu_jpl, ++Pu_kpl, ++Pu_jpkpl, ++Pu_dmjpl, ++Pu_dmkpl, ++Pu_dpdmjmkpl, ++Pu_dmjpkpl, ++Pu_dmkpjpl, Pds_l_jpdk+=dd, Pds_l_kpdj+=dd, Pdsc_l_jpdk+=dd, Pdsc_l_kpdj+=dd) {
     *Pds_l_jpdk = std::conj((*Pu_jpl)*(*Pu_kpl))*(*Pu_jpkpl) + std::conj((*Pu_dmjpl)*(*Pu_dmkpl))*(*Pu_dpdmjmkpl);
     *Pdsc_l_jpdk = std::conj((*Pu_dmjpl)*(*Pu_kpl))*(*Pu_dmjpkpl) + std::conj((*Pu_dmkpl)*(*Pu_jpl))*(*Pu_dmkpjpl);
     *Pds_l_kpdj = *Pds_l_jpdk;
     *Pdsc_l_kpdj = *Pdsc_l_jpdk;
    }
   }
  }

  for (j=d, Pdssdz_j=&dssdz[0], Pds_j_k=&ds[0], Pdsc_j_k=&dsc[0]; j>0; --j, ++Pdssdz_j) {
   *Pdssdz_j = 0;
   for (k=dd, Ps_k=&s[0]; k>0; --k, ++Ps_k, ++Pds_j_k, ++Pdsc_j_k) {
    *Pdssdz_j += (*Pdsc_j_k)*(*Ps_k) + (*Pds_j_k)*std::conj(*Ps_k);
   }
  }

  a = 0;
  for (j=d, Pdssdz_j=&dssdz[0], Pux_j=Px, Puy_j=Px+d; j>0; --j, ++Pdssdz_j, ++Pux_j, ++Puy_j) {
   a += (*Pdssdz_j).real()*(*Pux_j) - (*Pdssdz_j).imag()*(*Puy_j);
  }

  a /= su*su;
  b = 2.0*dp1/su;

  for (j=d, Pdssdz_j=&dssdz[0], Pux_j=Px, Puy_j=Px+d, Pgradx=Pgrad, Pgrady=Pgrad+d; j>0; --j, ++Pdssdz_j, ++Pux_j, ++Puy_j, ++Pgradx, ++Pgrady) {
   *Pgradx =  b*( (*Pdssdz_j).real() - a*(*Pux_j) );
   *Pgrady = -b*( (*Pdssdz_j).imag() + a*(*Puy_j) );
  }

}


// Now impose Zauner symmetry:

// Faster version without "new"

template<class T> class funcZSICff: public funcBase<T> {
public:
    funcZSICff(unsigned n = 4);
    virtual ~funcZSICff();
    void setNumVars(unsigned n);
    void getFuncGrad(std::vector<T> &x, T &func, std::vector<T> &grad);
    void getFunc(std::vector<T> &x, T &func);
    void getFuncGradC(T* Px, T* Pfunc, T* Pgrad);
    void getFuncGradC2(T* Px, T* Pfunc, T* Pgrad);
private:
    int d;
    std::complex<T> u[480], s[25600], dssduc[160], dssdzc[160], ds[4096000], dsc[4096000], Zc[25600]; 
    int j, k, l;
    T su, a, b;
    std::complex<T> Sl;
    int twod, dp1, dm1, dd;
    T *Pux_j, *Puy_j, *Pgradx, *Pgrady; 
    std::complex<T> *Pu_j, *Pu_jpd, *Pu_jpdpd, *Pu_jpj, *Pu_l, *Pu_jpl, *Pu_jpjpl, *Pu_k, *Pu_jpk, *Pu_kpl, *Pu_jpkpl, *Ps_jpdj, *Ps_jpdk, *Ps_kpdj, *Pu_dmj, *Pu_dpdmjmj, *Pu_dmjpl, *Pu_dpdmjmjpl, *Pu_dmk, *Pu_dpdmjmk, *Pu_dmjpk, *Pu_dmkpj, *Pu_dmkpl, *Pu_dpdmjmkpl, *Pu_dmjpkpl, *Pu_dmkpjpl, *Pds__jpdj, *Pdsc__jpdj, *Pds_l_jpdj, *Pdsc_l_jpdj, *Pds__jpdk, *Pds__kpdj, *Pdsc__jpdk, *Pdsc__kpdj, *Pds_l_jpdk, *Pds_l_kpdj, *Pdsc_l_jpdk, *Pdsc_l_kpdj, *Pdssduc_j, *Pdssduc_k, *Pdssdzc_j, *Ps_k, *Pds_j_k, *Pdsc_j_k, *PZc_jk;
    void initialise();
    void clear();
};

template<class T> void funcZSICff<T>::initialise() {
     
    std::complex<T> Ipi;
    Ipi = std::complex<T>(0,3.14159265358979323846264338327950288419716939937510);
    
    for (k=0, PZc_jk=&Zc[0]; k<d; ++k) {
        for (j=0; j<d; j++, ++PZc_jk) {
            *PZc_jk = (std::exp(-Ipi*(T((2*k+(d+1)*j)*j)/T(d)+T(d-1)/T(12)))+std::exp(Ipi*(T((2*j+(d+1)*k)*k)/T(d)+T(d-1)/T(12))))/std::sqrt(T(d));
            if (j==k) {
                *PZc_jk += 1; 
            }
            *PZc_jk /= 3;
        }
    }
    
}

template<class T> void funcZSICff<T>::clear() {}

template<class T> funcZSICff<T>::funcZSICff(unsigned n) : funcBase<T>(n) {
    d = n/2;
    twod = 2*d;
    dp1 = d+1;
    dm1 = d-1;
    dd = d*d;
    initialise();
}

template<class T> funcZSICff<T>::~funcZSICff() { 
    clear(); 
}

template<class T> void funcZSICff<T>::setNumVars(unsigned n) {
    if (n != funcBase<T>::numvars) {
        funcBase<T>::numvars = n;
        clear();
        d = n/2;
        twod = 2*d;
        dp1 = d+1;
        dm1 = d-1;
        dd = d*d;
        initialise();
    }
}

template<class T> void funcZSICff<T>::getFuncGrad(std::vector<T> &x, T &func, std::vector<T> &grad) {}

template<class T> void funcZSICff<T>::getFunc(std::vector<T> &x, T &func) {
  
  typename std::vector<T>::iterator Ix=x.begin(), Iux_j, Iuy_j;  
  
  su = 0;
  for (k=d, Pu_k=&u[0], PZc_jk=&Zc[0]; k>0; --k, ++Pu_k) {
      *Pu_k = 0;
      for (j=d, Iux_j=Ix, Iuy_j=Ix+d; j>0; --j, ++Iux_j, ++Iuy_j, ++PZc_jk) {
          *Pu_k += (*PZc_jk)*(std::complex<T>(*Iux_j, *Iuy_j));
      }
      su += std::real(*Pu_k)*std::real(*Pu_k) + std::imag(*Pu_k)*std::imag(*Pu_k);
  }
  su = std::sqrt(su);
  for (j=0, Pu_j=&u[0], Pu_jpd=&u[0]+d, Pu_jpdpd=&u[0]+twod; j<d; ++j, ++Pu_j, ++Pu_jpd, ++Pu_jpdpd) {
   *Pu_j /= std::complex<T>(su);
   *Pu_jpd = *Pu_j;
   *Pu_jpdpd = *Pu_j;
  }  

  func = 0;
  for (j=0, Pu_j=&u[0], Pu_jpj=&u[0]; j<d; ++j, ++Pu_j, ++Pu_jpj, ++Pu_jpj) {
   Sl = 0;
   for (l=0, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_jpjpl=Pu_jpj; l<d; ++l, ++Pu_l, ++Pu_jpl, ++Pu_jpjpl) {
    Sl += std::conj((*Pu_jpl)*(*Pu_jpl))*(*Pu_l)*(*Pu_jpjpl);
   }
   func += Sl.real()*Sl.real() + Sl.imag()*Sl.imag();
   for (k=j+1, Pu_k=Pu_j+1, Pu_jpk=Pu_jpj+1; k<d; ++k, ++Pu_k, ++Pu_jpk) {
    Sl = 0;
    for (l=0, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_kpl=Pu_k, Pu_jpkpl=Pu_jpk; l<d; ++l, ++Pu_l, ++Pu_jpl, ++Pu_kpl, ++Pu_jpkpl) {
     Sl += std::conj((*Pu_jpl)*(*Pu_kpl))*(*Pu_l)*(*Pu_jpkpl);
    }
    func += 2.0*(Sl.real()*Sl.real()+Sl.imag()*Sl.imag());
   }
  }
  func *= d+1;
  func -= 2;

}

template<class T> void funcZSICff<T>::getFuncGradC(T* Px, T* Pfunc, T* Pgrad) {
  
  su = 0;
  for (k=d, Pu_k=&u[0], PZc_jk=&Zc[0]; k>0; --k, ++Pu_k) {
      *Pu_k = 0;
      for (j=d, Pux_j=Px, Puy_j=Px+d; j>0; --j, ++Pux_j, ++Puy_j, ++PZc_jk) {
          *Pu_k += (*PZc_jk)*(std::complex<T>(*Pux_j, *Puy_j));
      }
      su += std::real(*Pu_k)*std::real(*Pu_k) + std::imag(*Pu_k)*std::imag(*Pu_k);
  }
  su = std::sqrt(su);
  for (j=d, Pu_j=&u[0], Pu_jpd=&u[0]+d, Pu_jpdpd=&u[0]+twod; j>0; --j, ++Pu_j, ++Pu_jpd, ++Pu_jpdpd) {
   *Pu_j /= std::complex<T>(su);
   *Pu_jpd = *Pu_j;
   *Pu_jpdpd = *Pu_j;
  }  

  *Pfunc = 0;
  for (j=d, Pu_j=&u[0], Pu_jpj=&u[0], Ps_jpdj=&s[0]; j>0; --j, ++Pu_j, ++Pu_jpj, ++Pu_jpj, Ps_jpdj+=dp1) {
   Sl = 0;
   for (l=d, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_jpjpl=Pu_jpj; l>0; --l, ++Pu_l, ++Pu_jpl, ++Pu_jpjpl) {
    Sl += std::conj((*Pu_jpl)*(*Pu_jpl))*(*Pu_l)*(*Pu_jpjpl);
   }
   *Ps_jpdj = Sl;
   *Pfunc += Sl.real()*Sl.real() + Sl.imag()*Sl.imag();
   for (k=j-1, Pu_k=Pu_j+1, Pu_jpk=Pu_jpj+1, Ps_jpdk=Ps_jpdj+d, Ps_kpdj=Ps_jpdj+1; k>0; --k, ++Pu_k, ++Pu_jpk, Ps_jpdk+=d, ++Ps_kpdj) {
    Sl = 0;
    for (l=d, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_kpl=Pu_k, Pu_jpkpl=Pu_jpk; l>0; --l, ++Pu_l, ++Pu_jpl, ++Pu_kpl, ++Pu_jpkpl) {
     Sl += std::conj((*Pu_jpl)*(*Pu_kpl))*(*Pu_l)*(*Pu_jpkpl);
    }
    *Ps_jpdk = Sl;
    *Ps_kpdj = Sl;
    *Pfunc += 2.0*(Sl.real()*Sl.real()+Sl.imag()*Sl.imag());
   }
  }
  *Pfunc *= d+1;
  *Pfunc -= 2;

  for (j=d, Pu_j=&u[0], Pu_jpj=&u[0], Pu_dmj=&u[0]+d, Pu_dpdmjmj=&u[0]+twod, Pds__jpdj=&ds[0], Pdsc__jpdj=&dsc[0]; j>0; --j, ++Pu_j, ++Pu_jpj, ++Pu_jpj, --Pu_dmj, --Pu_dpdmjmj, --Pu_dpdmjmj, Pds__jpdj+=dp1, Pdsc__jpdj+=dp1) {
   for (l=d, Pu_l=&u[0], Pu_jpl=Pu_j, Pu_jpjpl=Pu_jpj, Pu_dmjpl=Pu_dmj, Pu_dpdmjmjpl=Pu_dpdmjmj, Pds_l_jpdj=Pds__jpdj, Pdsc_l_jpdj=Pdsc__jpdj; l>0; --l, ++Pu_l, ++Pu_jpl, ++Pu_jpjpl, ++Pu_dmjpl, ++Pu_dpdmjmjpl, Pds_l_jpdj+=dd, Pdsc_l_jpdj+=dd) {
    *Pds_l_jpdj = std::conj((*Pu_jpl)*(*Pu_jpl))*(*Pu_jpjpl) + std::conj((*Pu_dmjpl)*(*Pu_dmjpl))*(*Pu_dpdmjmjpl);
    *Pdsc_l_jpdj = std::complex<T>(2.0)*std::conj((*Pu_dmjpl)*(*Pu_jpl))*(*Pu_l);
   }
   for (k=j-1, Pu_k=Pu_j+1, Pu_jpk=Pu_jpj+1, Pu_dmk=Pu_dmj-1, Pu_dpdmjmk=Pu_dpdmjmj-1, Pu_dmjpk=&u[0]+dp1, Pu_dmkpj=&u[0]+dm1, Pds__jpdk=Pds__jpdj+d, Pds__kpdj=Pds__jpdj+1, Pdsc__jpdk=Pdsc__jpdj+d, Pdsc__kpdj=Pdsc__jpdj+1; k>0; --k, ++Pu_k, ++Pu_jpk, --Pu_dmk, --Pu_dpdmjmk, ++Pu_dmjpk, --Pu_dmkpj, Pds__jpdk+=d, ++Pds__kpdj, Pdsc__jpdk+=d, ++Pdsc__kpdj) {
    for (l=d, Pu_jpl=Pu_j, Pu_kpl=Pu_k, Pu_jpkpl=Pu_jpk, Pu_dmjpl=Pu_dmj, Pu_dmkpl=Pu_dmk, Pu_dpdmjmkpl=Pu_dpdmjmk, Pu_dmjpkpl=Pu_dmjpk, Pu_dmkpjpl=Pu_dmkpj, Pds_l_jpdk=Pds__jpdk, Pds_l_kpdj=Pds__kpdj, Pdsc_l_jpdk=Pdsc__jpdk, Pdsc_l_kpdj=Pdsc__kpdj; l>0; --l, ++Pu_jpl, ++Pu_kpl, ++Pu_jpkpl, ++Pu_dmjpl, ++Pu_dmkpl, ++Pu_dpdmjmkpl, ++Pu_dmjpkpl, ++Pu_dmkpjpl, Pds_l_jpdk+=dd, Pds_l_kpdj+=dd, Pdsc_l_jpdk+=dd, Pdsc_l_kpdj+=dd) {
     *Pds_l_jpdk = std::conj((*Pu_jpl)*(*Pu_kpl))*(*Pu_jpkpl) + std::conj((*Pu_dmjpl)*(*Pu_dmkpl))*(*Pu_dpdmjmkpl);
     *Pdsc_l_jpdk = std::conj((*Pu_dmjpl)*(*Pu_kpl))*(*Pu_dmjpkpl) + std::conj((*Pu_dmkpl)*(*Pu_jpl))*(*Pu_dmkpjpl);
     *Pds_l_kpdj = *Pds_l_jpdk;
     *Pdsc_l_kpdj = *Pdsc_l_jpdk;
    }
   }
  }

  //std::cout << std::endl;
  
  for (j=d, Pdssduc_j=&dssduc[0], Pds_j_k=&ds[0], Pdsc_j_k=&dsc[0]; j>0; --j, ++Pdssduc_j) {
   *Pdssduc_j = 0;
   for (k=dd, Ps_k=&s[0]; k>0; --k, ++Ps_k, ++Pds_j_k, ++Pdsc_j_k) {
    *Pdssduc_j += std::conj(*Pds_j_k)*(*Ps_k) + std::conj(*Pdsc_j_k)*std::conj(*Ps_k);
   }
   //std::cout << *Pdssduc_j << std::endl;
  }

  a = 0;
  for (j=d, Pdssduc_j=&dssduc[0], Pu_j=&u[0]; j>0; --j, ++Pdssduc_j, ++Pu_j) {
   a += (*Pdssduc_j).real()*(*Pu_j).real() + (*Pdssduc_j).imag()*(*Pu_j).imag();
  }

  //std::cout << std::endl;
  
  for (j=d, Pdssdzc_j=&dssdzc[0], Pu_j=&u[0], PZc_jk=&Zc[0]; j>0; --j, ++Pdssdzc_j, ++Pu_j) {
   *Pdssdzc_j = 0;
   for (k=d, Pdssduc_k=&dssduc[0]; k>0; --k, ++Pdssduc_k, ++PZc_jk) {
    *Pdssdzc_j += (*PZc_jk)*(*Pdssduc_k);
   }
   *Pdssdzc_j -= a*(*Pu_j);
   *Pdssdzc_j /= su;
   //std::cout << *Pdssdzc_j << std::endl;
  }
  
  b = 2.0*dp1;

  //std::cout << std::endl;
  
  for (j=d, Pdssdzc_j=&dssdzc[0], Pgradx=Pgrad, Pgrady=Pgrad+d; j>0; --j, ++Pdssdzc_j, ++Pgradx, ++Pgrady) {
   *Pgradx = b*(*Pdssdzc_j).real();
   *Pgrady = b*(*Pdssdzc_j).imag();
   //std::cout << *Pgradx << " " << *Pgrady << std::endl;
  }
  
  
  //std::cout << std::endl << *Pfunc << std::endl;

}


#endif	/* FUNC_H */

