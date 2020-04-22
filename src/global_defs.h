#ifndef _GLOBAL_DEFS_H
#define _GLOBAL_DEFS_H

#include<iostream>
#include<string>
using namespace std;

//#define NDEBUG         //<--dissables assert()

// (A constant named "G_DIM" is often included which specifies the number of
//  dimensions that atoms are allowed to move in (usually 3).  Any positive
//  integer is allowed, but it usually must be defined at compile time.)
#ifndef G_DIM
const int g_dim = 3; // <- The number of coordinates for each atom.
#else
const int g_dim = G_DIM;
#endif


typedef double Real;
typedef Real Vect[g_dim];
typedef Real const ConstVect[g_dim];


template<class _Real>
inline _Real
DotProduct(_Real const *A, _Real const *B)
{
  _Real AdotB = 0.0;
  for (int d=0; d < g_dim; ++d)
    AdotB += A[d]*B[d];
  return AdotB;
}

template<class _Real>
inline _Real
DistanceSqd(_Real const *A, _Real const *B)
{
  _Real  dsqr= 0.0;
  for (int d=0; d < g_dim; ++d) {
    Real xd = A[d]-B[d];
    dsqr += xd*xd;
  }
  return dsqr;
}


//periodic boundary conditions:
template<class X>
inline X PeriodicImage(X i, X a, X b) {
  if (i < a)
    return i + (b-a);
  if (i >= b)
    return i - (b-a);
  else 
    return i;
}


//non-periodic boundary conditions:
template<class X>
inline X SaturateImage(X i, X a, X b) {
  if (i < a) {
    cerr << i << " not in [" << a << "," << b << "]" << endl;
    return a;
  }
  if (i >= b) {
    cerr << i << " not in [" << a << "," << b << "]" << endl;
    return b-1;
  }
  else 
    return i;
}


class InputErr {
  string msg;
public:
  InputErr(const char *description):msg(description) {}
  InputErr(string description):msg(description) {}
  virtual const char *what() const throw() {return msg.c_str();}
};



class ArgParseErr: public InputErr {
public:
  ArgParseErr(const char *description):InputErr(description) {}
  ArgParseErr(string description):InputErr(description) {}
};



#endif //#ifndef _GLOBAL_DEFS_H
