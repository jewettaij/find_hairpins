/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/// @file     bins.h
/// @brief    Defines a class for binning points ("atoms") in arbitrary
///           dimensions.  This file borrows some code was from the
///           "neighbor.cpp" file distributed with LAMMPS.
/// @author   Andrew Jewett, Steve Plimpton, Pieter in 't Veld
/// @license  GPL-2.0

#ifndef _BINS_H
#define _BINS_H

#include <cmath>
#include <cassert>
using namespace std;
#include "global_defs.h"



typedef long Iatm;
typedef long Ibin;
typedef int  Icrd;
const   int  UNINITIALIZED = -1;
const   Icrd ICRD_UNINITIALIZED = -255; //any large negative integer will do



class Bins {

 public:
  Icrd  aNbins[g_dim]; // number of bins in each dimension
  Ibin  num_bins;      // total number of bins.
                       //  (Originally I used to set it to
                       //   = aNbins[0] * aNbins[1] * aNbins[2] ...
                       //  ...but to save space, I currently just 
                       //  set it equal to the number of atoms.)
  Vect  bboxlo;        // Boundary box minimum
  Vect  bboxhi;        // Boundary box maximum

  Iatm *aHeadFromIbin; // ptr to 1st atom in each bin
  Iatm *aNextFromIatm; // ptr to next atom in each bin
                       // These two arrays implement a linked-list 
                       // of atoms for each bin (terminated by -1).
  Ibin *aBinFromIatm;  // which bin is this atom in? ( =  IbinFromX(aaX[ia]))

            // Example: In 3-dimensions, the first atom in the bin at
            // location ix,iy,iz (integer coordinates) is stored in:
            // The next atom (in this bin) is at aNextFromAtomID[ia] where
            // ia = aHeadFromIbin[ib]    and, in the naive implementation:
            // ib = aNbins[1]*aNbins[0]*iz + aNbins[0]*iy + ix

 public:

  // infer bin array size and box size from atom coordinates and distance cutoff
  Bins(ConstVect *aaX,
       Iatm num_atoms,
       Real rcut); 


  // manually specify bin array size
  Bins(ConstVect *aaX,
       Iatm num_atoms,
       Real rcut,
       Icrd const *aSetNumBins);


  // manually specify bin array size and box boundaries
  Bins(Icrd const *aSetNumBins, 
       Iatm set_max_num_atoms,
       ConstVect set_bboxlo,
       ConstVect set_bboxhi,
       Real rcut);


  ~Bins();


  inline void
  IcrdFromX(Icrd *aIcrd, ConstVect aX) const
  {
    assert(aX && aIcrd);
    for (int d=0; d < g_dim; ++d) {
      aIcrd[d] = static_cast<Icrd>(floor((aX[d]-bboxlo[d])*bininv[d]));
      if (aIcrd[d] >= aNbins[d])
        aIcrd[d] = aNbins[d]-1;
    }
  }


  inline Ibin 
  IbinFromIcrd(Icrd const *aIcrd) const
  {
    // In 3-D, this function simply returns 
    //   aNbins[1]*aNbins[0]*iz + aNbins[0]*iy + iz
    // where ix = aIcrd[0]
    //       iy = aIcrd[1]
    //       iz = aIcrd[2]
    assert(aIcrd);

    Ibin ib = aIcrd[g_dim-1];
    assert((0 <= ib) && (ib < aNbins[g_dim-1]));
    //cerr << aIcrd[0] << " " << aIcrd[1] << " " << aIcrd[2] << endl;

    for(int d = g_dim-2; d >= 0; --d)
    {
      //cerr << "d="<< d <<", ib="<<ib << endl;
      assert((0 <= aIcrd[d]) && (aIcrd[d] < aNbins[d]));
      ib *= aNbins[d];
      ib += aIcrd[d];
    }

    return ib;
    //return ib % num_bins;


    // Alternate method: precoumpute aNumBinsInDimension[] array and then use:
    //Icrd i = aIcrd[0];
    //for(int d=1; d < g_dim; ++d)
    //  i += aNumBinsInDimension[d] * aIcrd[d];
  }


  // IbinFromX(): Given the floating-point coordinates, aX[], find the bin #.
  //              aX should lie in the boundary box (bboxlo[] and bboxhi[])
  inline 
  Ibin IbinFromX(ConstVect aX) const
  {
    // In 3-D, this function simply returns 
    //   aNbins[1]*aNbins[0]*iz + aNbins[0]*iy + iz
    int  d = g_dim - 1;
    Ibin ib = static_cast<Icrd>(floor((aX[d]-bboxlo[d])*bininv[d]));
    assert((0 <= ib) and (ib < aNbins[d]));
    for (d=g_dim-2; d >= 0; --d)
    {
      Icrd ix = static_cast<Icrd>(floor((aX[d]-bboxlo[d])*bininv[d]));
      assert((0 <= ix) and (ix < aNbins[d]));
      ib *= aNbins[d];
      ib += ix;
    }
    //assert((0 <= ib) && (ib < num_bins));
    return ib;
    //return ib % num_bins;
  }


  inline void 
  IcrdFromIbin(Icrd *aIcrd, Ibin ib) const
  {
    assert(aIcrd);
    Ibin next_ib;
    Icrd size;
    for(int d=0; d < (g_dim-1); ++d)
    {
      size = aNbins[d];
      next_ib  = ib / size;
      aIcrd[d] = ib % size;  
      //aIcrd[d] = ib - size*next_ib;  alternate version

      assert((0 <= aIcrd[d]) && (aIcrd[d] < aNbins[d]));
      ib = next_ib;
    }
    aIcrd[g_dim-1] = ib;
    assert((0 <= aIcrd[g_dim-1]) && (aIcrd[g_dim-1] < aNbins[g_dim-1]));
  }



  // BinAtoms():  Put all of the atoms into bins according to their position.
  // (This function was stolen from LAMMPS' Neighbor::bin_atoms())

  inline 
  void BinAtoms(ConstVect *aaX,  // <--> equivalent to "double (*x)[g_dim],"
                Iatm num_atoms)
  {
    assert(aHeadFromIbin && aNextFromIatm && aBinFromIatm);
    // Bin in reverse order so linked list will be in forward order.
    for (Iatm ia = num_atoms-1;  ia >= 0;  ia--) {
      // Original version assumes num_bins = aNbins[0]*aNbins[1]*aNbins[2]...
      // In that case, use this code:
      //
      //Ibin ib  = IbinFromX(aaX[ia]);
      //aNextFromIatm[ia]  = aHeadFromIbin[ib];
      //aHeadFromIbin[ib] = ia;
      //
      // New version tries to be more memory efficient (if occasionally slower)
      // by allowing multiple different lattice sites to use the same bin
      // in memory. To save memory, we do not allow num_bins to exceed num_atoms
      // and use modulo arithmetic to make sure the index ("ibb") into the
      // aHeadFromIbin[] array does not exceed num_bins.
      Ibin ib  = IbinFromX(aaX[ia]);
      Ibin ibb = ib % num_bins; 
      aNextFromIatm[ia]  = aHeadFromIbin[ibb];
      aHeadFromIbin[ibb] = ia;
      // Inverse lookup:
      aBinFromIatm[ia] = ib;
    }
  }



  // Find the coordinates of the next neighbor bin (voxel) adjacent to aIcrd.
  // Store the result in aJcrd.  Return true if the operation was successful.
  // (Return false if we have looped past the end of the available bins.)
  // If the coordinates of aJcrd have not yet been specified, then set the
  // first entry in aJcrd[0] argument to ICRD_UNINITIALIZED.
  inline bool
  NextNeighborBin(Icrd const *aIcrd,
		  Icrd *aJcrd,
		  bool first=false)
  {
    assert(aIcrd && aJcrd);

    if (aJcrd[0] == ICRD_UNINITIALIZED) {
      for(int d=0; d < g_dim; ++d)
        aJcrd[d] = aIcrd[d]-1;
    }
    // Now, "increment" the position of the cell (aJcrd[])
    // This is just like incrementing a g_dim - digit number:
    bool carry_the_1 = true;
    int d=0;
    while (carry_the_1 && (d < g_dim))
    {
      // For periodic boundary conditions, uncomment these lines:
      //aJcrd[d] = PeriodicImage(aJcrd[d]+1, 0, aNbins[d]);
      aJcrd[d] = aJcrd[d]+1;
      //if (aJcrd[d] == PeriodicImage(aIcrd[d]+2, 0, aNbins[d]))
      if (aJcrd[d] == aIcrd[d]+2)
      {
        //aJcrd[d] = PeriodicImage(aIcrd[d]-1, 0, aNbins[d]);
        aJcrd[d] = aIcrd[d]-1;
        ++d;
      }
      else
        carry_the_1 = false;
    }
    return (d < g_dim);
  } //NextNeighborBin()


#if 0
  // COMMENTING THIS VERSION OUT

  inline Ibin
  NextNeighborBin(Ibin ib, Ibin jb)
  {
    // There is slow funky to loop through the neighbor bins (voxels).
    // (It would be much easier to do this in 3-D instead of N-D.)
    // I could do this faster (without invoking IcrdFromIbin())
    // but the code would be uglier.  Perhaps I'll try this later.
    Icrd aIcrd[g_dim];
    Icrd aJcrd[g_dim];
    IcrdFromIbin(aIcrd, ib);

    if (jb == UNINITIALIZED) {
      for(int d=0; d < g_dim; ++d) {
        aJcrd[d] = aIcrd[d]-1;
        //aJcrd[d] = PeriodicImage(aIcrd[d]-1, 0, aNbins[d]);
      }
      return IbinFromIcrd(aJcrd);
    }
    else {
      IcrdFromIbin(aJcrd, jb);
      // Now, "increment" the position of the cell (aJcrd[])
      // This is just like incrementing a g_dim - digit number:
      bool carry_the_1 = true;
      int d=0;
      while (carry_the_1 && (d < g_dim))
      {
        //aJcrd[d] = PeriodicImage(aJcrd[d]+1, 0, aNbins[d]);
        aJcrd[d] = aJcrd[d]+1;
        //if (aJcrd[d] == PeriodicImage(aIcrd[d]+2, 0, aNbins[d]))
	if (aJcrd[d] == aIcrd[d]+2)
        {
          //aJcrd[d] = PeriodicImage(aIcrd[d]-1, 0, aNbins[d]);
          aJcrd[d] = aIcrd[d]-1;
          ++d;
        }
        else
          carry_the_1 = false;
      }
      if (d >= g_dim)
        return UNINITIALIZED;
      else {
        return IbinFromIcrd(aJcrd);
      }
    }
  } //NextNeighborBin()

#endif //#if 0


 private:

  void Alloc(Icrd const *aSetNumBins, Iatm set_max_num_atoms);
  void Dealloc();
  void FindMinMax(ConstVect *aaX, Iatm num_atoms);

  Ibin max_num_atoms;  //upper limit on the number of atoms which will be binned
  Vect bininv;         //  = aNBins[d] / (bboxhi[d]-bboxlo[d])

  // No need yet for a copy constructor.
  // (I think...  To make sure, dissable it.)
  Bins(const Bins &copy);
  Bins& operator = (const Bins &copy);

}; //class Bins


#endif //#ifndef _BINS_H
