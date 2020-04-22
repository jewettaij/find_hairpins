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
/// @file     bins.cpp
/// @brief    Implements a class for binning points ("atoms") in arbitrary
///           dimensions.  This file borrows some code was from the 
///           "neighbor.cpp" file distributed with LAMMPS.
/// @author   Andrew Jewett, Steve Plimpton, Pieter in 't Veld
/// @license  GPL-2.0


#include <sstream>
using namespace std;
#include "bins.h"



void Bins::
Alloc(Icrd const *aSetNumBins, 
      Iatm set_max_num_atoms)
{
  for (int d = 0; d < g_dim; d++)
    aNbins[d] = aSetNumBins[d] ;

  // The following works, but it uses up too much memory
  // when the atoms are sparsely populated:
  //num_bins=1;
  //for(int d=0; d < g_dim; ++d)
  //  num_bins *= aNbins[d];
  //assert(num_bins > 0);
  //if (num_bins < set_max_num_atoms)
  //  num_bins      = set_max_num_atoms;
  // Instead, just set num_bins = max_num_atoms, and allow multiple
  // lattice sites to share the same bin.
  num_bins      = set_max_num_atoms;
  max_num_atoms = set_max_num_atoms;

  aHeadFromIbin = new Iatm [num_bins];
  aNextFromIatm = new Iatm [max_num_atoms];
  aBinFromIatm  = new Iatm [max_num_atoms];


  if (! aHeadFromIbin) {
    stringstream err_msg;
    err_msg << "Error in memory allocation: file \""
            <<__FILE__<<"\":"<<__LINE__<<"\n"
            << "  Unable to allocate an array of size "<<num_bins<<"\n"<<flush;
    throw InputErr(err_msg.str());
  }


  if (! aNextFromIatm) {
    stringstream err_msg;
    err_msg << "Error in memory allocation: file \""
            <<__FILE__<<"\":"<<__LINE__<<"\n"
            << "  Unable to allocate an array of size "<<max_num_atoms<<"\n"<<flush;
    throw InputErr(err_msg.str());
  }


  for (Ibin ib=0; ib < num_bins; ++ib)
    aHeadFromIbin[ib] = UNINITIALIZED;
  for (Iatm ia=0; ia < max_num_atoms; ++ia)
    aNextFromIatm[ia] = UNINITIALIZED;

} //Bins::Alloc()




void
Bins::Dealloc() {
  //delete [] bboxlo;
  //delete [] bboxhi;
  //delete [] aNbins;
  delete [] aHeadFromIbin;
  delete [] aNextFromIatm;
  delete [] aBinFromIatm;
}





void
Bins::FindMinMax(ConstVect *aaX, Iatm num_atoms)
{
  // Infer the bin array dimensions from atom coordinates and distance cutoff
  assert(num_atoms > 0);
  max_num_atoms = num_atoms;
  for (int ia = 0; ia < num_atoms; ++ia) {
    for (int d = 0; d < g_dim; d++) {
      if ((ia==0) || (aaX[ia][d] < bboxlo[d]))
        bboxlo[d] = aaX[ia][d];
      if ((ia==0) || (aaX[ia][d] > bboxhi[d]))
        bboxhi[d] = aaX[ia][d];
    }
  }
  
} //Bins::MinMax(ConstVect, Iatm)





// infer bin array size and box size from atomic coordinates
Bins::Bins(ConstVect *aaX,
           Iatm num_atoms,
           Real rcut)
{
  FindMinMax(aaX, num_atoms);

  for (int d = 0; d < g_dim; d++) {
    // Find the number of bins in each direction.
    // Add some padding to make sure box width is an integer multiple of rcut.
    //aNbins[d] = static_cast<Icrd>(ceil((bboxhi[d]-bboxlo[d]) / rcut));
    aNbins[d] = static_cast<Icrd>(ceil(0.01 + (bboxhi[d]-bboxlo[d]) / rcut));
    aNbins[d] += 2; // Make sure enough space for a bin on each side of each
                    // atom.  No atom should lie at the boundary of empty space
    Real padding = (rcut * aNbins[d]) - (bboxhi[d] - bboxlo[d]);
    bboxlo[d] -= 0.5 * padding;
    bboxhi[d] += 0.5 * padding;
    bininv[d] = aNbins[d] / (bboxhi[d] - bboxlo[d]);
  }

  // Alloc tables for the bins
  Alloc(aNbins, num_atoms);

  // Now assign each atom to a bin
  BinAtoms(aaX, num_atoms);
}




// manually specify bin array size
Bins::Bins(ConstVect *aaX,
           Iatm num_atoms,
           Real rcut,
           Icrd const *aSetNumBins)
{
  assert(aSetNumBins);
  FindMinMax(aaX, num_atoms);
  
  for (int d = 0; d < g_dim; d++) {
    // Find the number of bins in each direction.
    if ((aSetNumBins[0] > 0) &&
        ((bboxhi[d] - bboxlo[d])/aSetNumBins[d]) >= rcut)
      aNbins[d] = aSetNumBins[d];
    else {
      //aNbins[d] = static_cast<Icrd>(ceil((bboxhi[d]-bboxlo[d]) / rcut));
      aNbins[d] = static_cast<Icrd>(ceil(0.01 + (bboxhi[d]-bboxlo[d]) / rcut));
    }
    aNbins[d] += 2; // Make sure enough space for a bin on each side of each
                    // atom.  No atom should lie at the boundary of empty space
    // Add some padding to make sure box width is an integer multiple of rcut.
    Real padding = (rcut * aNbins[d]) - (bboxhi[d] - bboxlo[d]);
    if (padding > 0.0) {
      bboxlo[d] -= 0.5 * padding;
      bboxhi[d] += 0.5 * padding;
      bininv[d] = aNbins[d] / (bboxhi[d] - bboxlo[d]);
    }
  }


  // Alloc tables for the bins
  Alloc(aNbins, num_atoms);

  // Now assign each atom to a bin
  BinAtoms(aaX, num_atoms);
} //Bins::Bins(ConstVect, Iatm, Real)



// manually specify bin array size and box boundaries
// I actually don't use this constructor, so perhaps I should get rid of it.
Bins::Bins(Icrd const *aSetNumBins, 
           Iatm set_max_num_atoms,
           ConstVect set_bboxlo,
           ConstVect set_bboxhi,
	   Real rcut)
{
  assert(aSetNumBins && set_bboxlo && set_bboxhi);

  for (int d=0; d < g_dim; ++d)
  {
    bboxlo[d] = set_bboxlo[d];
    bboxhi[d] = set_bboxhi[d];
    if ((aSetNumBins[0] > 0) &&
        ((bboxhi[d] - bboxlo[d])/aSetNumBins[d]) >= rcut)
      aNbins[d] = aSetNumBins[d];
    else {
      //aNbins[d] = static_cast<Icrd>(ceil((bboxhi[d]-bboxlo[d]) / rcut));
      aNbins[d] = static_cast<Icrd>(ceil(0.01+(bboxhi[d]-bboxlo[d]) / rcut));
    }
    aNbins[d] += 2; // Make sure enough space for a bin on each side of each
                    // atom.  No atom should lie at the boundary of empty space
    // Add some padding to make sure box width is an integer multiple of rcut.
    Real padding = (rcut * aNbins[d]) - (bboxhi[d] - bboxlo[d]);
    if (padding > 0.0) {
      bboxlo[d] -= 0.5 * padding;
      bboxhi[d] += 0.5 * padding;
      bininv[d] = aNbins[d] / (bboxhi[d] - bboxlo[d]);
    }
  } //for (int d=0; d < g_dim; ++d)

  Alloc(aNbins, set_max_num_atoms);
}




Bins::~Bins() {
  Dealloc();
}



