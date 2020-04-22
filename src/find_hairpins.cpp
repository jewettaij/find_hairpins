/// @file     find_hairpins.cpp
/// @brief    Given a list of points in space which are assumed to be
///           distributed unifornly along a curve in space, find the
///           locations of every hairpin.  This program can search for
///           straight hairpins in O(n) time, and branched hairpins in O(n^2)
///           time.  (The general case with gaps requires O(n*n1*n2) time.)
/// @author   Andrew Jewett, Steve Plimpton
/// @license  GPL-2.0
///
/// @note  To change the number of dimensions in which the point cloud lives,
/// (ie, the number of coordinates on each line of the structure file)
/// then compile using following compiler flag "-DG_DIM=4" (4 dimensions).
/// (You can also edit the "global_defs.h" file.)


#include <cassert> //needed for assert()
#include <cstring> //needed for strcmp()
#include <cstdlib> //needed for atol()
#include <set>
#include <vector>
using namespace std;


#include "global_defs.h"
#include "bins.h"
#include "io.h"


// global variables
const char g_program_name[]   = "find_hairpins";
const char g_version_string[] = "0.0.2";
const char g_date_string[]    = "<2020-4-21>";




struct Settings
{

  enum CoordFileFormat {COORD_FORMAT_RAW,
                        COORD_FORMAT_XYZ};

  CoordFileFormat coord_file_format;
  Real rcut; //cutoff distance for whether two atoms are "in contact"
  Iatm n1;
  Iatm n2;
  Iatm ni;
  Icrd aUserNumBins[g_dim]; //let the user choose the number of bins

  Settings() {
    coord_file_format = COORD_FORMAT_RAW;
    rcut = 1.5;
    n2 = 1;
    n1 = 1;
    ni = 3;
    for (int d = 0; d < g_dim; d++)
      aUserNumBins[d] = UNINITIALIZED;
  }

}; //struct Settings



void PrintHairpins(const set<pair<Iatm,Iatm> > &hairpins,
                   ostream& out_file)
{
  for(auto pl = hairpins.begin(); pl != hairpins.end(); pl++)
    out_file << (*pl).first << " " << (*pl).second << "\n";
}



void FindHairpins(const vector<set<Iatm> >& m,
                  set<pair<Iatm,Iatm> > &hairpins,
                  Iatm n1 = 1,
                  Iatm n2 = 1,
                  Iatm ni = 3)
{
  if (n1 == 0)
    n1 = m.size(); // (This effectively sets it to infinity, which means that
                   //  bulges of any size on one side are allowed.  Use this if
                   //  you want to find the total length of a hairpins including
                   //  any branches.  Otherwise set n1 = n2 if only want to look
                   //  for unbrached hairpins.)
  for (Iatm i0 = 0; i0 < m.size(); i0++) {
    bool continue_iter = true;
    Iatm ia = i0;
    Iatm ib = i0;
    while (continue_iter) {
      Iatm sa = ia;
      Iatm sb = ib;
      continue_iter = false;

      // Suppose we know that a hairpin contains a pair of points (with indices
      // ia,ib where ia<=ib) which are in close proximity to eachother.  This 
      // means all of the points in the interval [ia, ib] belong to the hairpin.
      // Now expand the hairpin outwards to include points on either side of the
      // interval [ia,ib].  We cannot assume that the hairpin is ideal.  Even
      // if ia-1 and ib+1 also belong to the hairpin, we cannot assume they are
      // inclose proximity to each other because the curve on either side of
      // the hairpin could bulge or make temporary excursions away from its
      // partner.  (Later on, I refer to these as gaps, branches or bulges).
      // We have to look further outward.
      // 
      // Figure: Diagram showing possible pairs of points in the 
      //         point cloud considered during a typical iteration.
      //
      // (@ denotes the hairpin discovered so far
      //  * denotes the pairs searched in the iteration starting from ia,ib)
      //
      //   sb
      //  axis 
      //  /|\
      //   |          
      //   |          d  _ 
      //   |            |`. 
      //   |    ^          `.          * * * * *
      //   |    |            `.---> c  * * * * *
      //   |    |             |`.      * * * * *
      //   |    |             v  `.    * * * * *
      //   |    |             c    `   * * * * *
      //   |    n1                     * * * * *
      //   |    |     ^    * * * * * * * * * * *
      //   |    |     |    * * * * * * * * * * *
      //   |    |    n2    * * * * * * * * * * *
      //   |    |     |    * * * * * * * * * * *             current known
      //   |    v     v    * * * * * * * * * * @ <-- (ia,ib) hairpin interval
      //   |                                     @           
      //   |                                       @
      //   |                                         @                 
      //   |                                           @               hairpin
      //   |                                             @ <-- (i0,i0) tip here
      //   |
      //   |
      //   |
      //   |_________________________________________________________\ sb
      //                                                             / axis
      // Algorithm:
      // Considerthe pairs of points which might belong to the hairpin.
      // Starting from pair ia,ib, at the next iteration we search over
      // over the pairs in the "solid" region in sa,sb space (denoted with
      // "*" in the figure above).  If ANY of those pairs of points lie within
      // the cutoff distance from each other (ie. if they are in "contact"),
      // then the hairpin expands to include them (ia=sa, and ib=sb).
      // Then the iteration continues until no such pair exists, or we
      // run out of points.  The final pair denotes the interval [ia, ib]
      // demarcating the beginning and ending of the hairpin.
      //
      // Search order:
      //   I choose to give preference to pairs that:
      // 1) minimize the larger gap size on either side (max(ia-sa, sb-ib)
      // 2) minimize the difference in gap size |(ia-sa)-(sb-ib)|

      // d is distance along along the diagonal (see diagram)
      for (Iatm d = 1;
           (! continue_iter) && (d <= n1);
           d++)
      {
        // c is distance away from the diagonal (see diagram)
        for (Iatm c = std::max(d - n2, static_cast<Iatm>(0));
             (! continue_iter) && (c <= d);
             c++)
        {
          // Check for pairs of points on either side of the "d" diagonal.
          // First check to the left of the "d" axis.  (See diagram.)
          sa = ia - (d - c);
          sb = ib + d;
          if ((0 <= sa) && (sb < m.size())) {
            if (m[sa].find(sb) != m[sa].end()) {
              continue_iter = true;
              break;
            }
          }
          // Then check to the right of the "d" axis.  (See diagram.)
          // (Note: Choosing left before right was arbitrary.)
          sa = ia - d;
          sb = ib + (d - c);
          if ((0 <= sa) && (sb < m.size())) {
            if (m[sa].find(sb) != m[sa].end()) {
              continue_iter = true;
              break;
            }
          }
        } // loop over c
      } // loop over d
      if (continue_iter) {
        ia = sa; // get ready for the next iteration
        ib = sb;
      }
    } //while (continue_iter)
    if (1 + ib - ia >= ni)
      hairpins.insert(pair<Iatm,Iatm>(ia, ib));
    if (i0+1/100 != i0/100)
      cerr << "    progress: " << 100*(i0+1)/m.size() << "%" << endl;
  } //for (Iatm i0 = 0; i0 < m.size(); i0++)
} // FindHairpins




// THIS FUNCTION IS INEFFICIENT.  I SHOULD REWRITE IT
void 
CalcAdjacency(vector<set<Iatm> >& m,
              ConstVect *aaX,  // <--> equivalent to "double (*x)[g_dim],"
              Iatm num_atoms,
              Real rcut,
              Bins &bins)
{
  Icrd aIcrd[g_dim];
  Icrd aJcrd[g_dim];
  Real rcutsq = rcut*rcut;

  for (Iatm ia = 0;  ia < num_atoms;  ia++)
  {
    Ibin ib = bins.aBinFromIatm[ia]; //the bin-id of bin containing atom ia
    // Find the array of coordinates (aIcrd) for that bin (ib):
    bins.IcrdFromIbin(aIcrd, ib);

    aJcrd[0] = ICRD_UNINITIALIZED; // <-signal this is the first neighbor bin.
    while (bins.NextNeighborBin(aIcrd, aJcrd))
    {
      // loop over all of the atoms in that bin
      Ibin jb  = bins.IbinFromIcrd(aJcrd);
      // original version:
      //Iatm ja = bins.aHeadFromIbin[jb]; // "ja" is the atom-id of a neighbor
      //                                  // candidate located in bin "jb"
      // memory efficient version:
      Ibin jbb = jb % bins.num_bins;
      Iatm ja  = bins.aHeadFromIbin[jbb];

      while (ja >= 0)
      {
        if ((DistanceSqd(aaX[ia], aaX[ja]) < rcutsq) && (ia != ja))
          m[ia].insert(ja);
        ja = bins.aNextFromIatm[ja];
      }
    }
  } //for (Iatm ia = 0;  ia < num_atoms;  ia++)
} //CalcAdjacency()






void 
ParseArgs(int argc,
          char **argv,
          Settings &settings)
{

  stringstream explanation;
  explanation 
    <<"\n"
    "Explanation:\n"
    "\n"
    "  This program reads a coordinate file containing one (or more) structures\n"
    "  (in .raw or .xyz format), and generates a sparse matrix of contacts\n"
    "  i1 j2 count1\n"
    "  i2 j2 count2\n"
    "  i3 j3 count3\n"
    "   :  :   :\n"
    "  where \"i\" and \"j\" refer to monomers in the structure(s), and\n"
    "  \"count\" refers to the number of times those two monomers were found\n"
    "  to be in contact with eachother in this structure(s).\n"
    "    (\"in contact\" <--> to be spatially separated by a distance <= rcut)\n"
    "\n"
    "Typical Usage:\n"
    "\n"
    "\n"
    <<"  "<< g_program_name<<" -r rcut [optional arguments..] < coordinate_file > matrix_file\n"
    "\n"
    "\n"
    "Optional arguments:\n"
    "\n"
    "  -r rcut      <- specify the threshold contact distance (1.5 by default).\n"
    "\n"
    "  -n1 n1       <-can skip n1 points in either side of the loop\n"
    "                 without discarding the hairpin\n"
    "  -n2 n2       <-can skip n2 points in both sides of the loop\n"
    "                 without discarding the hairpin\n"
    "  -ni ni       <-ignore loops whose size is less than ni\n"
    "  -raw         <-The default input file format: 3-column space-\n"
    "                 delimited text file with blank lines separating frames.\n"
    "                 Note: \".RAW\" format is assumed by default.\n"
    "  -xyz         <-Instead, assume the input coordinate file uses \".XYZ\" format.\n"
    "                 (This feature is not well tested.  Hopefully it works.)\n"
    "  -bins N      <-Specify the number of bins (aka voxels) in each direction.\n"
    "                 Alternatley, you can specify 3 numbers (in 3 dimensions)\n"
    "                 if you want it to vary in different directions.\n"
    "                 (This option has not been tested and may not work. 2013-60-6)\n"
    "\n";

  try {
    int which_arg = 1;
    while (which_arg < argc)
    {
      if (strcmp(argv[which_arg], "-raw") == 0)
      {
        settings.coord_file_format = Settings::COORD_FORMAT_RAW;
        which_arg += 1;
      }
      else if (strcmp(argv[which_arg], "-xyz") == 0)
      {
        settings.coord_file_format = Settings::COORD_FORMAT_XYZ;
        which_arg += 1;
      }
      else if ((strcmp(argv[which_arg], "-r") == 0) ||
	       (strcmp(argv[which_arg], "-d") == 0) ||
	       (strcmp(argv[which_arg], "-rcut") == 0))
      {
        if (which_arg+1 == argc) { 
          stringstream errmsg;
          errmsg << "\nError: expected a number following \""<<argv[which_arg]
                 << "\"\n";
          throw ArgParseErr(errmsg.str());
        }
        settings.rcut = atof(argv[which_arg+1]);
        which_arg += 2;
      }
      else if (strcmp(argv[which_arg], "-n1") == 0)
      {
        if (which_arg+1 == argc) { 
          stringstream errmsg;
          errmsg << "\nError: expected an integer following \""<<argv[which_arg]
                 << "\"\n";
          throw ArgParseErr(errmsg.str());
        }
        settings.n1 = atoi(argv[which_arg+1]) + 1;
        which_arg += 2;
      }
      else if (strcmp(argv[which_arg], "-n2") == 0)
      {
        if (which_arg+1 == argc) { 
          stringstream errmsg;
          errmsg << "\nError: expected an integer following \""<<argv[which_arg]
                 << "\"\n";
          throw ArgParseErr(errmsg.str());
        }
        settings.n2 = atoi(argv[which_arg+1]) + 1;
        which_arg += 2;
      }
      else if (strcmp(argv[which_arg], "-ni") == 0)
      {
        if (which_arg+1 == argc) { 
          stringstream errmsg;
          errmsg << "\nError: expected an integer following \""<<argv[which_arg]
                 << "\"\n";
          throw ArgParseErr(errmsg.str());
        }
        settings.ni = atoi(argv[which_arg+1]);
        which_arg += 2;
      }
      else if (strcmp(argv[which_arg], "-bins") == 0)
      {
        if (which_arg+1 == argc) { 
          stringstream errmsg;
          errmsg << "\nError: expected at least one number following \""<<argv[which_arg]
                 << "\"\n";
          throw ArgParseErr(errmsg.str());
        }

        vector<Icrd> vUserNumBins; //<-Store numbers here while reading arg list
        int j = which_arg+1;
        while ((j < argc) && 
               (strlen(argv[j]) != 0) &&
               isdigit(argv[j][0]))
        {
          vUserNumBins.push_back(atoi(argv[j]));
          j++;
        }
        // Then copy the numbers into the settings.aUserNumBins[] array:
        if (vUserNumBins.size() == 1)
          // That means use the same number of bins in all 3 directions
          for (int d = 0; d < g_dim; ++d) 
            settings.aUserNumBins[d] = vUserNumBins[0];
        else if (vUserNumBins.size() == g_dim) {
          for (int d = 0; d < g_dim; ++d) 
            settings.aUserNumBins[d] = vUserNumBins[d];
        }
        else {
          stringstream errmsg;
          errmsg << "Error: Wrong number of integers ("
                 << vUserNumBins.size()
                 << ")\n"
                 << "       following \""<<argv[which_arg]<<"\" argument.\n"
                 << "       Expected either 1 or "<<g_dim<<" integers.\n";
          throw ArgParseErr(errmsg.str());
        }
        which_arg = j;
      }
      else
      {
        stringstream errmsg;
        errmsg << "\nError: Unrecognized argument: \""<<argv[which_arg]<<"\"\n";
        throw ArgParseErr(errmsg.str());
      }
    } //while (which_arg < argc)
    if ((settings.n2 > settings.n1) && (settings.n1 != 0)) {
      stringstream errmsg;
      errmsg << "\nError: The n2 parameter (currently "<<settings.n2
             <<") must not exceed n1 (currently "<<settings.n1<<").\n"
             << "(Use the \"-n1\" and \"-n2\" arguments to specify these parameters.)\n";
      throw ArgParseErr(errmsg.str());
    }
  }
  catch (ArgParseErr& e)
  {
    cerr << "--------------------------------------------------\n";
    cerr << "       Error in argument list: \n" << e.what() << endl;
    cerr << "--------------------------------------------------\n";
    cerr << "\n";
    //cerr << explanation.str() << endl;
    exit(-1);
  }

  //---- finished parsing the argument list ----

} //ParseArgs()






int
main(int argc, char **argv)
{

  //---- Load the coordinate file into a large buffer ----

  try 
  {
    long long line_count = 1; //counts all lines including blanks and comments

    vector<set<Iatm> > m_tot; //adjacency matrix

    cerr << g_program_name   << ", v"
         << g_version_string << " "
         << g_date_string    << "\n";

    // Process the argument list
    Settings settings;
    ParseArgs(argc, argv, settings);

    cerr << "  reading input file...\n";

    long long frame_counter = 0;

    // Read one structure from the file and calculate its adjacency matrix.
    while(cin)
    {
      Vect *aaX = NULL;

      // Read in the next frame
      Iatm num_atoms = UNINITIALIZED; //number of atoms in each snapshot/frame

      if (settings.coord_file_format == Settings::COORD_FORMAT_RAW)
        num_atoms = ReadCoordsRAW(cin, &aaX, line_count);
      else if (settings.coord_file_format == Settings::COORD_FORMAT_XYZ)
        num_atoms = ReadCoordsXYZ(cin, &aaX, line_count);
      else
        assert(0);

      if (num_atoms > 0)
        cerr << "    finished reading frame " << frame_counter+1 << endl;

      // We could be at a point in the file with trailing whitespace.  Check
      // to make sure that we did actually read something before we continue.

      if (num_atoms > 0) {

        Bins bins(aaX, 
                  num_atoms, 
                  settings.rcut,
                  settings.aUserNumBins);

        vector<set<Iatm> > m(num_atoms); // adjacency matrix (sparse)
        m_tot.resize(std::max(m.size(), m_tot.size()));

        //    Note: If you are unfamilliar with C++ "maps" or "pairs", 
        //          they are similar to python dictionaries and tuples.  See:
        //    http://www.cplusplus.com/reference/map/map/
        //    http://www.cplusplus.com/reference/utility/pair/

        CalcAdjacency(m,
                      aaX,
                      num_atoms,
                      settings.rcut,
                      bins);

        // Add these contacts to the total contact matrix m_tot
        for (long i = 0; i < num_atoms; i++)
          for (auto pj = m[i].begin(); pj != m[i].end(); pj++)
            m_tot[i].insert(*pj);

        frame_counter++;
        cerr << "    finished calculating adjacency matrix for frame " << frame_counter << endl;
      }
      else
        assert(! cin);

      delete [] aaX;

    } //while(cin)


    if (frame_counter > 0) {
      // Now print the results
      set<pair<Iatm,Iatm> > hairpins;
      cerr << "  Finding hairpins:" << endl;
      FindHairpins(m_tot,
                   hairpins,
                   settings.n1,
                   settings.n2,
                   settings.ni);

      cerr << "  Printing list of hairpins:" << endl;
      PrintHairpins(hairpins, cout);
    }

  } //try
  catch (InputErr& e)
  {
    cerr << "       Error in input file(s): \n" << e.what() << endl;
    exit(-1);
  }

} //main()

