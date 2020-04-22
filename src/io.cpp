#include <cassert>
#include <cstdlib>
using namespace std;
#include "io.h"

void
SkipWhiteSpace(istream &in,
               long long& line_count,
               char comment_char)
{
  if (! in)
    throw InputErr("(SkipWhiteSpaces)\nError reading file: File ends prematurely.");
  char c;
  while(in.get(c))
  {
    if ((! isspace(c)) && in)
    {
      in.putback(c);
      break;
    }
  }
}

/*
void
SkipWhiteSpace(istream& in,
               long long& line_count,
               char comment_char='#')
{
  char c;
  while(in.get(c))
  {
    if ((! isspace(c)) && in)
    {
      in.putback(c);
      break;
    }
  }
}
*/

/*
void
SkipWhiteSpace(istream& in,
               long long& line_count,
               char comment_char='#')
{
  char c='\0';
  bool commented_state = false;
  bool skipover_state = true;
  while(in.get(c))
  {
    if (c == comment_char)
      commented_state = true;
    if ((! isspace(c)) && in)
    {
      in.putback(c);
      break;
    }
  }
}
*/

/*
This version does not carry out the in.putback(c) command for some reason.
void SkipWhiteSpace(istream& in,
                    long long& line_count,
                    char comment_char='#')
{
  char c = ' ';
  bool commented_state = false;
  bool skipover_state = true;
  while (in.get(c) && skipover_state) {
    if (c == comment_char)
      commented_state = true;
    else if (c == '\n') {
      line_count++;
      commented_state = false;
    }
    else if (! (isspace(c) || commented_state)) {
      skipover_state = false;
      in.putback(c);
    }
  }
}
*/







// ReadCoordsRAW()
// This is a somewhat complicated function for reading in coordinates 
// from multi-column ascii files with the following format:
//
// x1 y1 z1   #
// x2 y2 z2   #
// x3 y3 z3   # frame 1
//  :  :  :   #
// xN yN zN   #
//
// x1 y1 z1   #
// x2 y2 z2   #
// x3 y3 z3   # frame 2
//  :  :  :   #
// xN yN zN   #
//
//  :  :  :
//
// I can't simply use use >> to read numbers from the file.
// It's complicated because it skips over single-line comments
// (beginning with a '#', by default).
// It can also read coordinate files with more/less than 3 numbers on each line.
// (I have some 2D simulation data.)
// Blank lines (or lines containing comments)
// are used to separate frames in the trajectory.
// (Multiple consecutive blank lines are ignored.) 
//
// Any unexpected text causes the the function to gracefully exit (hopefully).

void
ReadCoordsRAW(istream& in,
              vector< vector<Real> >& vCoords,
              long long& line_count,
              char comment_char)
{
  vCoords.clear();

  char c = ' ';
  // First, skip over whitespace (including newlines) and comments
  while (in && (isspace(c) || c == comment_char)) {
    while (in.get(c) && isspace(c)) {
      //if (c=='\n') line_count++; //keep track of which line we're on
    }
    // Skip over comments
    if (! in) return;
    if (c == comment_char) {
      while (in.get(c) && (c != '\n')) {}
      if (c == '\n')
        line_count++; //keep track of which line we're on
    }
  }

  // Read one line at a time.
  //
  // We typically use this fuction to read a single conformation from a 
  // trajectory files which contain many conformations of the molecule,
  // separated by blank lines and/or comments.
  // Consequently, blank lines, lines containing only comments,
  // or lines containing non-numeric, non-commented text signal
  // that we have finished reading the shape of this molecule.

  while (in && 
         (isdigit(c) ||
          (c == '.') ||
          (c == '-') ||
          (c == '+'))) // <- read in one line at a time.
  {
    // Each line contains a list of numbers.

    vector<Real> position; //a list of coordinates for an atom/particle
    
    while (in &&
           (isdigit(c) ||
            (c == '.') ||
            (c == '-') ||
            (c == '+'))) // <- read in one number at a time.
    {

      string x_str;
      x_str.push_back(c);

      while (in.get(c) &&
             (isdigit(c) ||
              (c == '.') ||
              (c == '-') ||
              (c == '+') ||
              (c == 'e') ||
              (c == 'E')))
        x_str.push_back(c);

      // We have read one number.  Convert it, and store it in the list.
      Real x = atof(x_str.c_str());
      position.push_back(x);
      in.putback(c);

      // Now skip past any trailing spaces (excluding newlines) and comments 
      while (in.get(c) && isspace(c) && (c != '\n')) {}

    } //while (...) <- read in one number at a time.


    // Skip over any comments at the end of this line.
    if (in && (c == comment_char))
      while (in.get(c) && (c != '\n')) {}

    if (in) {
      if (c != '\n') 
        // Something terminated numeric entry
        // (probably non-numeric text).  In that case, put it back.
        in.putback(c); 
      else
        line_count++; //by not putting it back, we skip over it.
    }
    // Note: we only skip over at most 1 newline character in this loop.


    // Check for weird input.

    //if (vCoords.size() == 0)
    //if (g_dim == UNINITIALIZED)
    //  g_dim = position.size();
    //else if (g_dim != position.size()) {
    //  stringstream errmsg;
    //  errmsg << "Error: on line "<<line_count<<":\n"
    //    "      The number of coordinates on each line is not consistent.\n";
    //  throw InputErr(errmsg.str());
    //}

    if (g_dim != position.size()) {
      stringstream errmsg;
      errmsg << "Error: Wrong number of coordinates ("<<position.size()<<") on line "<<line_count<<":\n"
        "      This program expects " << g_dim << " coordinates on each line.\n";
      throw InputErr(errmsg.str());
    }

    vCoords.push_back(position);

    // Finally, skip past any leading spaces on the next line
    while (in.get(c) && isspace(c) && (c != '\n')) {}

  } //while (...) <- read in one line at a time.



  // Put back whichever non-numeric character caused termination.
  // It could be important (especially if it's not a newline).
  in.putback(c);

} //ReadCoordsRAW()




// This version saves the coordinate data in a pointer to a fixed-sized arrays
// (C style).  This assumes that each atom has 3 coordinates (3-dimensions).
// It returns the number of lines of numeric data read (the number of atoms).
long
ReadCoordsRAW(istream& in,
              Vect **paCoords, 
              long long& line_count,
              char comment_char)
{
  assert( paCoords != NULL);
  assert(*paCoords == NULL);

  // First load everything into a 2-D vector coordinates 
  // (these vectors are resizeable).
  vector< vector<Real> > vCoords;

  ReadCoordsRAW(in, vCoords, line_count, comment_char);

  // Then, copy the contents of vCoords[] into (*paCoords)[]
  Vect 
    *aCoords = new Vect [vCoords.size()];
  assert(aCoords != NULL);
  for (int r=0; r < vCoords.size(); ++r)
    for (int d=0; d < g_dim; ++d)
      aCoords[r][d] = vCoords[r][d];
  *paCoords = aCoords;
  return vCoords.size();

} //ReadCoordsRAW()






void
ReadCoordsXYZ(istream& in,
              vector< vector<Real> >& vCoords,
              long long& line_count,
              char comment_char)
{
  vCoords.clear();

    
  // First, skip over whitespace (including newlines) and comments
  SkipWhiteSpace(in, line_count, comment_char);

  // The first non-comment line should contain one number: the number of atoms
  // in the current frame.
  string s;
  getline(in, s);
  line_count++;
  long natoms = atol(s.c_str());

  // Skip over the next line (a mandatory comment-line).
  getline(in, s);
  line_count++;

  char c;

  // Read one line at a time.
  //
  // We typically use this fuction to read a single conformation from a 
  // trajectory files which contain many conformations of the molecule,
  // separated by blank lines and/or comments.
  // Consequently, blank lines, lines containing only comments,
  // or lines containing non-numeric, non-commented text signal
  // that we have finished reading the shape of this molecule.
  long atom_count = 0;
  while (in && (atom_count < natoms))
  {
    // Each line contains a list of numbers.

    vector<Real> position; //a list of coordinates for an atom/particle

    // Read in (and ignore) the atom-name string
    while (in.get(c) && isspace(c) && (c != '\n')) {}
    while (in.get(c) && (! isspace(c))) {}
    while (in.get(c) && (isspace(c) || (c==',')) && (c != '\n')) {}


    long ncrd = 0;
    while (in &&
           (ncrd < 3) &&
           (isdigit(c) ||
            (c == '.') ||
            (c == '-') ||
            (c == '+'))) // <- read in one number at a time.
    {

      string x_str;
      x_str.push_back(c);

      while (in.get(c) &&
             (isdigit(c) ||
              (c == '.') ||
              (c == '-') ||
              (c == '+') ||
              (c == 'e') ||
              (c == 'E')))
        x_str.push_back(c);
      in.putback(c);

      // We have read one number.  Convert it, and store it in the list.
      Real x = atof(x_str.c_str());
      position.push_back(x);

      // Now skip past any trailing spaces (excluding newlines) and comments 
      while (in.get(c) && (isspace(c) || (c==',')) && (c != '\n')) {}
      if (c=='\n')
        in.putback(c);

      ncrd++;
    } //while (...) <- read in one number at a time.

    // Skip over any remaining text on this line
    getline(in, s);
    line_count++;

    vCoords.push_back(position);

    atom_count++;
  } //while (...) <- read in one line at a time.

  // Skip over any trailing blank text before the next frame:
  SkipWhiteSpace(in, line_count, comment_char);

} //ReadCoordsXYZ()





// This version saves the coordinate data in a pointer to a fixed-sized arrays
// (C style).  This assumes that each atom has 3 coordinates (3-dimensions).
// It returns the number of lines of numeric data read (the number of atoms).
long
ReadCoordsXYZ(istream& in,
              Vect **paCoords, 
              long long& line_count,
              char comment_char)
{
  assert( paCoords != NULL);
  assert(*paCoords == NULL);

  // First load everything into a 2-D vector coordinates 
  // (these vectors are resizeable).
  vector< vector<Real> > vCoords;

  ReadCoordsXYZ(in, vCoords, line_count, comment_char);

  // Then, copy the contents of vCoords[] into (*paCoords)[]
  Vect
    *aCoords = new Vect [vCoords.size()];
  assert(aCoords != NULL);
  for (int r=0; r < vCoords.size(); ++r)
    for (int d=0; d < g_dim; ++d)
      aCoords[r][d] = vCoords[r][d];
  *paCoords = aCoords;
  return vCoords.size();

} //ReadCoordsXYZ()



