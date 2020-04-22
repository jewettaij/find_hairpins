#ifndef _IO_H
#define _IO_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;
#include "global_defs.h"


void
SkipWhiteSpace(istream &in,
               long long& line_count,
               char comment_char='#');


void
ReadCoordsRAW(istream& in,
              vector< vector<Real> >& vCoords,
              long long& line_count,
              char comment_char='#');

long
ReadCoordsRAW(istream& in,
              Vect **paCoords, 
              long long& line_count,
              char comment_char='#');

void
ReadCoordsXYZ(istream& in,
              vector< vector<Real> >& vCoords,
              long long& line_count,
              char comment_char='#');

long
ReadCoordsXYZ(istream& in,
              Vect **paCoords, 
              long long& line_count,
              char comment_char='#');



#endif //#ifndef _IO_H
