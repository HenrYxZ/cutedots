// Copyright 2012 Esteban Hurtado
//
// This file is part of Cutedots.
//
// Cutedots is distributed under the terms of the Reciprocal Public License 1.5.
//
// Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.
//
// You should have received a copy of the Reciprocal Public License along with
// Cutedots. If not, see <http://http://opensource.org/licenses/rpl-1.5>.

#ifndef _HEADS_H
#define _HEADS_H

#include <vector>

using namespace std;

/// Indices for three points
////////////////////////////

struct triangle
{
  int i, j, k;       ///< indices
  //  float e1, e2, e3;  ///< edge lengths (ij, jk, ki)

  triangle(int i_, int j_, int k_) : //, float e1_, float e2_, float e3_) :
    i(i_), j(j_), k(k_) {} //, e1(e1_), e2(e2_), e3(e3_) {}
};


/// Possible heads information
//////////////////////////////

class heads
{
 public:
  heads() : _heads(2) {}

  void add_head(int hd, int i, int j, int k)
  { _heads[hd].push_back( triangle(i,j,k) ); }

  const triangle& get_head(int hd, int idx) const
  { return _heads[hd][idx]; }

  const size_t get_head_len(int hd) const
  { return _heads[hd].size(); }

  const int get_i(int hd, int idx) const
  { return _heads[hd][idx].i; }

  const int get_j(int hd, int idx) const
  { return _heads[hd][idx].j; }

  const int get_k(int hd, int idx) const
  { return _heads[hd][idx].k; }

  void head_center(int hd, int idx, float* data, float* center) const;

  void clear()
  { _heads[0].clear(); _heads[1].clear(); }

 private:
  vector< vector<triangle> > _heads;
};


void find_heads(float *data, int numPoints, heads* head_info);
extern float head_pattern[6];


#endif // _HEADS_H
