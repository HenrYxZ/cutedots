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

#include <iostream>
#include <set>
#include <cmath>
#include <iomanip>
#include "heads.h"

using namespace std;

const float EPS3D = 10.0;

float head_pattern[6];


/// Indices for an edge
///////////////////////

struct edge
{
  int i, j;         ///< indices
  float length;     ///< edge length

  edge(int i_, int j_, float length_) :
    i(i_), j(j_), length(length_) {}
};


/// Edge comparison functor
struct ltedge
{
  bool operator()(const edge &e1, const edge &e2) const
  { return e1.i < e1.j; }
};

typedef set<edge, ltedge> edge_set;
typedef vector<edge> edge_vec;



/// Distance computation
////////////////////////

inline float*
get_pt(float* data, int i)
{
  return data + (3*i);
}

inline float
calc_sq_dist(float* data, int i, int j)
{
  float *a = get_pt(data, i);
  float *b = get_pt(data, j);
  float dx = b[0] - a[0];
  float dy = b[1] - a[1];
  float dz = b[2] - a[2];
  return dx*dx + dy*dy + dz*dz;
}


/// Trajectory
//////////////

class trajectory
{
public:
  trajectory(size_t first_frame) :
    _first_frame(first_frame), _data() {}

  void add(float x, float y, float z)
  { _data.push_back(x);
    _data.push_back(y);
    _data.push_back(z); }

  const vector<float>& get_data() const
  { return _data; }

  size_t size() const
  { return _data.size() / 3; }

private:
  size_t _first_frame;
  vector<float> _data;
};


/// Head center
///////////////

void
heads::head_center(int hd, int idx, float* data, float* center) const
{
  const triangle& head = _heads[hd][idx];

  float* p1 = &data[3*head.i];
  float* p2 = &data[3*head.j];
  float* p3 = &data[3*head.k];

  center[0] = (p1[0] + p2[0] + p3[0]) / 3.0;
  center[1] = (p1[1] + p2[1] + p3[1]) / 3.0;
  center[2] = (p1[2] + p2[2] + p3[2]) / 3.0;
}

/// Find heads
//////////////

void find_heads(float *data, int numPoints, heads* head_info)
{
  vector<edge_vec> edges(6);

  // find edges
  for (size_t i=0; i<numPoints; i++) {
    for (size_t j=i+1; j<numPoints; j++) {

      for (size_t k=0; k<6; k++) {
	float diff = fabs(sqrt(head_pattern[k]) - sqrt(calc_sq_dist(data,i,j)));

	if ( diff < EPS3D )
	  edges[k].push_back( edge(i, j, diff) );
      }
    }
  }


  head_info->clear();

  // for the two heads
  for (size_t hd=0; hd<2; hd++) {
    int edge1 = 3*hd;
    int edge2 = edge1 + 1;
    int edge3 = edge1 + 2;

    // edge 1
    for (size_t e1=0; e1<edges[edge1].size(); e1++) {
      int a1 = edges[edge1][e1].i;
      int b1 = edges[edge1][e1].j;

      // edge 2
      for (size_t e2=0; e2<edges[edge2].size(); e2++) {
	int a2 = edges[edge2][e2].i;
	int b2 = edges[edge2][e2].j;

	if (a1 == a2) {
	  // edge 3
	  for (size_t e3=0; e3<edges[edge3].size(); e3++) {
	    int a3 = edges[edge3][e3].i;
	    int b3 = edges[edge3][e3].j;

	    if ( (b1 == a3 && b2 == b3) ||
		 (b1 == b3 && b2 == a3) )
	      head_info->add_head(hd, a1, a3, b3);
	  }
	}

	else if (a1 == b2) {
	  // edge 3
	  for (size_t e3=0; e3<edges[edge3].size(); e3++) {
	    int a3 = edges[edge3][e3].i;
	    int b3 = edges[edge3][e3].j;

	    if ( (b1 == a3 && a2 == b3) ||
		 (b1 == b3 && a2 == a3) )
	      head_info->add_head(hd, a1, a3, b3);
	  }
	}

	else if (b1 == a2) {
	  // edge 3
	  for (size_t e3=0; e3<edges[edge3].size(); e3++) {
	    int a3 = edges[edge3][e3].i;
	    int b3 = edges[edge3][e3].j;

	    if ( (a1 == a3 && b2 == b3) ||
		 (a1 == b3 && b2 == a3) )
	      head_info->add_head(hd, b1, a3, b3);
	  }
	}

	else if (b1 == b2) {
	  // edge 3
	  for (size_t e3=0; e3<edges[edge3].size(); e3++) {
	    int a3 = edges[edge3][e3].i;
	    int b3 = edges[edge3][e3].j;

	    if ( (a1 == a3 && a2 == b3) ||
		 (a1 == b3 && a2 == a3) ) {
	      head_info->add_head(hd, b1, a3, b3);
	    }
	  }
	}
      }
    }
  }
}


