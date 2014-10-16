// Copyright 2012 Esteban Hurtado

// This file is part of Cutedots.

// Cutedots is distributed under the terms of the Reciprocal Public License 1.5.

// Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.

// You should have received a copy of the Reciprocal Public License along with
// Cutedots. If not, see <http://opensource.org/licenses/rpl-1.5>.

#includde"pointlist.h"
#include <cmath>

using namespace std;


PointArray::PointArray() :
    _data(),
{}


void
PointArray::addPoint(const float* coords)
{
    // add point
    _data.push_back(coords[0]);
    _data.push_back(coords[1]);
    _data.push_back(coords[2]);
}


float*
PointArray::getPoint(size_t index)
{
    return &( _data[index*3] );
}
