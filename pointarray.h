// Copyright 2012 Esteban Hurtado

// This file is part of Cutedots.

// Cutedots is distributed under the terms of the Reciprocal Public License 1.5.

// Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.

// You should have received a copy of the Reciprocal Public License along with
// Cutedots. If not, see <http://opensource.org/licenses/rpl-1.5>.

#include <vector>
using namespace std;

class PointArray
{
public:
    PointList();
    void addPoint(const float* coords);
    float* getPoint(size_t index);
private:
    vector<float> _data;
};
