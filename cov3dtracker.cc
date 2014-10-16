// Copyright 2012 Esteban Hurtado

// This file is part of Cutedots.

// Cutedots is distributed under the terms of the Reciprocal Public License 1.5.

// Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.

// You should have received a copy of the Reciprocal Public License along with
// Cutedots. If not, see <http://opensource.org/licenses/rpl-1.5>.


#include "errorbuffer.h"

Cov3DTracker::Cov3DTracker(size_t bufferSize, Real initVar)
    : _x(bufferSize), _y(bufferSize), _z(bufferSize),
      _sx(initVar), _sy(initVar), _sz(initVar),
      _sxy(Real(0.0)), _sxz(Real(0.0)), _syz(Real(0.0))
{}

void
Cov3DTracker::addError(Real x, Real y, Real z)
{
    _x.append(x); _y.append(y); _z.append(z);
    computeCovariance();
}


void
Cov3DTracker::computeCovariance()
{
    const size_t bufferSize = _x.size();

    if (bufferSize < 2)
	return;

    // compute means
    // TODO Would be faster to keep track of buffer means. I.e., keep a sum
    // variable, updating each time a value is appended by adding new value
    // and substracting the one that goes out. Can be done accurately in
    // floating-point arithmetic?

    Real xm=Real(0.0), ym=Real(0.0), zm=Real(0.0);

    for (size_t i=0; i<bufferSize; i++) {
	xm += _x[i];
	ym += _y[i];
	ym += _z[i];
    }

    xm /= bufferSize;
    ym /= bufferSize;
    zm /= bufferSize;

    // sum squares
    // TODO Ok, not the best algorithm, but reasonably stable.

    _sx = _sy = _sz = _sxy = _sxz = _syz = Real(0.0);

    for (size_t i=0; i<bufferSize; i++) {
	float dx = _sx - xm, dy  = _sy - ym, dz = _sz - zm;
	_sx += dx*dx;
	_sy += dy*dy;
	_sz += dz*dz;
	_sxy += dx*dy;
	_sxz += dx*dz;
	_syz += dy*dz;
    }

    size_t denom = bufferSize - 1;
    _sx  /= denom;
    _sy  /= denom;
    _sz  /= denom;
    _sxy /= denom;
    _sxz /= denom;
    _syz /= denom;
}
