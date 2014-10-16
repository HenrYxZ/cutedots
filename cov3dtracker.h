// Copyright 2012 Esteban Hurtado

// This file is part of Cutedots.

// Cutedots is distributed under the terms of the Reciprocal Public License 1.5.

// Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.

// You should have received a copy of the Reciprocal Public License along with
// Cutedots. If not, see <http://opensource.org/licenses/rpl-1.5>.

/// @file
/// Defines class Cov3DTracker for computing covariance matrix for a list
/// of recently given points.

#include "ringbuffer.h"

typedef float Real;
typedef RingBuffer<Real> Buffer;

/// Cov3DTracker keeps track of recently given 3D points, and
/// estimates a covariance matrix for them.
class Cov3DTracker
{
public:

    /// Constructor
    /// @param bufferSize Indicates how many recent 3D points to keep track of.
    /// @param initVar Initial value for variances (covariances start at zero).
    ///        If you need to inialize covariances, assign directly after construction.
    Cov3DTracker(size_t bufferSize, Real initVar);

    /// Adds a new point to tracking and updates covariance matrix.
    /// If buffer already full, oldest 3D point is lost.
    /// @param x X value.
    /// @param y Y value.
    /// @param z Z value.
    void addPoint(Real x, Real y, Real z);

    Real sx;    ///< X variance
    Real sy;    ///< Y variance
    Real sz;    ///< Z variance
    Real sxy;   ///< XY covariance
    Real sxz;   ///< XZ covariance
    Real syz;   ///< YZ covariance

protected:
    void computeCovariance();

private:
    Buffer _x;  // Storage for recent points.
    Buffer _y;
    Buffer _z;
};
