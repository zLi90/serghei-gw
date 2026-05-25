#ifndef _GEOMETRY_
#define _GEOMETRY_

#include "SArray.h"

#define _X 0
#define _Y 1
#define _Z 2

#define _SOUTH 0
#define _EAST 1
#define _NORTH 2
#define _WEST 3
#define N_SPATIAL_DIM 2
#define N_CELL_EDGES 4

namespace geometry
{
	typedef SArray<real, N_SPATIAL_DIM> point;

	typedef SArray<real, 3> point3D; // ParticleTracking // [FROM CODE2] needed for particle tracking

	KOKKOS_INLINE_FUNCTION real distance(const point &p, const point &q)
	{
		point d;
		real ss = 0.;
		for (int i = 0; i < N_SPATIAL_DIM; i++)
		{
			d(i) = p(i) - q(i);
			ss += d(i) * d(i);
		}
		return (mysqrt(ss));
	};

	inline point crossProduct(const point &a, const point &b)
	{
		point c;
		c(0) = a(1) * b(2) - a(2) * b(1);
		c(1) = a(2) * b(0) - a(0) * b(2);
		c(2) = a(0) * b(1) - a(1) * b(0);
		return (c);
	};

	inline geometry::point normalToLineInPlane(const point &a, const point &b)
	{
		point c;
		c(0) = a(1) - b(1);
		c(1) = b(0) - a(0);
		return (c);
	};

	// returns 1 if the coordinate is inside, 0 otherwise
	int isInsidePoly(int np, realArr &xPoly, realArr &yPoly, real &xCoord, real &yCoord)
	{
		int count;
		int i;
		double xInter;
		real x1, x2, y1, y2;

		count = 0;
		for (i = 0; i < np; i++)
		{
			x1 = xPoly(i);
			y1 = yPoly(i);
			x2 = xPoly((i + 1) % np);
			y2 = yPoly((i + 1) % np);
			if ((yCoord > fmin(y1, y2)) && (yCoord <= fmax(y1, y2)) && (xCoord <= fmax(x1, x2)))
			{
				if (y1 != y2)
				{
					xInter = x1 + (yCoord - y1) * (x2 - x1) / (y2 - y1);
					if (x1 == x2 || xCoord <= xInter)
					{
						count++;
					}
				}
			}
		}

		if (count % 2 == 0)
		{
			return 0;
		}
		else
		{
			return 1;
		}
	}

	// [FROM CODE1 - MERGED] Improved isInsidePoly3D with inverse distance weighting for z_bottom
	// This version uses a proper ray-casting algorithm for XY projection,
	// then uses inverse distance weighting (IDW) to compute the bottom elevation z_bottom,
	// instead of the simple min/max comparison in the original code2 version.
	// This is critical for the source/sink (root water uptake) feature where the 3D polygon
	// defines the domain of subsurface cells affected by the source/sink term.
	int isInsidePoly3D(int np, realArr &xPoly, realArr &yPoly, realArr &zPoly, real &xCoord, real &yCoord, real &zCoord)
	{
		int count;
		int i;
		double xInter;
		real x1, x2, y1, y2;

		count = 0;
		// 1. Determine if inside the 2D polygon using ray-casting
		for (i = 0; i < np; i++)
		{
			x1 = xPoly(i);
			y1 = yPoly(i);
			x2 = xPoly((i + 1) % np);
			y2 = yPoly((i + 1) % np);
			if ((yCoord > fmin(y1, y2)) && (yCoord <= fmax(y1, y2)) && (xCoord <= fmax(x1, x2)))
			{
				if (y1 != y2)
				{
					xInter = x1 + (yCoord - y1) * (x2 - x1) / (y2 - y1);
					if (x1 == x2 || xCoord <= xInter)
					{
						count++;
					}
				}
			}
		}

		// If count is even, the point is outside the polygon
		if (count % 2 == 0)
		{
			return 0;
		}
		else
		{
			// 2. Inside the XY polygon - compute the bottom elevation (z_bottom) 
			// using Inverse Distance Weighting (IDW) based on polygon vertices
			double sum_z = 0.0;
			double sum_w = 0.0;
			double z_bottom = 0.0;
			int exact_match = 0;

			for (i = 0; i < np; i++)
			{
				double dx = xCoord - xPoly(i);
				double dy = yCoord - yPoly(i);
				double dist2 = dx * dx + dy * dy;

				// If the grid center coincides with a polygon vertex
				if (dist2 < 1e-12)
				{
					z_bottom = zPoly(i);
					exact_match = 1;
					break;
				}

				double w = 1.0 / dist2;
				sum_w += w;
				sum_z += zPoly(i) * w;
			}

			if (exact_match == 0)
			{
				z_bottom = sum_z / sum_w;
			}

			// 3. Check if the grid Z coordinate is within range
			if (zCoord >= z_bottom)
			{
				return 1; // Within the 3D domain controlled by source/sink
			}
			else
			{
				return 0; // Below the bottom surface
			}
		}
	}

};
#endif