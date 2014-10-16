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

#ifndef _DOTS_DRAWING_H
#define _DOTS_DRAWING_H

void setup_drawing(int width, int height);

void draw_wall(float width, float height, int hdivs, int vdivs);

void draw_side_wall(float tx, float ty, float tz, float rx, float ry, float rz,
		    float w, float h, int hdivs, int vdivs);

void drawFloor(void);

void drawWalls(void);

void drawChairs(void);

void draw_single_head(float x1, float y1, float z1,
		      float x2, float y2, float z2,
		      float x3, float y3, float z3);

void draw_torso(float xa, float ya, float za, float xb, float yb, float zb,
		float xc, float yc, float zc, float xd, float yd, float zd);


#endif

