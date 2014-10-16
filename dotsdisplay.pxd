cdef extern from "dots_drawing.h":
    cdef void setup_drawing(int width, int height)

    cdef void drawSolidSphere(double radius, int slices, int stacks)
    cdef void drawSolidCube(double size)

    cdef void draw_wall(float width, float height, int hdivs, int vdivs)
    cdef void draw_side_wall(float tx, float ty, float tz, float rx, float ry, float rz, \
                             float w, float h, int hdivs, int vdivs)
    cdef void drawFloor()
    cdef void drawWalls()
    cdef void drawChairs()

    cdef void draw_single_head(float x1, float y1, float z1, \
		               float x2, float y2, float z2, \
		               float x3, float y3, float z3)

    cdef void draw_torso(float xa, float ya, float za, float xb, float yb, float zb, \
                         float xc, float yc, float zc, float xd, float yd, float zd)


cdef class DotsDisplay:
    cdef public object model
    cdef public object textRenderFunc
    cdef void drawText(DotsDisplay self, str s)
    cdef void subjColor(DotsDisplay self, str subj, float alpha)
    cdef void drawElement(DotsDisplay self, traj)
    cdef void drawElements(DotsDisplay self, int frame)
    cdef void drawHeads(DotsDisplay self, int frame)
    cdef void drawSelSquare(DotsDisplay self, int frame)
    cdef void drawLine(DotsDisplay self, str tr1, str tr2)
    cdef void drawLines(DotsDisplay self)
    cdef void locateCamera(DotsDisplay self, int cam)
    cdef void rotateScene(DotsDisplay self)
    cdef void draw(DotsDisplay self)
