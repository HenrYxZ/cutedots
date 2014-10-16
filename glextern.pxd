# Copyright 2012 Esteban Hurtado
#
# This file is part of Cutedots.
#
# Cutedots is distributed under the terms of the Reciprocal Public License 1.5.
#
# Cutedots is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the Reciprocal Public License 1.5 for more details.
#
# You should have received a copy of the Reciprocal Public License along with
# Cutedots. If not, see <http://http://opensource.org/licenses/rpl-1.5>.

cdef extern from "GL/gl.h":
     ctypedef int GLint
     ctypedef unsigned int GLenum
     ctypedef float GLfloat
     ctypedef double GLdouble
     ctypedef float GLclampf
     ctypedef unsigned int GLbitfield
     ctypedef unsigned char GLboolean
     ctypedef int GLsizei

     enum:
         GL_MODELVIEW = 0x1700
         GL_LIGHTING = 0x0B50
         GL_LIGHT0 = 0x4000
         GL_QUADS = 0x0007
         GL_LINES = 0x0001
         GL_LINE_LOOP = 0x0002
         GL_LINE_STRIP = 0x0003
         GL_QUAD_STRIP = 0x0008
         GL_POLYGON = 0x0009
         GL_FALSE = 0x0
         GL_TRUE = 0x1
         GL_DEPTH_BUFFER_BIT = 0x00000100
         GL_COLOR_BUFFER_BIT = 0x00004000
         GL_PROJECTION = 0x1701
         GL_DEPTH_TEST = 0x0B71
         GL_SMOOTH = 0x1D01
         GL_COLOR_MATERIAL = 0x0B57
         GL_POSITION = 0x1203
         GL_AMBIENT = 0x1200
         GL_DIFFUSE = 0x1201
         GL_SPECULAR = 0x1202
         GL_SHININESS = 0x1601
         GL_LINEAR_ATTENUATION = 0x1208
         GL_FRONT_AND_BACK = 0x0408
         GL_BLEND = 0x0BE2
         GL_SRC_ALPHA = 0x0302
         GL_ONE_MINUS_SRC_ALPHA = 0x0303

     cdef void glMatrixMode( GLenum mode )
     cdef void glLoadIdentity()
     cdef void glRotatef( GLfloat angle, GLfloat x, GLfloat y, GLfloat z )
     cdef void glPushMatrix()
     cdef void glPopMatrix()
     cdef void glEnable( GLenum cap )
     cdef void glDisable( GLenum cap )
     cdef void glColor3f( GLfloat red, GLfloat green, GLfloat blue )
     cdef void glColor4f( GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha )
     cdef void glRasterPos2f( GLfloat x, GLfloat y )
     cdef void glRasterPos3f( GLfloat x, GLfloat y, GLfloat z )
     cdef void glTranslatef( GLfloat x, GLfloat y, GLfloat z )
     cdef void glLineWidth( GLfloat width )
     cdef void glBegin( GLenum mode )
     cdef void glEnd()
     cdef void glVertex3f( GLfloat x, GLfloat y, GLfloat z )
     cdef void glVertex2f( GLfloat x, GLfloat y )
     cdef void glColorMask( GLboolean red, GLboolean green, GLboolean blue, GLboolean alpha )
     cdef void glClear( GLbitfield mask )
     cdef void glFlush()
     cdef void glShadeModel( GLenum mode )
     cdef void glClearColor( GLclampf red, GLclampf green, GLclampf blue, GLclampf alpha )
     cdef void glLightf( GLenum light, GLenum pname, GLfloat param )
     cdef void glLightfv( GLenum light, GLenum pname, GLfloat *params )
     cdef void glMaterialf( GLenum face, GLenum pname, GLfloat param )
     cdef void glMaterialfv( GLenum face, GLenum pname, GLfloat *params )
     cdef void glBlendFunc( GLenum sfactor, GLenum dfactor )
     cdef void glViewport( GLint x, GLint y, GLsizei width, GLsizei height )
     cdef void glOrtho( GLdouble left, GLdouble right, \
                        GLdouble bottom, GLdouble top, \
                        GLdouble near_val, GLdouble far_val )

cdef extern from "GL/glu.h":
     ctypedef double GLdouble
     cdef void gluLookAt (GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ, GLdouble centerX, \
     	  GLdouble centerY, GLdouble centerZ, GLdouble upX, GLdouble upY, GLdouble upZ)
     cdef void gluPerspective (GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)


