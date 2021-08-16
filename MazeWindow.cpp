/************************************************************************
     File:        MazeWindow.cpp

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for the MazeWindow class. The MazeWindow is
						the window in which the viewer's view of the maze is displayed.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "MazeWindow.h"
#include <Fl/math.h>
#include <Fl/gl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <iostream>

float* ProjectionMatrix;
float* ViewMatrix;

//*************************************************************************
//
// * Constructor
//=========================================================================
MazeWindow::
MazeWindow(int x, int y, int width, int height, const char *label,Maze *m)
	: Fl_Gl_Window(x, y, width, height, label)
//=========================================================================
{
	maze = m;

	// The mouse button isn't down and there is no key pressed.
	down = false;
	z_key = 0;
}


//*************************************************************************
//
// * Set the maze. Also causes a redraw.
//=========================================================================
void MazeWindow::
Set_Maze(Maze *m)
//=========================================================================
{
	// Change the maze
	maze = m;

	// Force a redraw
	redraw();
}


//*************************************************************************
//
// * draw() method invoked whenever the view changes or the window
//   otherwise needs to be redrawn.
//=========================================================================
void MazeWindow::
draw(void)
//=========================================================================
{
	float   focal_length;

	if ( ! valid() ) {
		// The OpenGL context may have been changed
		// Set up the viewport to fill the window.
		glViewport(0, 0, w(), h());

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		// We are using orthogonal viewing for 2D. This puts 0,0 in the
		// middle of the screen, and makes the image size in view space
		// the same size as the window.
		gluOrtho2D(-w() * 0.5, w() * 0.5, -h() * 0.5, h() * 0.5);

		// Sets the clear color to black.
		glClearColor(0.0, 0.0, 0.0, 1.0);
	}

	// Clear the screen.
	glClear(GL_COLOR_BUFFER_BIT);

	glBegin(GL_QUADS);
		// Draw the "floor". It is an infinite plane perpendicular to
		// vertical, so we know it projects to cover the entire bottom
		// half of the screen. Walls of the maze will be drawn over the top
		// of it.
		glColor3f(0.2f, 0.2f, 0.2f);
		glVertex2f(-w() * 0.5f, -h() * 0.5f);
		glVertex2f( w() * 0.5f, -h() * 0.5f);
		glVertex2f( w() * 0.5f, 0.0       );
		glVertex2f(-w() * 0.5f, 0.0       );

		// Draw the ceiling. It will project to the entire top half
		// of the window.
		glColor3f(0.4f, 0.4f, 0.4f);
		glVertex2f( w() * 0.5f,  h() * 0.5f);
		glVertex2f(-w() * 0.5f,  h() * 0.5f);
		glVertex2f(-w() * 0.5f, 0.0       );
		glVertex2f( w() * 0.5f, 0.0       );
	glEnd();


	if ( maze ) {
		// Set the focal length. We can do this because we know the
		// field of view and the size of the image in view space. Note
		// the static member function of the Maze class for converting
		// radians to degrees. There is also one defined for going backwards.
		focal_length = w()
						 / (float)(2.0*tan(Maze::To_Radians(maze->viewer_fov)*0.5));

		glClear(GL_DEPTH_BUFFER_BIT);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		float aspect = (float)w() / h();
		/*gluPerspective(maze->viewer_fov, aspect, 0.01, 200);*/

		ProjectionMatrix = Perspective(maze->viewer_fov, aspect, 0.01, 200);

		//glLoadMatrixf(ProjectionMatrix);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		float viewer_pos[3] = { maze->viewer_posn[Maze::Y], 0.0f, maze->viewer_posn[Maze::X] };
		/*gluLookAt(viewer_pos[Maze::X], viewer_pos[Maze::Y], viewer_pos[Maze::Z],
			viewer_pos[Maze::X] + sin(Maze::To_Radians(maze->viewer_dir)),
			viewer_pos[Maze::Y],
			viewer_pos[Maze::Z] + cos(Maze::To_Radians(maze->viewer_dir)),
			0.0, 1.0, 0.0);*/

		
		ViewMatrix = LookAt(viewer_pos[Maze::X], viewer_pos[Maze::Y], viewer_pos[Maze::Z],
							 viewer_pos[Maze::X] + sin(Maze::To_Radians(maze->viewer_dir)), viewer_pos[Maze::Y], viewer_pos[Maze::Z] + cos(Maze::To_Radians(maze->viewer_dir)),
							 0.0, 1.0, 0.0);

		//glLoadMatrixf(ViewMatrix);
	
		// Draw the 3D view of the maze (the visible walls.) You write this.
		// Note that all the information that is required to do the
		// transformations and projection is contained in the Maze class,
		// plus the focal length.
		maze->Draw_View(focal_length);
	}
}

float* MazeWindow::
Perspective(float fovy, float aspect, float zNear, float zFar)
{
	float* perspectiveMatrix = new float[16]{0};
	float f = (float)1.0 / tanf(Maze::To_Radians(fovy / 2));

	float top = (float)zNear * tanf(Maze::To_Radians(fovy / 2));
	float buttom = (float)-top;
	float right = (float)top * aspect;
	float left = (float)-right;

	perspectiveMatrix[0] = (float)2 * zNear / (right - left);
	perspectiveMatrix[2] = (float)(right + left) / (right - left);
	perspectiveMatrix[5] = (float)2 * zNear / (top - buttom);
	perspectiveMatrix[6] = (float)(top + buttom) / (top - buttom);
	perspectiveMatrix[10] = (float)(zNear + zFar) / (zNear - zFar);
	perspectiveMatrix[11] = (float)2 * zNear * zFar / (zNear - zFar);
	perspectiveMatrix[14] = (float)-1.0;

	perspectiveMatrix = TransformationMatrix(perspectiveMatrix);

	return perspectiveMatrix;
}

float* MazeWindow::
LookAt(float eyeX, float eyeY, float eyeZ, float centerX, float centerY, float centerZ, float upX, float upY, float upZ)
{
	float cameraPosition[] = { eyeX, eyeY, eyeZ };
	float target[] = { centerX, centerY, centerZ };
	float up[] = { upX, upY, upZ };


	float* viewMatrix;
	float* translateMatrix = new float[16]{ 1, 0, 0, (float)-eyeX,
											0, 1, 0, (float)-eyeY,
											0, 0, 1, (float)-eyeZ,
											0, 0, 0, 1 };

	float* zAxis = Normalize(SubtractVectors(cameraPosition, target));
	float* xAxis = Normalize(Cross(up, zAxis));
	float* yAxis = Normalize(Cross(zAxis, xAxis));

	float* rotateMatrix = new float[16] { xAxis[0], xAxis[1], xAxis[2], 0,
										  yAxis[0], yAxis[1], yAxis[2], 0,
										  zAxis[0], zAxis[1], zAxis[2], 0,
										  0,		   0,		 0,		   1 };

	viewMatrix = MultMatrix(rotateMatrix, translateMatrix);

	viewMatrix = TransformationMatrix(viewMatrix);

	return viewMatrix;
}

float* MazeWindow::
SubtractVectors(float* a, float* b)
{
	float* vector = new float [3]{ 0 };

	vector[0] = a[0] - b[0];
	vector[1] = a[1] - b[1];
	vector[2] = a[2] - b[2];

	return vector;
}

//Make it into a unit vector
float* MazeWindow::
Normalize(float* v)
{
	float* vector = new float[3]{ 0 };

	float length = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	if (length > 0)
	{
		vector[0] = (float)v[0] / length;
		vector[1] = (float)v[1] / length;
		vector[2] = (float)v[2] / length;
	}

	return vector;
}

float* MazeWindow::
Cross(float* a, float* b)
{
	float* vector = new float [3] { a[1] * b[2] - a[2] * b[1],
									a[2] * b[0] - a[0] * b[2],
									a[0] * b[1] - a[1] * b[0] };

	return vector;
}

float* MazeWindow::
MultMatrix(float* a, float* b)
{
	float* matrix = new float[16]{ 0 };

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			for (int k = 0; k < 4; ++k)
				matrix[i * 4 + j] += a[i * 4 + k] * b[k * 4 + j];

	return matrix;
}

float* MazeWindow::
TransformationMatrix(float* matrix)
{
	float* m = new float[16]{ 0 };

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			m[j * 4 + i] = matrix[i * 4 + j];

	return m;
}

//*************************************************************************
//
// *
//=========================================================================
bool MazeWindow::
Drag(float dt)
//=========================================================================
{
	float   x_move, y_move, z_move;

	if ( down ) {
		int	dx = x_down - x_last;
		int   dy = y_down - y_last;
		float dist;

		// Set the viewing direction based on horizontal mouse motion.
		maze->Set_View_Dir(d_down + 360.0f * dx / (float)w());

		// Set the viewer's linear motion based on a speed (derived from
		// vertical mouse motion), the elapsed time and the viewing direction.
		dist = 10.0f * dt * dy / (float)h();
		x_move = dist * (float)cos(Maze::To_Radians(maze->viewer_dir));
		y_move = dist * (float)sin(Maze::To_Radians(maze->viewer_dir));
	}
	else {
		x_move = 0.0;
		y_move = 0.0;
	}

	// Update the z posn
	z_move = z_key * 0.1f;
	z_key = 0;

	// Tell the maze how much the view has moved. It may restrict the motion
	// if it tries to go through walls.
	maze->Move_View_Posn(x_move, y_move, z_move);

	return true;
}


//*************************************************************************
//
// *
//=========================================================================
bool MazeWindow::
Update(float dt)
//=========================================================================
{
	// Update the view

	if ( down || z_key ) // Only do anything if the mouse button is down.
		return Drag(dt);

	// Nothing changed, so no need for a redraw.
	return false;
}


//*************************************************************************
//
// *
//=========================================================================
int MazeWindow::
handle(int event)
//=========================================================================
{
	if (!maze)
		return Fl_Gl_Window::handle(event);

	// Event handling routine.
	switch ( event ) {
		case FL_PUSH:
			down = true;
			x_last = x_down = Fl::event_x();
			y_last = y_down = Fl::event_y();
			d_down = maze->viewer_dir;
			return 1;
		case FL_DRAG:
			x_last = Fl::event_x();
			y_last = Fl::event_y();
			return 1;
			case FL_RELEASE:
			down = false;
			return 1;
		case FL_KEYBOARD:
			/*
			if ( Fl::event_key() == FL_Up )	{
				z_key = 1;
				return 1;
			}
			if ( Fl::event_key() == FL_Down ){
				z_key = -1;
				return 1;
			}
			*/
			return Fl_Gl_Window::handle(event);
		case FL_FOCUS:
		case FL_UNFOCUS:
			return 1;
	}

	// Pass any other event types on the superclass.
	return Fl_Gl_Window::handle(event);
}


