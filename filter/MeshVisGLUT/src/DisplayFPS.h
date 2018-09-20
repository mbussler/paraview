#pragma once

#ifndef DISPLAY_FPS_H
#define DISPLAY_FPS_H

/**
 *  \brief FPS Counter with OpenGL output
 * 
 *  Code originally from lighthouse3D
 *
 */

namespace FPS
{

	// frame should be incremented each time a frame was drawn
	static int counter=0;
	static int iterations=0;
	static int totalParticles;
	static int activeParticles;
	static float timestep;
	static float hops[3];
	static size_t mem_free;
	static size_t mem_total;
	static bool record=false;

	static vector<double> fps_recorder;
	static vector<double> iter_recorder;

	// window parameter for the ortographic projection
	// !!! should be set as well !!!
	static int window_width=800;
	static int window_height=600;

	// time measurement
	static int ellapsed_time,timebase=0;

	// font output
	static void* font= GLUT_BITMAP_8_BY_13;
	static char sFps[50];
	static char sParticleCount[40];
	static char sFrames[40];
	static char sHops[60];
	static char sTimestep[20];
	static char sMem[50];

	static void displayFps();
	static void setOrthographicProjection();
	static void resetPerspectiveProjection();
	static void renderBitmapString(float x, float y, void *font,char *string);
	static void drawFramecount(int currentFrame, int maxFrame);
	static void drawWriteMessage(int currentFile, int maxFile);

	static void displayFps()
	{
		ellapsed_time=glutGet(GLUT_ELAPSED_TIME);

		if (ellapsed_time - timebase > 500.0f)
		{
			float fps = counter*500.0f/(ellapsed_time-timebase);
			float sps = iterations*500.0f/(ellapsed_time-timebase);

			sprintf(sFps,"%.2f frames/s, %.2f steps/s", fps, sps );
			sprintf(sParticleCount,"Particles: %d", totalParticles );
			sprintf(sTimestep, "Timestep: %.3f", timestep);
			sprintf(sHops, "Traversed Cells avg/max: %.3f/%.3f", hops[2], hops[1] );
			sprintf(sMem, "Memory free/total [MB]: %d/%d", mem_free, mem_total );

			timebase = ellapsed_time;
			counter = 0;
			iterations = 0;

			if( record )
			{
				fps_recorder.push_back(fps);
				iter_recorder.push_back(sps);
			}
		}

		glPushMatrix();
		glLoadIdentity();

		setOrthographicProjection();

		renderBitmapString(30, 110, font, sFps);
		renderBitmapString(30, 90, font, sTimestep );
		renderBitmapString(30, 73, font, sParticleCount);
		renderBitmapString(30, 56, font, sHops);
		renderBitmapString(30, 39, font, sMem);

		resetPerspectiveProjection();

		glPopMatrix();

	};

	static void drawFramecount(int currentFrame, int maxFrame)
	{
		sprintf(sFrames, "Frame %d", currentFrame );
		sprintf(sFrames, "%s / ", sFrames);
		sprintf(sFrames, "%s%d.", sFrames, maxFrame);

		glPushMatrix();
		glLoadIdentity();

		setOrthographicProjection();

		renderBitmapString(30,50, font,sFrames);

		resetPerspectiveProjection();

		glPopMatrix();

	};

	static void drawWriteMessage( int currentFile, int maxFiles )
	{
		char buffer[50];

		sprintf(buffer, "Writing File %d / %d.", currentFile, maxFiles );

		glPushMatrix();
		glLoadIdentity();

		setOrthographicProjection();

		renderBitmapString(30,50, font,buffer);

		resetPerspectiveProjection();

		glPopMatrix();

	};

	static void renderBitmapString(float x, float y, void *font,char *string)
	{
		glPushAttrib( GL_LIGHTING_BIT | GL_ENABLE_BIT );

		glColor3f(0.7f,1.0f,0.7f);
		glDisable(GL_LIGHTING);

		char *c;
		// set position to start drawing fonts
		glRasterPos2f(x, y);
		// loop all the characters in the string
		for (c=string; *c != '\0'; c++)
		{
			glutBitmapCharacter(font, *c);
		}

		glPopAttrib();
	}

	static void setOrthographicProjection() {

		// switch to projection mode
		glMatrixMode(GL_PROJECTION);
		// save previous matrix which contains the
		//settings for the perspective projection
		glPushMatrix();
		// reset matrix
		glLoadIdentity();
		// set a 2D orthographic projection
		gluOrtho2D(0, window_width, 0, window_height);
		// invert the y axis, down is positive
		//glScalef(1, -1, 1);
		// mover the origin from the bottom left corner
		// to the upper left corner
		//glTranslatef(0, -window_height, 0);
		glMatrixMode(GL_MODELVIEW);
	}
	static void resetPerspectiveProjection() {
		// set the current matrix to GL_PROJECTION
		glMatrixMode(GL_PROJECTION);
		// restore previous settings
		glPopMatrix();
		// get back to GL_MODELVIEW matrix
		glMatrixMode(GL_MODELVIEW);
	}
	
	static void printMedianFPS()
	{
		vector<double>::iterator midpoint;

		if(!fps_recorder.empty())
		{
			midpoint = fps_recorder.begin() + (fps_recorder.end() - fps_recorder.begin())/2;
			nth_element(fps_recorder.begin(), midpoint, fps_recorder.end());
			printf(" fps:%.2f", *midpoint);
		}

		if(!iter_recorder.empty())
		{
			midpoint = iter_recorder.begin() + (iter_recorder.end() - iter_recorder.begin())/2;
			nth_element(iter_recorder.begin(), midpoint, iter_recorder.end());
			printf(" itps:%.2f", *midpoint);
		}
	}

	static void resetRecorder()
	{
		fps_recorder.clear();
		iter_recorder.clear();
	}

} // end namespace

#endif
