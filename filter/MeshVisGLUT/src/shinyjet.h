// ShinyJet.cpp
// OpenGL SuperBible
// Demonstrates OpenGL Lighting
// Program by Richard S. Wright Jr.

#ifndef SHINY_JET_H
#define SHINY_JET_H

// includes, GL
#include <GL/glew.h>
#include <GL/glut.h>

#include <vista_frame/common/glh_linear.h>

using namespace glh;

// Called to draw scene
void RenderJet(void)
	{

	glh::vec3 vNormal;	// Storeage for calculated surface normal

	// Clear the window with current clearing color
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Save the matrix state and do the rotations
	glPushMatrix();

		glScalef(.01,.01,.01);

	// Nose Cone - Points straight down
    // Set material color
	glColor3ub(128, 128, 128);
	glBegin(GL_TRIANGLES);
        glNormal3f(0.0f, -1.0f, 0.0f);
		glNormal3f(0.0f, -1.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 60.0f);
		glVertex3f(-15.0f, 0.0f, 30.0f);
		glVertex3f(15.0f,0.0f,30.0f);
                
	
        // Verticies for this panel
        {
		vec3 vPoints0 = vec3( 15.0f, 0.0f,  30.0f);
		vec3 vPoints1 = vec3( 0.0f,  15.0f, 30.0f);
		vec3 vPoints2 = vec3( 0.0f,  0.0f,  60.0f);
			                                                         

        // Calculate the normal for the plane
        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
		glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }	


        {
		vec3 vPoints0 = vec3( 0.0f, 0.0f,  60.0f);
		vec3 vPoints1 = vec3( 0.0f,  15.0f, 30.0f);
		vec3 vPoints2 = vec3( -15.0f,  0.0f,  30.0f);

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
		glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
		}


        // Body of the Plane ////////////////////////
        {
        vec3 vPoints0 = vec3 ( -15.0f, 0.0f, 30.0f );
		vec3 vPoints1 = vec3 ( 0.0f, 15.0f, 30.0f );
		vec3 vPoints2 = vec3 ( 0.0f, 0.0f, -56.0f );

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
		glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                	
        {
        vec3 vPoints0 = vec3( 0.0f, 0.0f, -56.0f );
        vec3 vPoints1 = vec3( 0.0f, 15.0f, 30.0f );
        vec3 vPoints2 = vec3( 15.0f, 0.0f, 30.0f );
	
        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                		
    
		glNormal3f(0.0f, -1.0f, 0.0f);
		glVertex3f(15.0f,0.0f,30.0f);
		glVertex3f(-15.0f, 0.0f, 30.0f);
		glVertex3f(0.0f, 0.0f, -56.0f);
    
        ///////////////////////////////////////////////
        // Left wing
        // Large triangle for bottom of wing
        {
       vec3 vPoints0 = vec3( 0.0f,2.0f,27.0f );
		vec3 vPoints1 = vec3( -60.0f, 2.0f, -8.0f );
		vec3 vPoints2 = vec3( 60.0f, 2.0f, -8.0f );

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                
                        	
        {
        vec3 vPoints0 = vec3( 60.0f, 2.0f, -8.0f );
		vec3 vPoints1 = vec3( 0.0f, 7.0f, -8.0f );
		vec3 vPoints2 = vec3( 0.0f,2.0f,27.0f );
                
        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                
        {
        vec3 vPoints0 = vec3( 60.0f, 2.0f, -8.0f );
		vec3 vPoints1 = vec3( -60.0f, 2.0f, -8.0f );
		vec3 vPoints2 = vec3( 0.0f,7.0f,-8.0f );

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                
        {
        vec3 vPoints0 = vec3( 0.0f,2.0f,27.0f);
        vec3 vPoints1 = vec3(0.0f, 7.0f, -8.0f);
        vec3 vPoints2 = vec3(-60.0f, 2.0f, -8.0f);

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                
                        
        // Tail section///////////////////////////////
        // Bottom of back fin
		glNormal3f(0.0f, -1.0f, 0.0f);
		glVertex3f(-30.0f, -0.50f, -57.0f);
		glVertex3f(30.0f, -0.50f, -57.0f);
		glVertex3f(0.0f,-0.50f,-40.0f);

        {
        vec3 vPoints0 = vec3( 0.0f,-0.5f,-40.0f );
        vec3 vPoints1 = vec3( 30.0f, -0.5f, -57.0f);
        vec3 vPoints2 = vec3( 0.0f, 4.0f, -57.0f );

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                
                        
        {
       vec3 vPoints0 = vec3( 0.0f, 4.0f, -57.0f );
       vec3 vPoints1 = vec3(  -30.0f, -0.5f, -57.0f );
       vec3 vPoints2 = vec3( 0.0f,-0.5f,-40.0f );

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }

        {
         vec3 vPoints0 = vec3( 30.0f,-0.5f,-57.0f );
		 vec3 vPoints1 = vec3( -30.0f, -0.5f, -57.0f );
		 vec3 vPoints2 = vec3( 0.0f, 4.0f, -57.0f );

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                
        {
        vec3 vPoints0 = vec3( 30.0f,0.5f,-40.0f );
		vec3 vPoints1 = vec3( 33.0f, 0.5f, -57.0f );
		vec3 vPoints2 = vec3( 30.0f, 25.0f, -65.0f );

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                
                        
        {
        vec3 vPoints0 = vec3( 0.0f, 25.0f, -65.0f );
		vec3 vPoints1 = vec3( -3.0f, 0.5f, -57.0f);
		vec3 vPoints2 = vec3( 0.0f,0.5f,-40.0f );

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                
        {
        vec3 vPoints0 = vec3( 3.0f,0.5f,-57.0f );
		vec3 vPoints1 = vec3( -3.0f, 0.5f, -57.0f );
		vec3 vPoints2 = vec3( 0.0f, 25.0f, -65.0f );

        vNormal = vPoints0.cross( vPoints1 ).cross( vPoints2);
        glNormal3fv(vNormal);
		glVertex3fv(vPoints0);
		glVertex3fv(vPoints1);
		glVertex3fv(vPoints2);
        }
                
                
        glEnd();
                
    	// Restore the matrix state
	glPopMatrix();

}

// This function does any needed initialization on the rendering
// context. 
void SetupRC()
    {
    // Light values and coordinates
    GLfloat  ambientLight[] = { 0.3f, 0.3f, 0.3f, 1.0f };
    GLfloat  diffuseLight[] = { 0.7f, 0.7f, 0.7f, 1.0f };
    GLfloat  specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat  specref[] = { 1.0f, 1.0f, 1.0f, 1.0f };

    glEnable(GL_DEPTH_TEST);	// Hidden surface removal
    glFrontFace(GL_CCW);		// Counter clock-wise polygons face out
    glEnable(GL_CULL_FACE);		// Do not calculate inside of jet

    // Enable lighting
    glEnable(GL_LIGHTING);

    // Setup and enable light 0
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientLight);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuseLight);
    glLightfv(GL_LIGHT0,GL_SPECULAR, specular);
    glEnable(GL_LIGHT0);

    // Enable color tracking
    glEnable(GL_COLOR_MATERIAL);
	
    // Set Material properties to follow glColor values
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    // All materials hereafter have full specular reflectivity
    // with a high shine
    glMaterialfv(GL_FRONT, GL_SPECULAR, specref);
    glMateriali(GL_FRONT, GL_SHININESS, 128);
    
    // Light blue background
   
    glEnable(GL_NORMALIZE);

    }

#endif