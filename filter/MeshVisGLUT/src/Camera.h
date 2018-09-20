#ifndef CAMERA_H
#define CAMERA_H

#include <GL/freeglut.h>
#include "Vector.h"
#include "Matrix.h"

class Camera
{
public:
	float angle_y,angle_z;
	int ar,br;
	Vector Position;
	Vector LookAt;
	Vector Up;
	
	Camera()
	{
		Position= Vector(0.0f,   0.0f,   0.0f);
		LookAt  = Vector(0.0f,   0.0f,  -1.0f);
		Up      = Vector(0.0f,   1.0f,   0.0f);
	}

	Camera(const Vector& p, const Vector& l, const Vector& u)
	{
		Position= p;
		LookAt  = l;
		Up      = u;
	}

	void CameraMove(float speed)
	{
		Vector vVector = LookAt-Position;	// Get the view vector
	
		// forward positive camera speed and backward negative camera speed.
		Position.x  = Position.x  + vVector.x * speed;
		Position.z  = Position.z  + vVector.z * speed;

		LookAt.x = LookAt.x + vVector.x * speed;
		LookAt.z = LookAt.z + vVector.z * speed;

		Position.y = Position.y + vVector.y * speed;
		LookAt.y   = LookAt.y + vVector.y * speed;
	}
	void CameraMoveUp(float speed)
	{
		Position.y = Position.y + speed;
		  LookAt.y =   LookAt.y + speed;
	}
	void CameraStrafe(float speed)
	{
		Vector vVector = VectorNormalize(LookAt - Position);	// Get the view vector
		Vector vOrthoVector;              // Orthogonal vector for the view vector

		vOrthoVector.x = -vVector.z;
		vOrthoVector.z =  vVector.x;

		// left negative -cameraspeed and right positive +cameraspeed.
		Position.x  = Position.x  - vOrthoVector.x * speed;
		Position.z  = Position.z  - vOrthoVector.z * speed;
		LookAt.x = LookAt.x - vOrthoVector.x * speed;
		LookAt.z = LookAt.z - vOrthoVector.z * speed;
	}
	
	
	void CameraRotateViewY(float speed)
	{
		Vector vVector = VectorNormalize(LookAt - Position);	// Get the view vector

		LookAt.z = (float)(Position.z + sin(speed)*vVector.x + cos(speed)*vVector.z);
		LookAt.x = (float)(Position.x + cos(speed)*vVector.x - sin(speed)*vVector.z);
		Up = VectorNormalize(MatrixRotationY(speed)*Up);
		


		//upvector?

	}

	void CameraRotateViewX(float speed)
	{
		
		Vector vVector = VectorNormalize(LookAt - Position);
		LookAt.z = (float)(Position.z + sin(speed)*vVector.y + cos(speed)*vVector.z);
		LookAt.y = (float)(Position.y + cos(speed)*vVector.y - sin(speed)*vVector.z);
		Up = VectorNormalize(MatrixRotationX(speed)*Up);
	}

	void CameraMove(int dx, int dy)
	{
		if( (dx == 0) && (dy == 0) ) return;

		//float angle_y  = 0.0f;				
		//float angle_z  = 0.0f;		

		Vector Temp = LookAt;
// 		CameraRotateViewX((float) (dx / 1000));		
		angle_y = (float)( (dx) ) / 100;	//500	
		angle_z = (float)( (-dy) ) / 100;	//500

		LookAt.y += angle_z * 0.25f;
		CameraRotateViewY(angle_y); // Rotate
		LookAt = Position+(VectorNormalize(LookAt-Position));
	}
	
	void CameraMouseMove(int wndWidth, int wndHeight, int x, int y)
	{
		
		Vector Temp = LookAt;
		int mid_x = wndWidth/2;	
		int mid_y = wndHeight/2;	
		
		float angle_y  = 0.0f;				
		float angle_z  = 0.0f;		

		if( (x == mid_x) && (y == mid_y) ) return;
						
		angle_y = (float)( (x-mid_x) ) / 1000;	//500	
		angle_z = (float)( (mid_y-y) ) / 1000;	//500

		// The higher the value is the faster the camera looks around.
		LookAt.y += angle_z * 0.25f;
		CameraRotateViewY(angle_y); // Rotate
		//CameraRotateViewX(angle_z);
		
		LookAt = Position+(VectorNormalize(LookAt-Position));
		
		
		// limit the rotation around the x-axis
		//if((LookAt.y - Position.y) > 8)   LookAt.y = Position.y - 8;
		//if((LookAt.y - Position.y) <-8)   LookAt.y = Position.y + 8;
	
		
		
		
		CameraRotateViewY(angle_y); // Rotate
		//CameraRotateViewX(angle_z);
		
		LookAt = Position+(VectorNormalize(LookAt-Position));

		//glutWarpPointer( mid_x, mid_y);
              
	}

};


#endif
