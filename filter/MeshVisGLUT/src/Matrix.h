//No Matrix.cpp necessary
#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
public:
	union
	{
		struct
		{
			float m11, m21, m31, m41,	//OpenGL Matrix Convention
			      m12, m22, m32, m42,
				  m13, m23, m33, m43,
				  m14, m24, m34, m44;
		};

			float m[4][4];	//OpenGL Matrix Convention
			float n[16];	//OpenGL Matrix Convention
						
	};
	Matrix(){};
	~Matrix(){};
	Matrix(const Matrix& m) : m11(m.m11),m12(m.m12),m13(m.m13),m14(m.m14),
							  m21(m.m21),m22(m.m22),m23(m.m23),m24(m.m24),
							  m31(m.m31),m32(m.m32),m33(m.m33),m34(m.m34),
							  m41(m.m41),m42(m.m42),m43(m.m43),m44(m.m44){};
	Matrix(float a11, float a12, float a13, float a14,
		   float a21, float a22, float a23, float a24,
		   float a31, float a32, float a33, float a34,
		   float a41, float a42, float a43, float a44)
		   : m11(a11), m12(a12), m13(a13), m14(a14),
			 m21(a21), m22(a22), m23(a23), m24(a24),
			 m31(a31), m32(a32), m33(a33), m34(a34),
			 m41(a41), m42(a42), m43(a43), m44(a44){};

	Matrix& operator =(const Matrix& m)//normal C Matrix Convention
	{
		m11 = m.m11; m12=m.m12; m13=m.m13; m14 = m.m14;
		m21 = m.m21; m22=m.m22; m23=m.m23; m24 = m.m24;
		m31 = m.m31; m32=m.m32; m33=m.m33; m34 = m.m34;
		m41 = m.m41; m42=m.m42; m43=m.m43; m44 = m.m44;

		return(*this);
	};

	float& operator()(int row, int column)//normal C Matrix Convention
	{
		return m[column - 1][row - 1];
	};
	
	float operator()(int row, int column) const
	{
		return m[column - 1][row - 1];
	};

};

	inline Matrix operator * (const Matrix& a, const Matrix& b)
	{
		return Matrix(b.m11*a.m11 + b.m21*a.m12 + b.m31*a.m13 + b.m41*a.m14,
				      b.m12*a.m11 + b.m22*a.m12 + b.m32*a.m13 + b.m42*a.m14,
					  b.m13*a.m11 + b.m23*a.m12 + b.m33*a.m13 + b.m43*a.m14,
					  b.m14*a.m11 + b.m24*a.m12 + b.m34*a.m13 + b.m44*a.m14,

					  b.m11*a.m21 + b.m21*a.m22 + b.m31*a.m23 + b.m41*a.m24,
					  b.m12*a.m21 + b.m22*a.m22 + b.m32*a.m23 + b.m42*a.m24,
					  b.m13*a.m21 + b.m23*a.m22 + b.m33*a.m23 + b.m43*a.m24,
					  b.m14*a.m21 + b.m24*a.m22 + b.m34*a.m23 + b.m44*a.m24,

					  b.m11*a.m31 + b.m21*a.m32 + b.m31*a.m33 + b.m41*a.m34,
					  b.m12*a.m31 + b.m22*a.m32 + b.m32*a.m33 + b.m42*a.m34,
					  b.m13*a.m31 + b.m23*a.m32 + b.m33*a.m33 + b.m43*a.m34,
					  b.m14*a.m31 + b.m24*a.m32 + b.m34*a.m33 + b.m44*a.m34,
			
					  b.m11*a.m41 + b.m21*a.m42 + b.m31*a.m43 + b.m41*a.m44,
					  b.m12*a.m41 + b.m22*a.m42 + b.m32*a.m43 + b.m42*a.m44,
					  b.m13*a.m41 + b.m23*a.m42 + b.m33*a.m43 + b.m43*a.m44,
					  b.m14*a.m41 + b.m24*a.m42 + b.m34*a.m43 + b.m44*a.m44);
	};



	inline Matrix operator * (const Matrix& m, const float f)
	{
		return Matrix(m.m11*f, m.m12*f, m.m13*f, m.m14*f,
					  m.m21*f, m.m22*f, m.m23*f, m.m24*f,
					  m.m31*f, m.m32*f, m.m33*f, m.m34*f,
					  m.m41*f, m.m42*f, m.m43*f, m.m44*f);
	};
	inline Matrix operator * (const float f, const Matrix& m)
	{
		return Matrix(m.m11*f, m.m12*f, m.m13*f, m.m14*f,
					  m.m21*f, m.m22*f, m.m23*f, m.m24*f,
					  m.m31*f, m.m32*f, m.m33*f, m.m34*f,
					  m.m41*f, m.m42*f, m.m43*f, m.m44*f);
	};

	inline Vector operator * (const Matrix& m, const Vector& v)
	{	
	
		Vector Result = Vector(m.m11*v.x + m.m12 * v.y + m.m13 * v.z + m.m14,
							   m.m21*v.x + m.m22 * v.y + m.m23 * v.z + m.m24,
							   m.m31*v.x + m.m32 * v.y + m.m33 * v.z + m.m34);

		float w = m.m41 * v.x + m.m42*v.y + m.m43*v.z + m.m44;
		if(w!=1.0f) Result /= w;
		return Result;
	};
	inline Vector operator * (const Vector& v, const Matrix& m)
	{	//do not use this method.. only if MatrixTranslate, MatrixDet and MatrixInvert is changed
	
		Vector Result = Vector(m.m11*v.x + m.m21 * v.y + m.m31 * v.z + m.m41,
							   m.m12*v.x + m.m22 * v.y + m.m32 * v.z + m.m42,
							   m.m13*v.x + m.m23 * v.y + m.m33 * v.z + m.m43);

		float w = m.m14 * v.x + m.m24*v.y + m.m34*v.z + m.m44;
		if(w!=1.0f) Result /= w;
		return Result;
	};
	inline Matrix operator / (const Matrix& a, const float f)
	{
		return Matrix(a.m11/f , a.m12/f, a.m13/f, a.m14/f,
			          a.m21/f , a.m22/f, a.m23/f, a.m24/f,
					  a.m31/f , a.m32/f, a.m33/f, a.m34/f,
					  a.m41/f , a.m42/f, a.m43/f, a.m44/f);
	};
	inline Matrix operator + (const Matrix& a, const Matrix& b)
	{
		return Matrix(a.m11+b.m11, a.m12+b.m12, a.m13+b.m13, a.m14+b.m14,
					  a.m21+b.m21, a.m22+b.m22, a.m23+b.m23, a.m24+b.m24,
					  a.m31+b.m31, a.m32+b.m32, a.m33+b.m33, a.m34+b.m34,
					  a.m41+b.m41, a.m42+b.m42, a.m43+b.m43, a.m44+b.m44);
	};
	inline Matrix operator - (const Matrix& a, const Matrix& b)
	{
		return Matrix(a.m11-b.m11, a.m12-b.m12, a.m13-b.m13, a.m14-b.m14,
					  a.m21-b.m21, a.m22-b.m22, a.m23-b.m23, a.m24-b.m24,
					  a.m31-b.m31, a.m32-b.m32, a.m33-b.m33, a.m34-b.m34,
					  a.m41-b.m41, a.m42-b.m42, a.m43-b.m43, a.m44-b.m44);
	};

	inline Matrix MatrixIdentity()
	{
		return Matrix(1.0f,0.0f,0.0f,0.0f,
					  0.0f,1.0f,0.0f,0.0f,
					  0.0f,0.0f,1.0f,0.0f,
					  0.0f,0.0f,0.0f,1.0f);
	};

	inline Matrix MatrixTranslation(const Vector& v)
	{
		return Matrix(1.0f, 0.0f, 0.0f,  v.x,
					  0.0f, 1.0f, 0.0f,  v.y,
					  0.0f, 0.0f, 1.0f,  v.z,
					  0.0f, 0.0f, 0.0f, 1.0f);

	};



	inline Matrix MatrixRotationX(const float f)
	{
		Matrix Result;
		Result.m11 = 1.0f;
		Result.m12 = 0.0f;
		Result.m13 = 0.0f;
		Result.m14 = 0.0f;

		Result.m21 = 0.0f;
		Result.m22 = cosf(f);
		Result.m23 = -sinf(f);
		Result.m24 = 0.0f;

		Result.m31 = 0.0f;
		Result.m32 = -Result.m23;
		Result.m33 =  Result.m22;
		Result.m34 = 0.0f;

		Result.m41 = 0.0f;
		Result.m42 = 0.0f;
		Result.m43 = 0.0f;
		Result.m44 = 1.0f;

		return Result;
	};

	inline Matrix MatrixRotationY(const float f)
	{
	Matrix Result;
	Result.m11 =  cosf(f);
	Result.m12 = 0.0f;
	Result.m13 =  sinf(f);
	Result.m14 = 0.0f;

	Result.m21 = 0.0f;
	Result.m22 = 1.0f;
	Result.m23 = 0.0f;
	Result.m24 = 0.0f;

	Result.m31 = -Result.m13;
	Result.m32 = 0.0f;
	Result.m33 =  Result.m11;
	Result.m34 = 0.0f;

	Result.m41 = 0.0f;
	Result.m42 = 0.0f;
	Result.m43 = 0.0f;
	Result.m44 = 1.0f;

	return Result;
};

	inline Matrix MatrixRotationZ(const float f)
	{
		Matrix Result;
		Result.m11 =  cosf(f);
		Result.m12 = -sinf(f);
		Result.m13 = 0.0f;
		Result.m14 = 0.0f;

		Result.m21 = -Result.m12;
		Result.m22 =  Result.m11;
		Result.m23 = 0.0f;
		Result.m24 = 0.0f;

		Result.m31 = 0.0f;
		Result.m32 = 0.0f;
		Result.m33 = 1.0f;
		Result.m34 = 0.0f;

		Result.m41 = 0.0f;
		Result.m42 = 0.0f;
		Result.m43 = 0.0f;
		Result.m44 = 1.0f;

		return Result;
	};



	inline Matrix MatrixRotation(const float x, const float y, const float z)
	{
		return MatrixRotationZ(z)*MatrixRotationY(y)*MatrixRotationX(x);
	};

	inline Matrix MatrixRotation(const Vector& v)
	{
		return MatrixRotationZ(v.z)*MatrixRotationY(v.y)*MatrixRotationX(v.x);
	};
	inline Matrix MatrixScaling(const float x, const float y, const float z)
	{
		return Matrix(   x, 0.0f, 0.0f, 0.0f,
					  0.0f,    y, 0.0f, 0.0f,
					  0.0f, 0.0f,    z, 0.0f,
					  0.0f, 0.0f, 0.0f, 1.0f);
	};
	inline Matrix MatrixScaling(const Vector& v)
	{
		return Matrix( v.x, 0.0f, 0.0f, 0.0f,
					  0.0f,  v.y, 0.0f, 0.0f,
					  0.0f, 0.0f,  v.z, 0.0f,
					  0.0f, 0.0f, 0.0f, 1.0f);
	};

	inline float  MatrixDet(const Matrix& m)
	{
		return m.m11 * ( m.m22 * m.m33 - m.m23 * m.m32) -
			   m.m12 * ( m.m21 * m.m33 - m.m23 * m.m31) +
			   m.m13 * ( m.m21 * m.m32 - m.m22 * m.m31) ;
	};


	inline Matrix MatrixInvert(const Matrix& m)
	{
		float InvDet = MatrixDet(m);
		if(InvDet == 0.0f) return MatrixIdentity();
		InvDet = 1.0f / InvDet;

		Matrix Result;
		Result.m11 =  InvDet * ( m.m22 * m.m33 - m.m23 * m.m32);
		Result.m12 = -InvDet * ( m.m12 * m.m33 - m.m13 * m.m32);
		Result.m13 =  InvDet * ( m.m12 * m.m23 - m.m13 * m.m22);
		Result.m41 =  0.0f;

		Result.m21 = -InvDet * ( m.m21 * m.m33 - m.m23 * m.m31);
		Result.m22 =  InvDet * ( m.m11 * m.m33 - m.m13 * m.m31);
		Result.m23 = -InvDet * ( m.m11 * m.m23 - m.m13 * m.m21);
		Result.m42 =  0.0f;
	
		Result.m31 =  InvDet * ( m.m21 * m.m32 - m.m22 * m.m31);
		Result.m32 = -InvDet * ( m.m11 * m.m32 - m.m12 * m.m31);
		Result.m33 =  InvDet * ( m.m11 * m.m22 - m.m12 * m.m21);
		Result.m43 =  0.0f;

		Result.m14 = -(m.m14 * Result.m11 + m.m24 * Result.m12 + m.m34 * Result.m13);
		Result.m24 = -(m.m14 * Result.m21 + m.m24 * Result.m22 + m.m34 * Result.m23);
		Result.m34 = -(m.m14 * Result.m31 + m.m24 * Result.m32 + m.m34 * Result.m33);
		Result.m44 = 1.0f;

		return Result;
	};
	inline Matrix MatrixTranspose(const Matrix& m)
	{
		return Matrix(m.m11, m.m21, m.m31, m.m41,
				    m.m12, m.m22, m.m32, m.m42,
					  m.m13, m.m23, m.m33, m.m43,
					  m.m14, m.m24, m.m34, m.m44);
	};
	inline Matrix MatrixCamera(const Vector& c, const Vector& d, const Vector& up)
	{	
		Vector r = VectorCross(d,up);
		Vector u = VectorCross(r,d);

		return Matrix( r.x,  r.y,  r.z, -1*VectorDot(r,c),
			           u.x,  u.y,  u.z, -1*VectorDot(u,c),
					  -d.x, -d.y, -d.z,  VectorDot(d,c),
					  0.0f, 0.0f, 0.0f,            1.0f);
	};
	inline Matrix MatrixFrustum(const float fovx, const float fovy, const float n, const float f)
	{
		float t =  n*tanf(fovy/2.0f);
		float b = -t;
		float r =  n*tanf(fovx/2.0f);
		float l = -r;
		return Matrix((2*n)/(r-l) ,          0 ,   (r+l)/(r-l)  ,               0 ,
			                    0 , (2*n)/(t-b),   (t+b)/(t-b)  ,               0 ,
						        0 ,          0 , -((f+n)/(f-n)) ,-((2*f*n)/(f-n)) ,
								0 ,          0 ,             -1 ,               0 );
	};

	inline Matrix MatrixInvertFrustum(const Matrix& m)
	{
		return Matrix(1/m.m11,		  0,		0,		m.m13/m.m11,
			                0,  1/m.m22,        0,      m.m23/m.m22,
							0,        0,        0,               -1,
							0,        0,  1/m.m34,      m.m33/m.m34);

	};
#endif