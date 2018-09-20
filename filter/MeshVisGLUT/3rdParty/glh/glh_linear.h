/*
    glh - is a platform-indepenedent C++ OpenGL helper library 


    Copyright (c) 2000 Cass Everitt
	Copyright (c) 2000 NVIDIA Corporation
    All rights reserved.

    Redistribution and use in source and binary forms, with or
	without modification, are permitted provided that the following
	conditions are met:

     * Redistributions of source code must retain the above
	   copyright notice, this list of conditions and the following
	   disclaimer.

     * Redistributions in binary form must reproduce the above
	   copyright notice, this list of conditions and the following
	   disclaimer in the documentation and/or other materials
	   provided with the distribution.

     * The names of contributors to this software may not be used
	   to endorse or promote products derived from this software
	   without specific prior written permission. 

       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
	   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
	   REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
	   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
	   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
	   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
	   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
	   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
	   POSSIBILITY OF SUCH DAMAGE. 


    Cass Everitt - cass@r3.nu
*/

/*
glh_linear.h
*/

// Author:  Cass W. Everitt

#ifndef GLH_LINEAR_H
#define GLH_LINEAR_H

#include <memory.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <ostream>

// only supports float for now...
#define GLH_REAL_IS_FLOAT

#ifdef GLH_REAL_IS_FLOAT
# define GLH_REAL float
//# define GLH_REAL_NAMESPACE ns_float
#endif

#ifdef _WIN32
# define TEMPLATE_FUNCTION
#else
# define TEMPLATE_FUNCTION <>
#endif

#define     GLH_QUATERNION_NORMALIZATION_THRESHOLD  64

#define     GLH_RAD_TO_DEG      GLH_REAL(57.2957795130823208767981548141052)
#define     GLH_DEG_TO_RAD      GLH_REAL(0.0174532925199432957692369076848861)
#define     GLH_ZERO            GLH_REAL(0.0)
#define     GLH_ONE             GLH_REAL(1.0)
#define     GLH_TWO             GLH_REAL(2.0)
#define     GLH_EPSILON         GLH_REAL(10e-6)
#define     GLH_PI              GLH_REAL(3.1415926535897932384626433832795)    

#define     equivalent(a,b)     (((a < b + GLH_EPSILON) && (a > b - GLH_EPSILON)) ? true : false)
/*
namespace glh{
	template <int N, class T>	class vec;}

template<int N, class T>
inline std::ostream& 
operator<<(std::ostream& os, const glh::vec<N, T>& v) {
	for(int i=0; i<N-1; i++) os << v[i] << " ";
	os << v[N-1];
	return os;
}
*/

/// glh - is a platform-indepenedent C++ OpenGL helper library.
/// The following classes are provided:
///
/// class vec<N,T> vector class with type T of dimension N
///
/// class line
///
///	class plane
///
/// class matrix4
///
///	class quaternion
///   
///
///    Copyright (c) 2000 Cass Everitt
///	Copyright (c) 2000 NVIDIA Corporation
///    All rights reserved.
///
///    Redistribution and use in source and binary forms, with or
///	without modification, are permitted provided that the following
///	conditions are met:
///
///     * Redistributions of source code must retain the above
///	   copyright notice, this list of conditions and the following
///	   disclaimer.
///
///     * Redistributions in binary form must reproduce the above
///	   copyright notice, this list of conditions and the following
///	   disclaimer in the documentation and/or other materials
///	   provided with the distribution.
///
///     * The names of contributors to this software may not be used
///	   to endorse or promote products derived from this software
///	   without specific prior written permission. 
///
///       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
///	   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
///	   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
///	   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
///	   REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
///	   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
///	   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
///	   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
///	   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
///	   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
///	   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
///	   POSSIBILITY OF SUCH DAMAGE. 
///
///
///    Cass Everitt - cass@r3.nu
namespace glh
{

	inline GLH_REAL to_degrees(GLH_REAL radians) { return radians*GLH_RAD_TO_DEG; }
	inline GLH_REAL to_radians(GLH_REAL degrees) { return degrees*GLH_DEG_TO_RAD; }

  /// Vector class	
	template <int N, class T>	
	class vec
	{				
    public:
    /// Returns vector dimensions.
		int size() const { return N; }
		
		vec(const T & t = T()) 
		{ for(int i = 0; i < N; i++) v[i] = t; }
		vec(const T * tp)
		{ for(int i = 0; i < N; i++) v[i] = tp[i]; }
		
		vec(const T t1,const T t2) {v[0] = t1;v[1] = t2;}
		vec(const T t1,const T t2,const T t3) {v[0] = t1;v[1] = t2;v[2] = t3;}
		vec(const T t1,const T t2,const T t3,const T t4) {v[0] = t1;v[1] = t2;v[2] = t3;v[3] = t4;}
		
		

    /// Returns pointer to vector components.
		const T * get_value() const
		{ return v; }
		
	
    /// Dot product. Vec lhs.dot(rhs) returns result = lhs * rhs.
    /// \param rhs right-hand-side of dot product.
    /// \return results of the dot product.
		T dot( const vec<N,T> & rhs ) const
		{ 
			T r = 0;
			for(int i = 0; i < N; i++) r += v[i]*rhs[i];
			return r;
		}
	
    /// Returns length of vector.
		T length() const
		{
			T r = 0;
			for(int i = 0; i < N; i++) r += v[i]*v[i]; 
			return T(sqrt(r));
		}	
	
    /// Returns the squared length of the vector.
		T square_norm() const
		{
			T r = 0;
			for(int i = 0; i < N; i++) r += v[i]*v[i]; 
			return r;
		}	
	
    /// Negates all components of the vector, vector v becomes -v.
		void  negate()
		{ for(int i = 0; i < N; i++) v[i] = -v[i]; }
		
	  /// Normalizes the vector, hence the length of the result vector is equal to 1.
    /// \return Returns the length of the original (not normalized) vector.
		T normalize() 
		{ 
			T sum(0);
			for(int i = 0; i < N; i++) 
                sum += v[i]*v[i];
			sum = T(sqrt(sum));
            if (sum > GLH_EPSILON)
			    for(int i = 0; i < N; i++) 
                    v[i] /= sum;
			return sum;
		}
		
		/// Sets the vector component.
    /// \param rhs vector components.
    /// \return Returns a reference to itself.
		vec<N,T> & set_value( const T * rhs )
		{ for(int i = 0; i < N; i++) v[i] = rhs[i]; return *this; }
		

		/// cast to array of scalars of type T
		inline operator T*() { return v; }

		/// cast to const T array
		inline operator const T*() const { return v; }

    /// Component wise access (read/write), 
    /// For vector v the x-component v.x is accessed by v[0].
		T & operator [] ( int i )
		{ return v[i]; }
    
    /// Component wise access (read-only). 
    /// For vector v the x-component v.x is accessed by v[0].
		const T & operator [] ( int i ) const
		{ return v[i]; }

		vec<N,T> & operator *= ( T d )
		{ for(int i = 0; i < N; i++) v[i] *= d; return *this;}
		
		vec<N,T> & operator *= ( const vec<N,T> & u )
		{ for(int i = 0; i < N; i++) v[i] *= u[i]; return *this;}
		
		vec<N,T> & operator /= ( T d )
		{ if(d == 0) return *this; for(int i = 0; i < N; i++) v[i] /= d; return *this;}
		
		vec<N,T> & operator += ( const vec<N,T> & u )
		{ for(int i = 0; i < N; i++) v[i] += u[i]; return *this;}
		
		vec<N,T> & operator -= ( const vec<N,T> & u )
		{ for(int i = 0; i < N; i++) v[i] -= u[i]; return *this;}
		
		/// Unary minus operator.
		vec<N,T> operator - () const
		{ vec<N,T> rv = v; rv.negate(); return rv; }
	
    /// Plus operator.
		vec<N,T> operator + ( const vec<N,T> &v) const
		{ vec<N,T> rt(*this); return rt += v; }
	
    /// Minus operator.
		vec<N,T> operator - ( const vec<N,T> &v) const
		{ vec<N,T> rt(*this); return rt -= v; }
	
    /// Multiply with scalar operator.
		vec<N,T> operator * ( T d) const
		{ vec<N,T> rt(*this); return rt *= d; }
		
		//friend bool operator == TEMPLATE_FUNCTION ( const vec<N,T> &v1, const vec<N,T> &v2 );
		//friend bool operator != TEMPLATE_FUNCTION ( const vec<N,T> &v1, const vec<N,T> &v2 );				
		
	protected:
		T v[N];
	};
	
	
	
	// vector friend operators
	
	template <int N, class T> inline
		vec<N,T> operator * ( const vec<N,T> & b, T d )
	{
		vec<N,T> rt(b);
		return rt *= d;
	}

	template <int N, class T> inline
		vec<N,T> operator * ( T d, const vec<N,T> & b )
	{ return b*d; }
	
	template <int N, class T> inline
		vec<N,T> operator * ( const vec<N,T> & b, const vec<N,T> & d )
	{
		vec<N,T> rt(b);
		return rt *= d;
	}

	template <int N, class T> inline
		vec<N,T> operator / ( const vec<N,T> & b, T d )
	{ vec<N,T> rt(b); return rt /= d; }
	
	template <int N, class T> inline
		vec<N,T> operator + ( const vec<N,T> & v1, const vec<N,T> & v2 )
	{ vec<N,T> rt(v1); return rt += v2; }
	
	template <int N, class T> inline
		vec<N,T> operator - ( const vec<N,T> & v1, const vec<N,T> & v2 )
	{ vec<N,T> rt(v1); return rt -= v2; }
	
  /// Componentwise equal operator.	
	template <int N, class T> inline
		bool operator == ( const vec<N,T> & v1, const vec<N,T> & v2 )
	{
		for(int i = 0; i < N; i++)
			if(v1[i] != v2[i])
				return false;
			return true;
	}
  
  /// Componentwise not-equal operator.	
	template <int N, class T> inline
		bool operator != ( const vec<N,T> & v1, const vec<N,T> & v2 )
	{ return !(v1 == v2); }
	

	typedef vec<3,unsigned char> vec3ub;
	typedef vec<4,unsigned char> vec4ub;

	typedef vec<2,int> vec2i;
	typedef vec<3,int> vec3i;
	typedef vec<4,int> vec4i;




//	namespace GLH_REAL_NAMESPACE
//	{
	typedef GLH_REAL real;

	class line;
	class plane;
	class matrix4;
	class quaternion;
	typedef quaternion rotation; 
 
  /// 2D Vector
	class vec2 : public vec<2,real>
	{
    public:
		vec2(const real & t = real()) : vec<2,real>(t)
		{}
		vec2(const vec<2,real> & t) : vec<2,real>(t)
		{}
		vec2(const real * tp) : vec<2,real>(tp)
		{}
	
    /// Constructor with x and y components.
		vec2(real x, real y )
		{ v[0] = x; v[1] = y; }
		
		void get_value(real & x, real & y) const
		{ x = v[0]; y = v[1]; }
		
		vec2 & set_value( const real & x, const real & y)
		{ v[0] = x; v[1] = y; return *this; }
		
	};
	
  /// 3D Vector	
	class vec3 : public vec<3,real>
	{
    public:
		vec3(const real & t = real()) : vec<3,real>(t)
		{}
		vec3(const vec<3,real> & t) : vec<3,real>(t)
		{}
		vec3(const real * tp) : vec<3,real>(tp)
		{}
	
    /// Constructor with x,y,z components.
		vec3(real x, real y, real z)
		{ v[0] = x; v[1] = y; v[2] = z; }
		
		void get_value(real & x, real & y, real & z) const
		{ x = v[0]; y = v[1]; z = v[2]; }
	
    /// Cross product. Leaves this vector unchanged.
    /// \param rhs right-hand-side of cross product.
    /// \return Returns the cross product of this vector with rhs
		vec3 cross( const vec3 &rhs ) const
		{
			vec3 rt;
			rt[0] = v[1]*rhs[2]-v[2]*rhs[1];
			rt[1] = v[2]*rhs[0]-v[0]*rhs[2];
			rt[2] = v[0]*rhs[1]-v[1]*rhs[0];	
			return rt;
		}
		
		vec3 & set_value( const real & x, const real & y, const real & z)
		{ v[0] = x; v[1] = y; v[2] = z; return *this; }
		
	};

  	
    /// 4D Vector
    class vec4 : public vec<4,real>
    {
    public:
        vec4(const real & t = real()) : vec<4,real>(t)
        {}
        vec4(const vec<4,real> & t) : vec<4,real>(t)
        {}

        vec4(const vec<3,real> & t, real fourth)

        { v[0] = t[0]; v[1] = t[1]; v[2] = t[2]; v[3] = fourth; }
        vec4(const real * tp) : vec<4,real>(tp)
        {}
        /// Constructor with x,y,z,w components.
        vec4(real x, real y, real z, real w)
        { v[0] = x; v[1] = y; v[2] = z; v[3] = w; }

        void get_value(real & x, real & y, real & z, real & w) const
        { x = v[0]; y = v[1]; z = v[2]; w = v[3]; }
  
        vec4 & set_value( const real & x, const real & y, const real & z, const real & w)
        { v[0] = x; v[1] = y; v[2] = z; v[3] = w; return *this; }
    };

    /// Homogenizes 4D vector, i.e. divides x,y,z components by w
    /// \param v vector to be homogenized.
    /// \return 3D vector r=(v.x/v.w,v.y/v.w,v.z/v.w)
    inline
    vec3 homogenize(const vec4 & v)
    {
      vec3 rt;
      assert(v[3] != GLH_ZERO);
      rt[0] = v[0]/v[3];
      rt[1] = v[1]/v[3];
      rt[2] = v[2]/v[3];
      return rt;
    }
  

    /// Line in 3D, represented by two points
    class line
    {
    public:
  
        line()
        { set_value(vec3(0,0,0),vec3(0,0,1)); }

        /// Constructor with start and end point
        /// \param p0 start point of line.
        /// \param p1 end point of line
        line( const vec3 & p0, const vec3 &p1)
        { set_value(p0,p1); }

        void set_value( const vec3 &p0, const vec3 &p1)
        {
          position = p0;
          direction = p1-p0;
          direction.normalize();
        }
 
        /// Compute the closest points on two lines. Uses the algorithm of Brian Jean.
        /// \param line2 other line
        /// \param pointOnThis will contain the closest point on this line
        /// \param pointOnThat will contain the closest point on line2
        /// \return false if the lines are parallel, otherwise true
        bool get_closest_points(const line &line2, 
					          vec3 &pointOnThis,
					          vec3 &pointOnThat)
        {
  
          // quick check to see if parallel -- if so, quit.
          if(fabs(direction.dot(line2.direction)) == 1.0)
	          return 0;
          line l2 = line2;
  
          // Algorithm: Brian Jean
          // 
          register real u;
          register real v;
          vec3 Vr = direction;
          vec3 Vs = l2.direction;
          register real Vr_Dot_Vs = Vr.dot(Vs);
          register real detA = real(1.0 - (Vr_Dot_Vs * Vr_Dot_Vs));
          vec3 C = l2.position - position;
          register real C_Dot_Vr =  C.dot(Vr);
          register real C_Dot_Vs =  C.dot(Vs);
  
          u = (C_Dot_Vr - Vr_Dot_Vs * C_Dot_Vs)/detA;
          v = (C_Dot_Vr * Vr_Dot_Vs - C_Dot_Vs)/detA;
  
          pointOnThis = position;
          pointOnThis += direction * u;
          pointOnThat = l2.position;
          pointOnThat += l2.direction * v;
  
          return 1;
        }
 
        /// Computes the closest point between a line and a point.
        /// \param point other point 
        /// \return closest point on this line to point
        vec3 get_closest_point(const vec3 &point)
        {
          vec3 np = point - position;
          vec3 rp = direction*direction.dot(np)+position;
          return rp;
        }
 
        /// Returns start point of line.
        const vec3 & get_position() const {return position;}

        /// Returns direction of line.
        const vec3 & get_direction() const {return direction;}
  
    //protected:
        vec3 position;
        vec3 direction;
    };
  
  














  
  











  /// 4x4 Matrix
  class matrix4
  {
    
  public:
    /// Standard constructor creates identity matrix.    
    matrix4() { make_identity(); }
   
    /// Constructor that creates a matrix where all components are the same.
    /// Hence, to create a zero matrix: matrix4 m(0)
    /// \param r Scalar value that will be assigned to every components.
	matrix4( real r ) 
	{ set_value(r); }

  /// Creates a matrix according to an array of scalars.
  /// \param m scalar array, size must at least be 16=4*4.
	matrix4( real * m )
	{ set_value(m); }
   
  //// Creates a matrix according to paramter list.
    matrix4( real a00, real a01, real a02, real a03,
	       real a10, real a11, real a12, real a13,
		   real a20, real a21, real a22, real a23,
		   real a30, real a31, real a32, real a33 )
	{
		element(0,0) = a00;
		element(0,1) = a01;
		element(0,2) = a02;
		element(0,3) = a03;
		
		element(1,0) = a10;
		element(1,1) = a11;
		element(1,2) = a12;
		element(1,3) = a13;
		
		element(2,0) = a20;
		element(2,1) = a21;
		element(2,2) = a22;
		element(2,3) = a23;
		
		element(3,0) = a30;
		element(3,1) = a31;
		element(3,2) = a32;
		element(3,3) = a33;
	}
            
    /// Fills an array with the matrix components.
    /// \todo determine column or row major
    /// \param mp Matrix components will be stored here.
    void get_value( real * mp ) const
	{
		int c = 0;
		for(int j=0; j < 4; j++)
			for(int i=0; i < 4; i++)
				mp[c++] = element(i,j);
	}
    
    /// Returns a pointer to the beginning of the matric components as array. 
    const real * get_value() const
	{ return m; }
    
	void set_value( real * mp)
	{
		int c = 0;
		for(int j=0; j < 4; j++)
			for(int i=0; i < 4; i++)
				element(i,j) = mp[c++];
	}
    
	void set_value( real r ) 
	{
		for(int i=0; i < 4; i++)
			for(int j=0; j < 4; j++)
				element(i,j) = r;
	}
   
  /// Creates identity matrix.
    void make_identity()
	{
		element(0,0) = 1.0;
		element(0,1) = 0.0;
		element(0,2) = 0.0; 
		element(0,3) = 0.0;
		
		element(1,0) = 0.0;
		element(1,1) = 1.0; 
		element(1,2) = 0.0;
		element(1,3) = 0.0;
		
		element(2,0) = 0.0;
		element(2,1) = 0.0;
		element(2,2) = 1.0;
		element(2,3) = 0.0;
		
		element(3,0) = 0.0; 
		element(3,1) = 0.0; 
		element(3,2) = 0.0;
		element(3,3) = 1.0;
	}
	
	
    static matrix4 identity()
	{
		static matrix4 mident (
			1.0, 0.0, 0.0, 0.0,
			0.0, 1.0, 0.0, 0.0,
			0.0, 0.0, 1.0, 0.0,
			0.0, 0.0, 0.0, 1.0  );
		return mident;
	}
    
  
    /// Sets the diagonal elements of the upper 3x3 matrix to a fixed value.
    /// \param s New scaling factor.
    void set_scale( real s )
	{
		element(0,0) = s;
		element(1,1) = s;
		element(2,2) = s;
	}
   
    /// Sets the diagonal elements of the upper 3x3 matrix according to the components in s.
    /// \param s New scaling factors.
    void set_scale( const vec3 & s )
	{
		element(0,0) = s[0];
		element(1,1) = s[1];
		element(2,2) = s[2];
	}
    
    /// Sets the 4th column of the matrix (translation) to a fixed vector.
    /// \param t New Translation.
    void set_translate( const vec3 & t )
	{
		element(0,3) = t[0];
		element(1,3) = t[1];
		element(2,3) = t[2];
	}
    
    /// Sets a row of the matrix to a fixed vector.
    /// \param r Row that will be changed
    /// \param t New vector that will be assinged to r-th row.
	void set_row(int r, const vec4 & t)
	{
		element(r,0) = t[0];
		element(r,1) = t[1];
		element(r,2) = t[2];
		element(r,3) = t[3];
	}

	/// Sets a column of the matrix to a fixed vector.
    /// \param c Column that will be changed
    /// \param t New vector that will be assinged to c-th column.
void set_column(int c, const vec4 & t)
	{
		element(0,c) = t[0];
		element(1,c) = t[1];
		element(2,c) = t[2];
		element(3,c) = t[3];
	}

  /// Retrieves a row of the matix.
  /// \param r Row of the matrix that will be retrieved.
  /// \param t Vector that will store the r-th row.
	void get_row(int r, vec4 & t) const
	{
		t[0] = element(r,0);
		t[1] = element(r,1);
		t[2] = element(r,2);
		t[3] = element(r,3);
	}

	/// Retrieves a row of the matix.
  /// \param r Row of the matrix that will be retrieved.
  /// \return Returns vector that will store the r-th row.
vec4 get_row(int r) const
	{
		vec4 v; get_row(r, v);
		return v;
	}

	/// Retrieves a column of the matix.
  /// \param c Column of the matrix that will be retrieved.
  /// \param t Vector that will store the c-th column.
void get_column(int c, vec4 & t) const
	{
		t[0] = element(0,c);
		t[1] = element(1,c);
		t[2] = element(2,c);
		t[3] = element(3,c);
	}

	/// Retrieves a column of the matix.
  /// \param c Column of the matrix that will be retrieved.
  /// \return Returns vector that will store the c-th column.
vec4 get_column(int c) const
	{
		vec4 v; get_column(c, v);
		return v;
	}

    /// Computes and returns the inverse of the matrix. Leaves this matrix unchanged.
    matrix4 inverse() const
	{
		matrix4 minv;
		
		real r1[8], r2[8], r3[8], r4[8];
		real *s[4], *tmprow;
		
		s[0] = &r1[0];
		s[1] = &r2[0];
		s[2] = &r3[0];
		s[3] = &r4[0];
		
		register int i,j,p,jj;
		for(i=0;i<4;i++)
		{
			for(j=0;j<4;j++)
			{
				s[i][j] = element(i,j);
				if(i==j) s[i][j+4] = 1.0;
				else     s[i][j+4] = 0.0;
			}
		}
		real scp[4];
		for(i=0;i<4;i++)
		{
			scp[i] = real(fabs(s[i][0]));
			for(j=1;j<4;j++)
				if(real(fabs(s[i][j])) > scp[i]) scp[i] = real(fabs(s[i][j]));
				if(scp[i] == 0.0) return minv; // singular matrix!
		}
		
		int pivot_to;
		real scp_max;
		for(i=0;i<4;i++)
		{
			// select pivot row
			pivot_to = i;
			scp_max = real(fabs(s[i][i]/scp[i]));
			// find out which row should be on top
			for(p=i+1;p<4;p++)
				if(real(fabs(s[p][i]/scp[p])) > scp_max)
				{ scp_max = real(fabs(s[p][i]/scp[p])); pivot_to = p; }
				// Pivot if necessary
				if(pivot_to != i)
				{
					tmprow = s[i];
					s[i] = s[pivot_to];
					s[pivot_to] = tmprow;
					real tmpscp;
					tmpscp = scp[i];
					scp[i] = scp[pivot_to];
					scp[pivot_to] = tmpscp;
				}
				
				real mji;
				// perform gaussian elimination
				for(j=i+1;j<4;j++)
				{
					mji = s[j][i]/s[i][i];
					s[j][i] = 0.0;
					for(jj=i+1;jj<8;jj++)
						s[j][jj] -= mji*s[i][jj];
				}
		}
		if(s[3][3] == 0.0) return minv; // singular matrix!
		
		//
		// Now we have an upper triangular matrix.
		//
		//  x x x x | y y y y
		//  0 x x x | y y y y 
		//  0 0 x x | y y y y
		//  0 0 0 x | y y y y
		//
		//  we'll back substitute to get the inverse
		//
		//  1 0 0 0 | z z z z
		//  0 1 0 0 | z z z z
		//  0 0 1 0 | z z z z
		//  0 0 0 1 | z z z z 
		//
		
		real mij;
		for(i=3;i>0;i--)
		{
			for(j=i-1;j > -1; j--)
			{
				mij = s[j][i]/s[i][i];
				for(jj=j+1;jj<8;jj++)
					s[j][jj] -= mij*s[i][jj];
			}
		}
		
		for(i=0;i<4;i++)
			for(j=0;j<4;j++)
				minv(i,j) = s[i][j+4] / s[i][i];
			
			return minv;
	}
    
    /// Computes and returns the transpose of the matrix. Leaves this matrix unchanged.
    matrix4 transpose() const
	{
		matrix4 mtrans;
		
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				mtrans(i,j) = element(j,i);		
		return mtrans;
	}
    /// Right-Multiplies with a matrix and stores the result in this matrix M.
    /// \param b matrix 
    /// \return M * b
    matrix4 & mult_right( const matrix4 & b )
	{
		matrix4 mt(*this);
		set_value(real(0));

		for(int i=0; i < 4; i++)
			for(int j=0; j < 4; j++)
				for(int c=0; c < 4; c++)
					element(i,j) += mt(i,c) * b(c,j);
		return *this;
	}    

    /// Left-Multiplies with a matrix and stores the result in this matrix M.
    /// \param b matrix 
    /// \return b * M.
matrix4 & mult_left( const matrix4 & b )
	{
		matrix4 mt(*this);
		set_value(real(0));

		for(int i=0; i < 4; i++)
			for(int j=0; j < 4; j++)
				for(int c=0; c < 4; c++)
					element(i,j) += b(i,c) * mt(c,j);
		return *this;
	}
	
	/// Multiplies matrix with vector: dst = M * src
  /// \param src vector on the right-hand-side
  /// \param dst stores result
    void mult_matrix_vec( const vec3 &src, vec3 &dst ) const
	{
		real w = (
			src[0] * element(3,0) +
			src[1] * element(3,1) + 
			src[2] * element(3,2) +
			element(3,3)          );
        
        assert(w != GLH_ZERO);

        dst[0]  = (
			src[0] * element(0,0) +
			src[1] * element(0,1) +
			src[2] * element(0,2) +
			element(0,3)          ) / w;
		dst[1]  = (
			src[0] * element(1,0) +
			src[1] * element(1,1) +
			src[2] * element(1,2) +
			element(1,3)          ) / w;
		dst[2]  = (
			src[0] * element(2,0) +
			src[1] * element(2,1) + 
			src[2] * element(2,2) +
			element(2,3)          ) / w;
	}
    
	/// Multiplies matrix with vector: dst = M * src
  /// \param src_and_dst vector on the right-hand-side and stores result
void mult_matrix_vec( vec3 & src_and_dst) const
	{ mult_matrix_vec(vec3(src_and_dst), src_and_dst); }


	  /// Multiplies vector with matrix : dst = src * M 
  /// \param src vector on the left-hand-side
  /// \param dst stores result
  void mult_vec_matrix( const vec3 &src, vec3 &dst ) const
	{
		real w = (
			src[0] * element(0,3) +
			src[1] * element(1,3) +
			src[2] * element(2,3) +
			element(3,3)          );
        
        assert(w != GLH_ZERO);

		dst[0]  = (
			src[0] * element(0,0) +
			src[1] * element(1,0) + 
			src[2] * element(2,0) + 
			element(3,0)          ) / w;
		dst[1]  = (
			src[0] * element(0,1) +
			src[1] * element(1,1) +
			src[2] * element(2,1) +
			element(3,1)          ) / w;
		dst[2]  = (
			src[0] * element(0,2) +
			src[1] * element(1,2) +
			src[2] * element(2,2) +
			element(3,2)          ) / w;
	}
        

  /// Multiplies vector with matrix : dst = src * M 
  /// \param src_and_dst vector on the left-hand-side and stores result
	void mult_vec_matrix( vec3 & src_and_dst) const
	{ mult_vec_matrix(vec3(src_and_dst), src_and_dst); }

	///  dst = M * src
    void mult_matrix_vec( const vec4 &src, vec4 &dst ) const
	{
        dst[0]  = (
			src[0] * element(0,0) +
			src[1] * element(0,1) +
			src[2] * element(0,2) +
			src[3] * element(0,3));
		dst[1]  = (
			src[0] * element(1,0) +
			src[1] * element(1,1) +
			src[2] * element(1,2) +
			src[3] * element(1,3));
		dst[2]  = (
			src[0] * element(2,0) +
			src[1] * element(2,1) + 
			src[2] * element(2,2) +
			src[3] * element(2,3));
		dst[3] = (
			src[0] * element(3,0) +
			src[1] * element(3,1) + 
			src[2] * element(3,2) +
			src[3] * element(3,3));
	}
    
	///  dst = M * src
	void mult_matrix_vec( vec4 & src_and_dst) const
	{ mult_matrix_vec(vec4(src_and_dst), src_and_dst); }


    /// dst = src * M
    void mult_vec_matrix( const vec4 &src, vec4 &dst ) const
	{
		dst[0]  = (
			src[0] * element(0,0) +
			src[1] * element(1,0) + 
			src[2] * element(2,0) + 
			src[3] * element(3,0));
		dst[1]  = (
			src[0] * element(0,1) +
			src[1] * element(1,1) +
			src[2] * element(2,1) +
			src[3] * element(3,1));
		dst[2]  = (
			src[0] * element(0,2) +
			src[1] * element(1,2) +
			src[2] * element(2,2) +
			src[3] * element(3,2));
		dst[3] = (
			src[0] * element(0,3) +
			src[1] * element(1,3) +
			src[2] * element(2,3) +
			src[3] * element(3,3));
	}
        

  /// dst = src * M
	void mult_vec_matrix( vec4 & src_and_dst) const
	{ mult_vec_matrix(vec4(src_and_dst), src_and_dst); }

    
    /// dst = M * src
    /// \todo check this function
    void mult_matrix_dir( const vec3 &src, vec3 &dst ) const
	{
		dst[0]  = (
			src[0] * element(0,0) +
			src[1] * element(0,1) +
			src[2] * element(0,2) ) ;
		dst[1]  = ( 
			src[0] * element(1,0) +
			src[1] * element(1,1) +
			src[2] * element(1,2) ) ;
		dst[2]  = ( 
			src[0] * element(2,0) +
			src[1] * element(2,1) + 
			src[2] * element(2,2) ) ;
	}
        

  /// dst = M * src
  /// \todo check this function
	void mult_matrix_dir( vec3 & src_and_dst) const
	{ mult_matrix_dir(vec3(src_and_dst), src_and_dst); }


	/// dst = src * M
    void mult_dir_matrix( const vec3 &src, vec3 &dst ) const
	{
		dst[0]  = ( 
			src[0] * element(0,0) +
			src[1] * element(1,0) +
			src[2] * element(2,0) ) ;
		dst[1]  = ( 
			src[0] * element(0,1) +
			src[1] * element(1,1) +
			src[2] * element(2,1) ) ;
		dst[2]  = (
			src[0] * element(0,2) +
			src[1] * element(1,2) + 
			src[2] * element(2,2) ) ;
	}
    
  /// dst = src * M  
	void mult_dir_matrix( vec3 & src_and_dst) const
	{ mult_dir_matrix(vec3(src_and_dst), src_and_dst); }


    /// Element access (read/write)
    /// \param row read/write that row
    /// \param col read/write that column
    real & operator () (int row, int col)
    { return element(row,col); }

    /// Element access (read-only)
    /// \param row rea that row
    /// \param col read that column
    const real & operator () (int row, int col) const
    { return element(row,col); }

	real & element (int row, int col)
    { return m[row | (col<<2)]; }

    const real & element (int row, int col) const
    { return m[row | (col<<2)]; }

    /// cast to T array
    inline operator real*() { return m; }

		/// cast to const T array
		inline operator const real*() const { return m; }

   matrix4 & operator *= ( const matrix4 & mat )
	{
		mult_right( mat );
		return *this;
	}
    
    matrix4 & operator *= ( const real & r )
	{
		for (int i = 0; i < 4; ++i)
        {
            element(0,i) *= r;
            element(1,i) *= r;
            element(2,i) *= r;
            element(3,i) *= r;
        }
		return *this;
	}

    matrix4 & operator += ( const matrix4 & mat )
	{
		for (int i = 0; i < 4; ++i)
        {
            element(0,i) += mat.element(0,i);
            element(1,i) += mat.element(1,i);
            element(2,i) += mat.element(2,i);
            element(3,i) += mat.element(3,i);
        }
		return *this;
	}

    friend const vec3 operator * ( const matrix4 & m,	const vec3 & right);
    friend matrix4 operator * ( const matrix4 & m1,	const matrix4 & m2 );
    friend bool operator == ( const matrix4 & m1, const matrix4 & m2 );
    friend bool operator != ( const matrix4 & m1, const matrix4 & m2 );
  

  template <typename Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar & m;
  }

  protected:
	  real m[16];
  };
  
  inline const vec3 operator * ( const matrix4 & m,	const vec3 & right)
  {
    vec3 dst;
    m.mult_matrix_vec(right,dst);
    return dst;
  }
  
  inline  
  matrix4 operator * ( const matrix4 & m1, const matrix4 & m2 )
  {
	  matrix4 product;
	  
	  product = m1;
	  product.mult_right(m2);
	  
	  return product;
  }
  
  inline
  bool operator ==( const matrix4 &m1, const matrix4 &m2 )
  {
	  return ( 
		  m1(0,0) == m2(0,0) &&
		  m1(0,1) == m2(0,1) &&
		  m1(0,2) == m2(0,2) &&
		  m1(0,3) == m2(0,3) &&
		  m1(1,0) == m2(1,0) &&
		  m1(1,1) == m2(1,1) &&
		  m1(1,2) == m2(1,2) &&
		  m1(1,3) == m2(1,3) &&
		  m1(2,0) == m2(2,0) &&
		  m1(2,1) == m2(2,1) &&
		  m1(2,2) == m2(2,2) &&
		  m1(2,3) == m2(2,3) &&
		  m1(3,0) == m2(3,0) &&
		  m1(3,1) == m2(3,1) &&
		  m1(3,2) == m2(3,2) &&
		  m1(3,3) == m2(3,3) );
  }
  
  inline
  bool operator != ( const matrix4 & m1, const matrix4 & m2 )
  { return !( m1 == m2 ); }  












    /// Quaternion 
    class quaternion
    {
    public:
   
    /// Creates identity quaternion
    quaternion()
    {
        *this = identity();
    }

    /// Constructs quaternion from a float array 
    quaternion( const real v[4] )
    {
        set_value( v );
    }

    /// Constructs quaternion from float values
    quaternion( real q0, real q1, real q2, real q3 )
    {
        set_value( q0, q1, q2, q3 );
    }

    /// Constructs quaternion from the upper 3x3 rotation part of a matrix
    /// \param m rotation matrix
    quaternion( const matrix4 & m )
    {
        set_value( m );
    }

    /// Constructs quaternion as a rotation around an axis 
    /// \param axis rotation axis
    /// \param rotation angle around axis
    quaternion( const vec3 &axis, real radians )
    {
        set_value( axis, radians );
    }

    /// Constructs quaternion such that it represents the rotation that
    /// transforms one direction into another 
    /// \param rotate_from source direction
    /// \param rotate_to target direction
    quaternion( const vec3 &rotate_from, const vec3 &rotate_to )
    {
        set_value( rotate_from, rotate_to );
    }

    /// Constructs quaternion such that is orients a view direction with an
    /// vector pointing upwards to a new view and up direction.
    /// \param from_look source view direction
    /// \param from_up source up vector
    /// \param to_look target view direction
    /// \param to_up target up direction
    quaternion( const vec3 & from_look, const vec3 & from_up,
		      const vec3 & to_look, const vec3& to_up)
    {
	    set_value(from_look, from_up, to_look, to_up);
    }

    /// Returns pointer to quaternion components
    const real * get_value() const
    {
        return  &q[0];
    }

    /// Retrieve quaternion components
    void get_value( real &q0, real &q1, real &q2, real &q3 ) const
    {
        q0 = q[0];
        q1 = q[1];
        q2 = q[2];
        q3 = q[3];
    }

    quaternion & set_value( real q0, real q1, real q2, real q3 )
    {
        q[0] = q0;
        q[1] = q1;
        q[2] = q2;
        q[3] = q3;
        counter = 0;
        return *this;
    }

    /// Retrieves the rotation axis and angle represented by the quaternion
    /// \param axis rotation axis
    /// \param radians angle in radians
    void get_value( vec3 &axis, real &radians ) const
    {
        radians = real(acos( q[3] ) * GLH_TWO);
        if ( radians == GLH_ZERO )
            axis = vec3( 0.0, 0.0, 1.0 );
        else
        {
            axis[0] = q[0];
            axis[1] = q[1];
            axis[2] = q[2];
            axis.normalize();
        }
    }

    /// Retrieves the rotation matrix represented by the quaternion
    /// \param m rotation matrix
    void get_value( matrix4 & m ) const
    {
        real s, xs, ys, zs, wx, wy, wz, xx, xy, xz, yy, yz, zz;

        real norm = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];

        s = (equivalent(norm,GLH_ZERO)) ? GLH_ZERO : ( GLH_TWO / norm );

        xs = q[0] * s;
        ys = q[1] * s;
        zs = q[2] * s;

        wx = q[3] * xs;
        wy = q[3] * ys;
        wz = q[3] * zs;

        xx = q[0] * xs;
        xy = q[0] * ys;
        xz = q[0] * zs;

        yy = q[1] * ys;
        yz = q[1] * zs;
        zz = q[2] * zs;

        m(0,0) = real( GLH_ONE - ( yy + zz ));
        m(1,0) = real ( xy + wz );
        m(2,0) = real ( xz - wy );

        m(0,1) = real ( xy - wz );
        m(1,1) = real ( GLH_ONE - ( xx + zz ));
        m(2,1) = real ( yz + wx );

        m(0,2) = real ( xz + wy );
        m(1,2) = real ( yz - wx );
        m(2,2) = real ( GLH_ONE - ( xx + yy ));

        m(3,0) = m(3,1) = m(3,2) = m(0,3) = m(1,3) = m(2,3) = GLH_ZERO;
        m(3,3) = GLH_ONE;
    }

    quaternion & set_value( const real * qp )
    {
        memcpy(q,qp,sizeof(real) * 4);

        counter = 0;
        return *this;
    }

    quaternion & set_value( const matrix4 & m )
    {
        real tr, s;
        int i, j, k;
        const int nxt[3] = { 1, 2, 0 };

        tr = m(0,0) + m(1,1) + m(2,2);

        if ( tr > GLH_ZERO )
        {
            s = real(sqrt( tr + m(3,3) ));
            q[3] = real ( s * 0.5 );
            s = real(0.5) / s;

            q[0] = real ( ( m(1,2) - m(2,1) ) * s );
            q[1] = real ( ( m(2,0) - m(0,2) ) * s );
            q[2] = real ( ( m(0,1) - m(1,0) ) * s );
        }
        else
        {
            i = 0;
            if ( m(1,1) > m(0,0) )
              i = 1;

            if ( m(2,2) > m(i,i) )
              i = 2;

            j = nxt[i];
            k = nxt[j];

            s = real(sqrt( ( m(i,j) - ( m(j,j) + m(k,k) )) + GLH_ONE ));

            q[i] = real ( s * 0.5 );
            s = real(0.5 / s);

            q[3] = real ( ( m(j,k) - m(k,j) ) * s );
            q[j] = real ( ( m(i,j) + m(j,i) ) * s );
            q[k] = real ( ( m(i,k) + m(k,i) ) * s );
        }

        counter = 0;
        return *this;
    }

    quaternion & set_value( const vec3 &axis, real theta )
    {
        real sqnorm = axis.square_norm();

        if (sqnorm <= GLH_EPSILON)
        {
            // axis too small.
            x = y = z = 0.0;
            w = 1.0;
        } 
        else 
        {
            theta *= real(0.5);
            real sin_theta = real(sin(theta));

            if (!equivalent(sqnorm,GLH_ONE)) 
              sin_theta /= real(sqrt(sqnorm));
            x = sin_theta * axis[0];
            y = sin_theta * axis[1];
            z = sin_theta * axis[2];
            w = real(cos(theta));
        }
        return *this;
    }

    quaternion & set_value( const vec3 & rotateFrom, const vec3 & rotateTo )
    {
        vec3 p1, p2;
        real alpha;

        p1 = rotateFrom; 
        p1.normalize();
        p2 = rotateTo;  
        p2.normalize();

        alpha = p1.dot(p2);

        if(equivalent(alpha,GLH_ONE))
        { 
            *this = identity(); 
            return *this; 
        }

        // ensures that the anti-parallel case leads to a positive dot
        if(equivalent(alpha,-GLH_ONE))
        {
            vec3 v;

            if(p1[0] != p1[1] || p1[0] != p1[2])
    	        v = vec3(p1[1], p1[2], p1[0]);
            else
    	        v = vec3(-p1[0], p1[1], p1[2]);

            v -= p1 * p1.dot(v);
            v.normalize();

            set_value(v, GLH_PI);
            return *this;
        }

        p1 = p1.cross(p2);  
        p1.normalize();
        set_value(p1,real(acos(alpha)));

        counter = 0;
        return *this;
    }

    quaternion & set_value( const vec3 & from_look, const vec3 & from_up,
		      const vec3 & to_look, const vec3 & to_up)
    {
	    quaternion r_look = quaternion(from_look, to_look);
	    
	    vec3 rotated_from_up(from_up);
	    r_look.mult_vec(rotated_from_up);
	    
	    quaternion r_twist = quaternion(rotated_from_up, to_up);
	    
	    *this = r_twist;
	    *this *= r_look;
	    return *this;
    }

    /// Concatenates two quaternions
    quaternion & operator *= ( const quaternion & qr )
    {
        quaternion ql(*this);
   
        w = ql.w * qr.w - ql.x * qr.x - ql.y * qr.y - ql.z * qr.z;
        x = ql.w * qr.x + ql.x * qr.w + ql.y * qr.z - ql.z * qr.y;
        y = ql.w * qr.y + ql.y * qr.w + ql.z * qr.x - ql.x * qr.z;
        z = ql.w * qr.z + ql.z * qr.w + ql.x * qr.y - ql.y * qr.x;

        counter += qr.counter;
        counter++;
        counter_normalize();
        return *this;
    }

    /// Normalizes the quaternion such that its length is one.
    void normalize()
    {
        real rnorm = GLH_ONE / real(sqrt(w * w + x * x + y * y + z * z));
        if (equivalent(rnorm, GLH_ZERO))
            return;
        x *= rnorm;
        y *= rnorm;
        z *= rnorm;
        w *= rnorm;
        counter = 0;
    }

    /// component wise equal operation
    friend bool operator == ( const quaternion & q1, const quaternion & q2 );      

    /// component wise not-equal operation
    friend bool operator != ( const quaternion & q1, const quaternion & q2 );

    friend quaternion operator * ( const quaternion & q1, const quaternion & q2 );

    bool equals( const quaternion & r, real /*tolerance */) const
    {
        real t;

        t = (
			(q[0]-r.q[0])*(q[0]-r.q[0]) +
            (q[1]-r.q[1])*(q[1]-r.q[1]) +
            (q[2]-r.q[2])*(q[2]-r.q[2]) +
            (q[3]-r.q[3])*(q[3]-r.q[3]) );
        if(t > GLH_EPSILON) 
            return false;
        return 1;
    }

    /// Conjugates (negates first three components) of the quaternion
    /// \return the conjugates quaternion
    quaternion & conjugate()
    {
        q[0] *= -GLH_ONE;
        q[1] *= -GLH_ONE;
        q[2] *= -GLH_ONE;
        return *this;
    }

    /// Same as conjugate(). Inverts first three components.
    quaternion & invert()
    {
        return conjugate();
    }

    /// Returns the conjugated quaternion without changing this quaternion
    quaternion inverse() const
    {
        quaternion r = *this;
        return r.invert();
    }

    /// Quaternion multiplication with cartesian vector
    /// v' = q*v*q(star)
    /// \param src source vector
    /// \param dst target vector which stores result
    void mult_vec( const vec3 &src, vec3 &dst ) const
    {
        real v_coef = w * w - x * x - y * y - z * z;                     
        real u_coef = GLH_TWO * (src[0] * x + src[1] * y + src[2] * z);  
        real c_coef = GLH_TWO * w;                                       

        dst[0] = v_coef * src[0] + u_coef * x + c_coef * (y * src[2] - z * src[1]);
        dst[1] = v_coef * src[1] + u_coef * y + c_coef * (z * src[0] - x * src[2]);
        dst[2] = v_coef * src[2] + u_coef * z + c_coef * (x * src[1] - y * src[0]);
    }

    /// Quaternion multiplication with cartesian vector. Stores the result of
    /// the computation in the source vector.  v' = q*v*q(star)
    /// \param src_and_dst source vector and target vector 
    void mult_vec( vec3 & src_and_dst) const
    {
        mult_vec(vec3(src_and_dst), src_and_dst);
    }

    /// Scales the rotation angle by a factor. Leaves the rotation axis unchanged.
    /// \param scaleFactor scaling factor.
    void scale_angle( real scaleFactor )
    {
        vec3 axis;
        real radians;

        get_value(axis, radians);
        radians *= scaleFactor;
        set_value(axis, radians);
    }

    /// Performs a spherical interpolation between two quaternions. 
    /// res = p * (1-alpha) + q * alpha
    /// \param p first quaternion
    /// \param q second quaternion
    /// \param alpha interpolation coefficient
    /// \return interpolated quaternion
    static quaternion slerp( const quaternion & p, const quaternion & q, real alpha )
    {
        quaternion r;

        real cos_omega = p.x * q.x + p.y * q.y + p.z * q.z + p.w * q.w;
        // if B is on opposite hemisphere from A, use -B instead
      
        int bflip;
        if ( ( bflip = (cos_omega < GLH_ZERO)) )
            cos_omega = -cos_omega;

        // complementary interpolation parameter
        real beta = GLH_ONE - alpha;     

        if(cos_omega >= GLH_ONE - GLH_EPSILON)
            return p;

        real omega = real(acos(cos_omega));
        real one_over_sin_omega = GLH_ONE / real(sin(omega));

        beta    = real(sin(omega*beta)  * one_over_sin_omega);
        alpha   = real(sin(omega*alpha) * one_over_sin_omega);

        if (bflip)
            alpha = -alpha;

        r.x = beta * p.q[0]+ alpha * q.q[0];
        r.y = beta * p.q[1]+ alpha * q.q[1];
        r.z = beta * p.q[2]+ alpha * q.q[2];
        r.w = beta * p.q[3]+ alpha * q.q[3];
        return r;
    }

    /// creates identity quaternion
    static quaternion identity()
    {
        static quaternion ident( vec3( 0.0, 0.0, 0.0 ), GLH_ONE );
        return ident;
    }

    /// Component wise access access (read/write)
    real & operator []( int i )
    {
        assert(i < 4);
        return q[i];
    }

    /// Component wise access access (read-only)
    const real & operator []( int i ) const
    {
        assert(i < 4);
        return q[i];
    }

    protected:

        void counter_normalize()
        {
            if (counter > GLH_QUATERNION_NORMALIZATION_THRESHOLD)
                normalize();
        }

        union 
        {
            struct 
            {
                real q[4];
            };
            struct 
            {
                real x;
                real y;
                real z;
                real w;
            };
        };

        // renormalization counter
        unsigned char counter;
    };

    inline
    bool operator == ( const quaternion & q1, const quaternion & q2 )
    {
        return (equivalent(q1.x, q2.x) &&
		        equivalent(q1.y, q2.y) &&
		        equivalent(q1.z, q2.z) &&
		        equivalent(q1.w, q2.w) );
    }

    inline
    bool operator != ( const quaternion & q1, const quaternion & q2 )
    { 
        return ! ( q1 == q2 ); 
    }

    inline
    quaternion operator * ( const quaternion & q1, const quaternion & q2 )
    {	
        quaternion r(q1); 
        r *= q2; 
        return r; 
    }
  
      
    





  
  /// Plane in 3D 
  class plane
  {
  public:
	  
	  plane()
      {
		  planedistance = 0.0;
		  planenormal.set_value( 0.0, 0.0, 1.0 );
      }
	  
	  
	  plane( const vec3 &p0, const vec3 &p1, const vec3 &p2 )
      {
		  vec3 v0 = p1 - p0;
		  vec3 v1 = p2 - p0;
		  planenormal = v0.cross(v1);  
		  planenormal.normalize();
		  planedistance = p0.dot(planenormal);
      }
	  
	  plane( const vec3 &normal, real distance )
      {
		  planedistance = distance;
		  planenormal = normal;
		  planenormal.normalize();
      }
	  
	  plane( const vec3 &normal, const vec3 &point )
      {
		  planenormal = normal;
		  planenormal.normalize();
		  planedistance = point.dot(planenormal);
      }
	  
	  void offset( real d )
      {
		  planedistance += d;
      }
	  
	  bool intersect( const line &l, vec3 &intersection ) const
      {
		  vec3 pos, dir;
		  vec3 pn = planenormal;
		  real pd = planedistance;
		  
		  pos = l.get_position();
		  dir = l.get_direction();
		  
		  if(dir.dot(pn) == 0.0) return 0;
		  pos -= pn*pd;
		  // now we're talking about a plane passing through the origin
		  if(pos.dot(pn) < 0.0) pn.negate();
		  if(dir.dot(pn) > 0.0) dir.negate();
		  vec3 ppos = pn * pos.dot(pn);
		  pos = dir*(ppos.length()/dir.dot(-pn));
		  intersection = l.get_position();
		  intersection += pos;
		  return 1;
      }
	  void transform( const matrix4 &matrix )
      {
		  matrix4 invtr = matrix.inverse();
		  invtr = invtr.transpose();
		  
		  vec3 pntOnplane = planenormal * planedistance;
		  vec3 newPntOnplane;
		  vec3 newnormal;
		  
		  invtr.mult_dir_matrix(planenormal, newnormal);
		  matrix.mult_vec_matrix(pntOnplane, newPntOnplane);
		  
		  newnormal.normalize();
		  planenormal = newnormal;
		  planedistance = newPntOnplane.dot(planenormal);
      }
	  
	  bool is_in_half_space( const vec3 &point ) const
      {
		  
		  if(( point.dot(planenormal) - planedistance) < 0.0)
			  return 0;
		  return 1;
      }
	  
	  
	  real distance( const vec3 & point ) const 
      {
		  return planenormal.dot(point - planenormal*planedistance);
      }
	  
	  const vec3 &get_normal() const
      {
		  return planenormal;
      }
	  
	  
	  real get_distance_from_origin() const
      {
		  return planedistance;
      }
	  
	  
	  friend bool operator == ( const plane & p1, const plane & p2 );
	  
	  
	  friend bool operator != ( const plane & p1, const plane & p2 );
	  
  //protected:
	  vec3 planenormal;
	  real planedistance;
  };
  
  inline
  bool operator == (const plane & p1, const plane & p2 )
  {
	  return (  p1.planedistance == p2.planedistance && p1.planenormal == p2.planenormal);
  }
  
  inline
  bool operator != ( const plane & p1, const plane & p2 )
  { return  ! (p1 == p2); }
  
  

  //} // "ns_##GLH_REAL"

  // make common typedefs...
#ifdef GLH_REAL_IS_FLOAT
  /// 2D float Vector
/*  typedef GLH_REAL_NAMESPACE::vec2 vec2f;
  /// 3D float Vector
  typedef GLH_REAL_NAMESPACE::vec3 vec3f;
  /// 4D float Vector
  typedef GLH_REAL_NAMESPACE::vec4 vec4f;
  /// Quaternion
  typedef GLH_REAL_NAMESPACE::quaternion quaternionf;
  typedef GLH_REAL_NAMESPACE::quaternion rotationf;
  /// Line in 3D
  typedef GLH_REAL_NAMESPACE::line linef;
  /// Plane in 3D
  typedef GLH_REAL_NAMESPACE::plane planef;
  /// 4x4 float Matrix
  typedef GLH_REAL_NAMESPACE::matrix4 matrix4f;
  */
#endif

  


}  // namespace glh

/// output a vector by printing its space-separated compontens
/*inline std::ostream& 
operator<<(std::ostream& os, const glh::vec3i& v) {
	for(int i=0; i<3-1; i++) os << v[i] << " ";
	os << v[3-1];
	return os;
}

*/

/// print to stream for vec<N,T>
template<int N, class T>
inline std::ostream& 
operator<<(std::ostream& os, const glh::vec<N, T>& v) {
	for(int i=0; i<N-1; i++) os << static_cast<T>(v[i]) << " ";
	os << static_cast<T>(v[N-1]);
	return os;
}

/// print to sream for matrix4
inline std::ostream& 
operator<<(std::ostream& os, const glh::matrix4& m) {
	os << "\n";
  for(int i=0; i<4; i++) 
  {
	  for(int j=0; j<4; j++) 
      os << m(i,j) << " ";
    os << "\n";
  }
	return os;
}

#define M_PI 3.141592654

/// print to sream for quaternion (as axis and angle)
inline std::ostream& 
operator<<(std::ostream& os, const glh::quaternion& q) {
  glh::vec3 axis;
  float angle_rad;
  q.get_value(axis,angle_rad);

  const float angle_deg = (angle_rad/M_PI) * 180.0f;
  os << axis << ", " << angle_deg;
	
  return os;
}





#endif

