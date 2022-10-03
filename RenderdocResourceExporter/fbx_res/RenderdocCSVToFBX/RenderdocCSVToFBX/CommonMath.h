#include <cassert>
#include <algorithm>
#include <ccomplex>
#ifndef _COMMONMATH_H_INCLUDE
#define _COMMONMATH_H_INCLUDE


#define FLT_HALFPI  (3.1415926f*0.5)
#define FLT_PI      3.1415926f
#define FLT_2PI		6.2831852f

#define FLT_DTOR( degree ) ((degree) * (FLT_PI / 180.0f))
#define FLT_RTOD( radian ) ((radian) * (180.0f / FLT_PI))

#define FLT_SQRT_OF_2 (1.4142135f)


#ifdef ANDROID
#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif
#endif
int		Argmin(float a[],int n);
float  Sqr(float a); 
float  Roundf(float a,float precision);
float  Interpolate(const float &f0,const float &f1,float alpha) ;

inline bool	IsNan(float x)
{ 
	return x != x;
}

template <class T>
void Swap(T &a,T &b)
{
	T tmp = a;
	a=b;
	b=tmp;
}



template <class T>
T Max(const T &a,const T &b)
{
	return (a>b)?a:b;
}

template <class T>
T Min(const T &a,const T &b) 
{
	return (a<b)?a:b;
}
//----------------------------------
class  VectorI2
{
public:
	int x;
	int y;

	VectorI2(){x = 0;y = 0;}

	VectorI2(int _x,int _y){x = _x;y = _y;}

	void	Set(int _x,int _y){x = _x;y = _y;}

	bool	operator ==(VectorI2& vec)
	{
		return x == vec.x && y == vec.y; 
	}

	bool	operator ==(const VectorI2& vec)
	{
		return x == vec.x && y == vec.y; 
	}

	VectorI2&	operator -=(VectorI2& vec)
	{
		x -= vec.x;
		y -= vec.y; 
		return *this;
	}

	VectorI2&	operator -=(const VectorI2& vec)
	{
		x -= vec.x;
		y -= vec.y; 
		return *this;
	}

	VectorI2&	operator +=(VectorI2& vec)
	{
		x += vec.x;
		y += vec.y; 
		return *this;
	}

	VectorI2&	operator +=(const VectorI2& vec)
	{
		x += vec.x;
		y += vec.y; 
		return *this;
	}

	VectorI2	operator +(VectorI2& vec)
	{ 
		return VectorI2(x +vec.x,y + vec.y);
	}

	VectorI2	operator +(const VectorI2& vec)
	{
		return VectorI2(x +vec.x,y + vec.y);
	}

	VectorI2	operator -(VectorI2& vec)
	{ 
		return VectorI2(x -vec.x,y - vec.y);
	}

	VectorI2	operator -(const VectorI2& vec)
	{
		return VectorI2(x -vec.x,y - vec.y);
	}

	bool	operator <(const VectorI2& vec)
	{
		return (x*x + y*y) < (vec.x*vec.x + vec.y *vec.y);
	}


	const int& operator[](int i) const {return (&x)[i];}
	int& operator[](int i) {return (&x)[i];}  
};

inline bool
	operator <(const VectorI2& vec0,const VectorI2& vec1)
{
	return (vec0.x*vec0.x + vec0.y *vec0.y) < (vec1.x*vec1.x + vec1.y *vec1.y);  
}

class  VectorI3  
{
public:
	int x,y,z;
	VectorI3(){};
	VectorI3(int _x,int _y, int _z){x=_x;y=_y;z=_z;}
	const int& operator[](int i) const {return (&x)[i];}
	int& operator[](int i) {return (&x)[i];}

	void	Set(int _x,int _y,int _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}
	bool	operator ==(VectorI3& vec)
	{
		return x == vec.x && y == vec.y && z == vec.z; 
	}

	bool	operator ==(const VectorI3& vec)
	{
		return x == vec.x && y == vec.y && z == vec.z; 
	}
};


//-------- 2D --------

class  Vector2
{
public:
	float x,y;
	Vector2(){x=0;y=0;};
	Vector2( const float *_in ){x = _in[0];y = _in[1];}
	Vector2(float _x,float _y){x=_x;y=_y;}
	float& operator[](int i) {assert(i>=0&&i<2);return ((float*)this)[i];}
	const float& operator[](int i) const {assert(i>=0&&i<2);return ((float*)this)[i];}

#ifdef _USING_DXX
	operator       struct D3DXVECTOR2* ()       { return (struct D3DXVECTOR2*) this;}
	operator const struct D3DXVECTOR2* () const { return (struct D3DXVECTOR2*) this;}
#endif
	Vector2&	operator =(const VectorI2& vec)
	{
		x = (float)vec.x;
		y = (float)vec.y;
		return *this;
	}

	Vector2 operator-( const Vector2& b ) const{return Vector2(x-b.x,y-b.y);}
	Vector2 operator+( const Vector2& b ) const{return Vector2(x+b.x,y+b.y);}

	__inline Vector2  operator-()		const                  
	{
		return Vector2( -x, -y );
	}
 

	__inline Vector2&  operator+=(const Vector2& b )
	{
		x += b.x;
		y += b.y;
 
		return *this;
	}
	__inline Vector2&  operator-=(const Vector2& b )
	{
		x -= b.x;
		y -= b.y;

		return *this;
	}

	__inline Vector2&  operator*=(float s )
	{
		x *= s;
		y *= s; 
		return *this;
	}


	__inline Vector2&  operator/=(float s )
	{
		float sinv = 1.0f / s;
		x *= sinv;
		y *= sinv; 
		return *this;
	}

	__inline Vector2  operator*(float s )      const
	{
		return Vector2( x*s, y*s);
	}

	int operator!=(const Vector2 &b ) 
	{
		return (x!=b.x || y!=b.y ); 
	}
	int operator!=(const Vector2 &b )  const
	{
		return (x!=b.x || y!=b.y); 
	}



	float  Magnitude()	const
	{
		return std::sqrt(Sqr(x) + Sqr( y));
	}

	float MagnitudeSq()		const
	{
		return Sqr(x) + Sqr(y);
	}

	void	SetNormalize()
	{
		float d = Magnitude();
		if (d==0)
		{
			printf("Cant normalize ZERO vector\n");
			assert(0);// yes this could go here
			d=0.1f;
		}
		d = 1/d;
		x*= d;
		y*= d;
	}

	Vector2  Normalize();
	Vector2  Safenormalize();

	__inline float   Dot(const Vector2& b )    const
	{
		return  x*b.x +  y*b.y; 
	}

	 

};

inline Vector2 Vector2::Normalize()
{
	// this routine, normalize, is ok, provided magnitude works!!
	float d = Magnitude();
	if (d==0)
	{
		printf("Cant normalize ZERO vector\n");
		assert(0);// yes this could go here
		d=0.1f;
	}
	d = 1/d;
	return Vector2(x*d,y*d);
}

inline Vector2 Vector2::Safenormalize()
{
	if(Magnitude()<=0.0f)
	{
		return Vector2(1,0);
	}
	return Normalize();
}

//--------- 3D ---------

class  Vector3 // 3D优化可能包括去除构造函数,而成为结构,来提高速度
{
public:
	float x,y,z;
	Vector3(){x=0;y=0;z=0;};
	Vector3( const float *_in ){x = _in[0];y = _in[1]; z= _in[2];}
	Vector3(float _x,float _y,float _z){x=_x;y=_y;z=_z;};
	//operator float *() { return &x;};
	float& operator[](int i) {assert(i>=0&&i<3);return ((float*)this)[i];}
	const float& operator[](int i) const {assert(i>=0&&i<3);return ((float*)this)[i];}

	void	ResetZero()
	{
		x=0;y=0;z=0;
	}
#ifdef _USING_DXX
	operator       struct D3DXVECTOR3* ()       { return (struct D3DXVECTOR3*) this;}
	operator const struct D3DXVECTOR3* () const { return (struct D3DXVECTOR3*) this;}

	operator       D3DVECTOR* ()       { return (D3DVECTOR*) this;}
	operator const D3DVECTOR* () const { return (D3DVECTOR*) this;}

#endif

	__inline Vector3  operator+(const Vector3& b )	const
	{
		return Vector3(x+b.x, y+b.y, z+b.z); 
	}
 

	__inline Vector3  operator-(const Vector3& b )	const
	{
		return Vector3( x-b.x, y-b.y, z-b.z );
	}

 

	__inline Vector3  operator-()		const                  
	{
		return Vector3( -x, -y, -z );
	}

 

 
	__inline Vector3  operator*(float s )      const
	{
		return Vector3( x*s, y*s, z*s );
	}
 

	__inline Vector3  operator/(float s )	const
	{ 
		return Vector3(x*(1.0f/s),y*(1.0f/s),z*(1.0f/s)); 
	}
 
	__inline float   Dot(const Vector3& b )    const
	{
		return  x*b.x +  y*b.y + z*b.z; 
	}
 
 
	Vector3   CMul(const Vector3& b )		const
	{
		return    Vector3(x*b.x, y*b.y, z*b.z); 
	}

 


	__inline Vector3  Cross(const Vector3& b )	const
	{
		return Vector3( y*b.z - z*b.y,
			z*b.x - x*b.z,
			x*b.y - y*b.x );
	}


	__inline Vector3&  operator+=(const Vector3& b )
	{
		x += b.x;
		y += b.y;
		z += b.z;
		return *this;
	}


	__inline Vector3&  operator-=(const Vector3& b )
	{
		x -= b.x;
		y -= b.y;
		z -= b.z;
		return *this;
	}


	__inline Vector3&  operator*=(float s )
	{
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}


	__inline Vector3&  operator/=(float s )
	{
		float sinv = 1.0f / s;
		x *= sinv;
		y *= sinv;
		z *= sinv;
		return *this;
	}

	int operator==(const Vector3 &b ) 
	{
		return (x==b.x && y==b.y && z==b.z); 
	}

	int operator==(const Vector3 &b )	const
	{
		return (x==b.x && y==b.y && z==b.z); 
	}

	int operator!=(const Vector3 &b ) 
	{
		return (x!=b.x || y!=b.y || z!=b.z); 
	}
	int operator!=(const Vector3 &b )  const
	{
		return (x!=b.x || y!=b.y || z!=b.z); 
	}


	float  Magnitude()	const
	{
		return std::sqrt(Sqr(x) + Sqr( y)+ Sqr(z));
	}
 

	float MagnitudeSq()	const
	{
		return Sqr(x) + Sqr(y) + Sqr(z);
	}
	void	SetNormalize()
	{
		float d = Magnitude();
		if (d==0)
		{
			//printf("Cant normalize ZERO vector\n");
			//assert(0);// yes this could go here
			d = 1.0f;
		}
		d = 1/d;
		x*= d;
		y*= d;
		z*= d;
	}
 

	Vector3 Normalize()	const;
	Vector3 Safenormalize()	const;
 
	// due to ambiguity and inconsistent standards ther are no overloaded operators for mult such as va*vb.
 
	Vector3 Round(float precision);

	Vector3  Vabs()
	{
		return Vector3(std::abs(x),std::abs(y),std::abs(z));
	}

	Vector3  Vabs()	const
	{
		return Vector3(std::abs(x),std::abs(y),std::abs(z));
	}

	inline bool  IsZero() const
	{
		return ((x*x) <= FLT_EPSILON) && ((y*y) <= FLT_EPSILON) && ((z*z) <= FLT_EPSILON );
	} 
};





inline Vector3 Vector3::Normalize()	const
{
	// this routine, normalize, is ok, provided magnitude works!!
	float d=Magnitude();
	if (d==0)
	{
		printf("Cant normalize ZERO vector\n");
		assert(0);// yes this could go here
		d=0.1f;
	}
	d = 1/d;
	return Vector3(x*d,y*d,z*d);
}

inline Vector3 Vector3::Safenormalize() const
{
	if(Magnitude()<=0.0f)
	{
		return Vector3(0,0,0);
	}
	return Normalize();
}

inline Vector3 Vector3::Round(float precision)
{
	return Vector3(Roundf(x,precision), Roundf(y,precision), Roundf(z,precision));
}


inline Vector3 operator*( float s, const Vector3 &v) 
{
	return Vector3(v.x*s,v.y*s,v.z*s);
}

 Vector3 Vec3Interpolate(const Vector3 &v0,const Vector3 &v1,float alpha);  

 Vector3 Vec3VectorMin(const Vector3 &a,const Vector3 &b);
 Vector3 Vec3VectorMax(const Vector3 &a,const Vector3 &b);

// the statement v1*v2 is ambiguous since there are 3 types
// of vector multiplication
//  - componantwise (for example combining colors)
//  - dot product
//  - cross product
// Therefore we never declare/implement this function.
// So we will never see:  Vector3 operator*(Vector3 a,Vector3 b) 



//-----------------Matrix33---------------
class  Matrix33
{
public:
	Vector3 x,y,z;  // the 3 rows of the Matrix
	Matrix33(){}
	Matrix33(float xx,float xy,float xz,float yx,float yy,float yz,float zx,float zy,float zz):x(xx,xy,xz),y(yx,yy,yz),z(zx,zy,zz){}
	Matrix33(Vector3 _x,Vector3 _y,Vector3 _z):x(_x),y(_y),z(_z){}
	Vector3&       operator[](int i)       {assert(i>=0&&i<3);return (&x)[i];}
	const Vector3& operator[](int i) const {assert(i>=0&&i<3);return (&x)[i];}
	float&        operator()(int r, int c)       {assert(r>=0&&r<3&&c>=0&&c<3);return ((&x)[r])[c];}
	const float&  operator()(int r, int c) const {assert(r>=0&&r<3&&c>=0&&c<3);return ((&x)[r])[c];}
 
	Vector3	GetColumn(int index) const;
	Vector3	GetRow(int index) const;
	__inline void SetIndetity()
	{
		x.x = 1.f;
		x.y = 0.0f;
		x.z = 0.0f; 

		y.x = 0.f;
		y.y = 1.f;
		y.z = 0.0f; 

		z.x = 0.f;
		z.y = 0.0f;
		z.z = 1.f; 
 
	}
	//------------ Matrix33 ---------------
	__inline float Determinant() //行列式
	{
		return  x.x*y.y*z.z + y.x*z.y*x.z + z.x*x.y*y.z
			-x.x*z.y*y.z - y.x*x.y*z.z - z.x*y.y*x.z ;
	}

	__inline Matrix33 Inverse()
	{
		Matrix33 b;
		float d=Determinant();
		assert(d!=0);
		for(int i=0;i<3;i++) 
		{
			for(int j=0;j<3;j++) 
			{
				int i1=(i+1)%3;
				int i2=(i+2)%3;
				int j1=(j+1)%3;
				int j2=(j+2)%3;
				// reverse indexs i&j to take transpose
				b[j][i] = ((&x)[i1][j1]*(&x)[i2][j2]-(&x)[i1][j2]*(&x)[i2][j1])/d;
			}
		}
		// Matrix check=a*b; // Matrix 'check' should be the identity (or close to it)
		return b;
	}


	__inline Matrix33 Transpose()
	{
		return Matrix33( Vector3(x.x,y.x,z.x),
			Vector3(x.y,y.y,z.y),
			Vector3(x.z,y.z,z.z));
	}
 
	__inline Matrix33 operator*(const float& s )  
	{ 
		return Matrix33(x*s, y*s ,z*s); 
	}
		
	__inline Matrix33 operator/(const float& s )  
	{ 
		float t=1/s;
		return Matrix33(x*t, y*t ,z*t); 
	}
		
	__inline Matrix33 operator+(const Matrix33& b )
	{
		return Matrix33(x+b.x, y+b.y, z+b.z);
	}
		
	__inline Matrix33 operator-(const Matrix33& b )
	{
		return Matrix33(x-b.x, y-b.y, z-b.z);
	}
		
	__inline Matrix33 &operator+=(const Matrix33& b )
	{
		x+=b.x;
		y+=b.y;
		z+=b.z;
		return *this;
	}
		
	__inline Matrix33 &operator-=(const Matrix33& b )
	{
		x-=b.x;
		y-=b.y;
		z-=b.z;
		return *this;
	}
		
	__inline Matrix33 &operator*=(const float& s )
	{
		x*=s;
		y*=s;
		z*=s;
		return *this;
	}
 
}; 


inline Vector3 operator*(const Vector3& v, const Matrix33 &m) {
	return Vector3((m.x.x*v.x + m.y.x*v.y + m.z.x*v.z),
		(m.x.y*v.x + m.y.y*v.y + m.z.y*v.z),
		(m.x.z*v.x + m.y.z*v.y + m.z.z*v.z));
}

inline Vector3 operator*(const Matrix33 &m, const Vector3& v) {
	return Vector3(m.x.Dot(v), m.y.Dot(v), m.z.Dot(v));
}


Matrix33 operator*(const Matrix33& a, const Matrix33& b);
//-------- 4D Math --------


class  Vector4
{
public:
	float x,y,z,w;
	Vector4(){x=0;y=0;z=0;w=0;};
	Vector4(float _x,float _y,float _z,float _w){x=_x;y=_y;z=_z;w=_w;}
	Vector4( const float *_in ){x = _in[0];y = _in[1]; z= _in[2];w = _in[3];}
	Vector4(const Vector3 &v,float _w){x=v.x;y=v.y;z=v.z;w=_w;}
	//operator float *() { return &x;};
	float& operator[](int i) {assert(i>=0&&i<4);return ((float*)this)[i];}
	const float& operator[](int i) const {assert(i>=0&&i<4);return ((float*)this)[i];}
	const Vector3& xyz() const { return *((Vector3*)this);}
	Vector3&       xyz()       { return *((Vector3*)this);}
 
	int operator==(const Vector4 &b ) 
	{
		return (x==b.x && y==b.y && z==b.z && w==b.w); 
	}


	operator       class Vector4* ()       { return (class Vector4*) this;}
	operator const class Vector4* () const { return (class Vector4*) this;}
	//  Dont implement m*v for now, since that might confuse us
	//  All our transforms are based on multiplying the "row" vector on the left
 
	__inline Vector4 operator*( float s)	const
	{
		return Vector4(x*s,y*s,z*s,w*s);
	}

	__inline Vector4  operator/(float s )	const
	{ 
		return Vector4(x*(1.0f/s),y*(1.0f/s),z*(1.0f/s),w*(1.0f/s)); 
	}
 

	__inline float   Dot(const Vector4& b )    const
	{
		return  x*b.x +  y*b.y + z*b.z + w*b.w; 
	}
 

	Vector4 CMul(const Vector4 &b)	const
	{
		return Vector4(x*b.x,y*b.y,z*b.z,w*b.w);
	}


 
	Vector4 operator+(const Vector4 &b)	const
	{
		return Vector4(x+b.x,y+b.y,z+b.z,w+b.w);
	}



	Vector4 operator-(const Vector4 &b)		const
	{
		return Vector4(x-b.x,y-b.y,z-b.z,w-b.w);
	}

	__inline Vector4&  operator/=(float s )
	{
		float sinv = 1.0f / s;
		x *= sinv;
		y *= sinv;
		z *= sinv;
		w *= sinv;
		return *this;
	}

	float	Magnitude()	const
	{
		return std::sqrt(Sqr(w)+Sqr(x)+Sqr(y)+Sqr(z));
	}

};

 

class  Matrix44
{
public:
	Vector4 x,y,z,w;  // the 4 rows
	Matrix44(){}
	Matrix44(const Vector4 &_x, const Vector4 &_y, const Vector4 &_z, const Vector4 &_w):x(_x),y(_y),z(_z),w(_w){}
	Matrix44(float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23, 
		float m30, float m31, float m32, float m33 )
		:x(m00,m01,m02,m03),y(m10,m11,m12,m13),z(m20,m21,m22,m23),w(m30,m31,m32,m33){}
	float&       operator()(int r, int c)       {assert(r>=0&&r<4&&c>=0&&c<4);return ((&x)[r])[c];}
	const float& operator()(int r, int c) const {assert(r>=0&&r<4&&c>=0&&c<4);return ((&x)[r])[c];}
	operator       float* ()       {return &x.x;}
	operator const float* () const {return &x.x;}

#ifdef _USING_DXX
	operator       struct D3DXMATRIX* ()       { return (struct D3DXMATRIX*) this;}
	operator const struct D3DXMATRIX* () const { return (struct D3DXMATRIX*) this;}

	operator		D3DMATRIX* ()       { return (D3DMATRIX*) this;}
	operator const  D3DMATRIX* () const { return (D3DMATRIX*) this;}
#endif
	Vector4	GetColumn(int index)  const;
	Vector4	GetRow(int index)  const;

	Matrix33		GetRotationPart()	const;
	void	SetRotationPart(Matrix33& mat33);

	Vector3		GetPositionPart();
	const Vector3		GetPositionPart()	const;

	void	SetPositionPart(Vector3& vecPos);

 
	void GetOpenGLMatrix(float *m) const 
	{
		m[0]  = float(x.x); 
		m[1]  = float(y.x);
		m[2]  = float(z.x);
		m[3]  = float(w.x); 
		m[4]  = float(x.y);
		m[5]  = float(y.y);
		m[6]  = float(z.y);
		m[7]  = float(w.y); 
		m[8]  = float(x.z); 
		m[9]  = float(y.z);
		m[10] = float(z.z);
		m[11] = float(w.z);
		m[12] = float(x.w);
		m[13] = float(y.w);
		m[14] = float(z.w);
		m[15] = float(w.w);
	}

	__inline void SetIndetity()
	{
		x.x = 1.f;
		x.y = 0.0f;
		x.z = 0.0f;
		x.w = 0.0f;

		y.x = 0.f;
		y.y = 1.f;
		y.z = 0.0f;
		y.w = 0.0f;

		z.x = 0.f;
		z.y = 0.0f;
		z.z = 1.f;
		z.w = 0.0f;

		w.x = 0.0f;
		w.y = 0.0f;
		w.z = 0.0f;
		w.w = 1.f;
	}
 
		
	
	__inline Matrix44 Transpose()	const
	{
		return Matrix44(
			x.x, y.x, z.x, w.x,
			x.y, y.y, z.y, w.y,
			x.z, y.z, z.z, w.z,
			x.w, y.w, z.w, w.w );
	}
 
	__inline Matrix44 Inverse()	const
	{
		Matrix44 d;
		float *dst = &d.x.x;
		float tmp[12]; /* temp array for pairs */
		float src[16]; /* array of transpose source matrix */
		float det; /* determinant */
		/* transpose matrix */
		for ( int i = 0; i < 4; i++) {
			src[i] = (*this)(i,0) ;
			src[i + 4] = (*this)(i,1);
			src[i + 8] = (*this)(i,2);
			src[i + 12] = (*this)(i,3); 
		}
		/* calculate pairs for first 8 elements (cofactors) */
		tmp[0]  = src[10] * src[15];
		tmp[1]  = src[11] * src[14];
		tmp[2]  = src[9] * src[15];
		tmp[3]  = src[11] * src[13];
		tmp[4]  = src[9] * src[14];
		tmp[5]  = src[10] * src[13];
		tmp[6]  = src[8] * src[15];
		tmp[7]  = src[11] * src[12];
		tmp[8]  = src[8] * src[14];
		tmp[9]  = src[10] * src[12];
		tmp[10] = src[8] * src[13];
		tmp[11] = src[9] * src[12];
		/* calculate first 8 elements (cofactors) */
		dst[0]  = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
		dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
		dst[1]  = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
		dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
		dst[2]  = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
		dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
		dst[3]  = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
		dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
		dst[4]  = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
		dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
		dst[5]  = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
		dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
		dst[6]  = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
		dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
		dst[7]  = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
		dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
		/* calculate pairs for second 8 elements (cofactors) */
		tmp[0]  = src[2]*src[7];
		tmp[1]  = src[3]*src[6];
		tmp[2]  = src[1]*src[7];
		tmp[3]  = src[3]*src[5];
		tmp[4]  = src[1]*src[6];
		tmp[5]  = src[2]*src[5];
		tmp[6]  = src[0]*src[7];
		tmp[7]  = src[3]*src[4];
		tmp[8]  = src[0]*src[6];
		tmp[9]  = src[2]*src[4];
		tmp[10] = src[0]*src[5];
		tmp[11] = src[1]*src[4];
		/* calculate second 8 elements (cofactors) */
		dst[8]  = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
		dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
		dst[9]  = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
		dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
		dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
		dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
		dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
		dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
		dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
		dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
		dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
		dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
		dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
		dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
		dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
		dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
		/* calculate determinant */
		det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];
		/* calculate matrix inverse */
		det = 1/det;
		for ( int j = 0; j < 16; j++)
			dst[j] *= det;
		return d;
	}
 

};


inline Vector4 operator*( float s, const Vector4 &v) 
{
	return Vector4(v.x*s,v.y*s,v.z*s,v.w*s);
}


//--------------- 4D ----------------

 Vector4   operator*( const Vector4& v, const Matrix44& m );

//Vector4   operator*( const Vector4& v, const Matrix44* m );

 Matrix44 operator*( const Matrix44& a, const Matrix44& b );


 Vector3*	Vec3TransformCoord(Vector3* vecOut,const Vector3* vecIn,const Matrix44* mat);
 
 Vector3*	Vec3TransformNormal(Vector3* vecOut,const Vector3* vecIn,const Matrix44* mat);

 
 Vector2*	Vec2TransformCoord(Vector2* vecOut,const Vector2* vecIn,const Matrix33* mat);

 Vector2*	Vec2TransformNormal(Vector2* vecOut,const Vector2* vecIn,const Matrix33* mat);

 inline  Vector4 Vec4Lerp(Vector4& a, Vector4& b, float alpha)
 {	//a + (b - a) * t.
	 return a * (1 - alpha) + b * alpha;
 }
 inline  Vector3 Vec3Lerp(Vector3& a, Vector3& b, float alpha)
 {
	 return a * (1 - alpha) + b * alpha;
 }

class CAABBBox;
 CAABBBox*	AABBBoxTransform(CAABBBox* vecOut,const CAABBBox* vecIn,const Matrix44* mat);

//-------- Quaternion ------------

class  Quaternion :public Vector4
{
public:
	Quaternion() { x = y = z = 0.0f; w = 1.0f; }
	Quaternion( Vector3 v, float t ) 
	{ 
		if(v.MagnitudeSq() == 0.0f) //invalid nan
		{
				x = 1;
				y = 0;
				z = 0;
				w = 0;
		}
		else
		{
			v = v.Normalize();
			w = std::cos(t/2.0f); 
			v = v*std::sin(t/2.0f); 
			x = v.x; y = v.y; z = v.z; 
		}

	}
	Quaternion(float _x, float _y, float _z, float _w){x=_x;y=_y;z=_z;w=_w;}
	float Angle() const { return std::acos(w)*2.0f; }
	Vector3 Axis() const { Vector3 a(x,y,z); if(std::abs(Angle())<0.0000001f) return Vector3(1,0,0); return a*(1/std::sin(Angle()/2.0f)); }
	Vector3 XDir() const { return Vector3( 1-2*(y*y+z*z),  2*(x*y+w*z),  2*(x*z-w*y) ); }
	Vector3 YDir() const { return Vector3(   2*(x*y-w*z),1-2*(x*x+z*z),  2*(y*z+w*x) ); }
	Vector3 ZDir() const { return Vector3(   2*(x*z+w*y),  2*(y*z-w*x),1-2*(x*x+y*y) ); }
	Matrix33 Getmatrix() const { return Matrix33( XDir(), YDir(), ZDir() ); }
	operator Matrix33() { return Getmatrix(); }
 
	const Vector3& xyz() const { return *((Vector3*)this);}
	Vector3&       xyz()       { return *((Vector3*)this);}

	Quaternion operator*(const Quaternion& b )	const
	{
		Quaternion c;
		//c.w = w*b.w - x*b.x - y*b.y - z*b.z;
		//c.x = w*b.x + x*b.w + y*b.z - z*b.y; 
		//c.y = w*b.y - x*b.z + y*b.w + z*b.x; 
		//c.z = w*b.z + x*b.y - y*b.x + z*b.w; 
		c.w = b.w*w - b.x*x - b.y*y - b.z*z;
		c.x = b.w*x + b.x*w + b.y*z - b.z*y; 
		c.y = b.w*y - b.x*z + b.y*w + b.z*x; 
		c.z = b.w*z + b.x*y - b.y*x + b.z*w; 
		return c;
	}

	Quaternion operator*(float b )	const
	{
		return Quaternion(x*b, y*b, z*b ,w*b);
	}


	Quaternion  Inverse()	const
	{
		return Quaternion(-x,-y,-z,w);
	}

	Quaternion& operator*=(const float s )
	{
		x *= s;
		y *= s;
		z *= s;
		w *= s;
		return *this;
	}

	void	SetNormalize()
	{
		float m = std::sqrt(Sqr(w)+Sqr(x)+Sqr(y)+Sqr(z));
		if(m<0.000000001f) {
			x = 0;
			y = 0;
			z = 0;
			w = 1; 
		}
		else
		{
			x*=(1.0f/m);
			y*=(1.0f/m);
			z*=(1.0f/m);
			w*=(1.0f/m);
		}

	}

	float	Magnitude()	const
	{
		return std::sqrt(Sqr(w)+Sqr(x)+Sqr(y)+Sqr(z));
	}

	float	MagnitudeSq()	const
	{
		return	Sqr(w)+Sqr(x)+Sqr(y)+Sqr(z);
	}

	Quaternion Normalize()	const
	{
		float m = std::sqrt(Sqr(w)+Sqr(x)+Sqr(y)+Sqr(z));
		if(m<0.000000001f) {
			return Quaternion(0,0,0,1);
		}
		return Quaternion(x*(1.0f/m),y*(1.0f/m),z*(1.0f/m),w*(1.0f/m));
	}

	Vector3 TransformVector3Coord(const Vector3& v )	const
	{
		// The following is equivalent to:
		//return (q.getmatrix() * v);
		float qx2 = x*x;
		float qy2 = y*y;
		float qz2 = z*z;

		float qxqy = x*y;
		float qxqz = x*z;
		float qxqw = x*w;
		float qyqz = y*z;
		float qyqw = y*w;
		float qzqw = z*w;
		return Vector3(
			(1-2*(qy2+qz2))*v.x + (2*(qxqy-qzqw))*v.y + (2*(qxqz+qyqw))*v.z ,
			(2*(qxqy+qzqw))*v.x + (1-2*(qx2+qz2))*v.y + (2*(qyqz-qxqw))*v.z ,
			(2*(qxqz-qyqw))*v.x + (2*(qyqz+qxqw))*v.y + (1-2*(qx2+qy2))*v.z  );
	}

	//Vector3 operator*(const Vector3& v )	const
	//{
	//	// The following is equivalent to:
	//	//return (q.getmatrix() * v);
	//	float qx2 = x*x;
	//	float qy2 = y*y;
	//	float qz2 = z*z;

	//	float qxqy = x*y;
	//	float qxqz = x*z;
	//	float qxqw = x*w;
	//	float qyqz = y*z;
	//	float qyqw = y*w;
	//	float qzqw = z*w;
	//	return Vector3(
	//		(1-2*(qy2+qz2))*v.x + (2*(qxqy-qzqw))*v.y + (2*(qxqz+qyqw))*v.z ,
	//		(2*(qxqy+qzqw))*v.x + (1-2*(qx2+qz2))*v.y + (2*(qyqz-qxqw))*v.z ,
	//		(2*(qxqz-qyqw))*v.x + (2*(qyqz+qxqw))*v.y + (1-2*(qx2+qy2))*v.z  );
	//}


	Quaternion operator+(const Quaternion& b )	const
	{
		return Quaternion(x+b.x, y+b.y, z+b.z, w+b.w);
	}

	float Dot(const Quaternion &b )
	{
		return  (w*b.w + x*b.x + y*b.y + z*b.z);
	}


	float Dot(const Quaternion &b )	const
	{
		return  (w*b.w + x*b.x + y*b.y + z*b.z);
	}
 

	void SetYawPitchRollInRadian( float yaw, float pitch, float roll ) 
	{
		//*this = Quaternion(Vector3(0.0f,0.0f,1.0f),yaw)*Quaternion(Vector3(1.0f,0.0f,0.0f),pitch)*Quaternion(Vector3(0.0f,1.0f,0.0f),roll);
		float halfYaw = float(yaw) * float(0.5);  
		float halfPitch = float(pitch) * float(0.5);  
		float halfRoll = float(roll) * float(0.5);  
		float cosYaw = std::cos(halfYaw);
		float sinYaw = std::sin(halfYaw);
		float cosPitch = std::cos(halfPitch);
		float sinPitch = std::sin(halfPitch);
		float cosRoll = std::cos(halfRoll);
		float sinRoll = std::sin(halfRoll);
		x = cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw;
		y = cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw;
		z = sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw;
		w = cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw;
	}

	Quaternion SetYawPitchRollInDegree( float yaw, float pitch, float roll ) 
	{
		roll  = FLT_DTOR(roll);
		yaw   = FLT_DTOR(yaw);
		pitch = FLT_DTOR(pitch);
		return Quaternion(Vector3(0.0f,0.0f,1.0f),yaw)*Quaternion(Vector3(1.0f,0.0f,0.0f),pitch)*Quaternion(Vector3(0.0f,1.0f,0.0f),roll);
	}

	float Yaw()
	{
		Vector3 v;
		v= YDir();
		return (v.y==0.0&&v.x==0.0) ? 0.0f: FLT_RTOD(std::atan2(-v.x,v.y));
	}

	float Pitch()
	{
		Vector3 v;
		v= YDir();
		return std::atan2(v.z,FLT_RTOD(std::sqrt(Sqr(v.x)+Sqr(v.y))));
	}

	float Roll()
	{
		Quaternion q = Quaternion(Vector3(0.0f,0.0f,1.0f),-FLT_RTOD(Yaw()))*  *this;
		q = Quaternion(Vector3(1.0f,0.0f,0.0f),-FLT_RTOD(Pitch()))  *q;
		return FLT_RTOD(std::atan2(-q.XDir().z,q.XDir().x));
	}

	float Yaw( const Vector3& v )
	{
		return (v.y==0.0&&v.x==0.0) ? 0.0f: FLT_RTOD(std::atan2(-v.x,v.y));
	}

	float Pitch( const Vector3& v )
	{
		return FLT_RTOD(std::atan2(v.z,std::sqrt(Sqr(v.x)+Sqr(v.y))));
	} 

};

inline Quaternion
	operator*(const Quaternion& q, const Vector3& w)
{
	return Quaternion( 
		q.w * w.x + q.y * w.z - q.z * w.y,
		q.w * w.y + q.z * w.x - q.x * w.z,
		q.w * w.z + q.x * w.y - q.y * w.x,
		-q.x * w.x - q.y * w.y - q.z * w.z); 
}

inline Quaternion
	operator*(const Vector3& w, const Quaternion& q)
{
	return Quaternion( 
		+w.x * q.w + w.y * q.z - w.z * q.y,
		+w.y * q.w + w.z * q.x - w.x * q.z,
		+w.z * q.w + w.x * q.y - w.y * q.x,
		-w.x * q.x - w.y * q.y - w.z * q.z); 
}
//------- Plane ----------

class  Plane
{
public:
	Vector3	normal;
	float	dist;   // distance below origin - the D from plane equasion Ax+By+Cz+D=0
	Plane(const Vector3 &n,float d):normal(n),dist(d){}
	Plane(const Vector3 &ptInPlane,const Vector3 &planeNormal);
	Plane(const Vector4 v):normal(v.xyz()),dist(v.w){}
	Plane(float a,float b,float c,float d):normal(a,b,c),dist(d){}
	Plane():normal(),dist(0){}
	Plane(const Vector3 &pt0,const Vector3& pt1,Vector3& pt2);
	void	Transform(const Vector3 &position, const Quaternion &orientation);
	void	Transform(const Matrix44& transMat);  	
		
	inline Plane PlaneFlip(const Plane &plane){return Plane(-plane.normal,-plane.dist);}
	inline int operator==(const Plane &b ) { return (normal==b.normal && dist==b.dist); }
	inline int coplanar(const Plane &b ) { return (*this==b || *this==PlaneFlip(b)); }

	operator       class Vector4* ()       { return (class Vector4*) this;}

	operator       float* ()       {return &normal.x;} 


	Plane operator*(float b )
	{
		return Plane(normal*b ,dist*b);
	}

	Plane operator*(float b )	const
	{
		return Plane(normal*b ,dist*b);
	}


	float Dot(const Vector3& v)
	{
		return normal.Dot(v) + dist;
	}

	float Dot(const Vector3& v)	const
	{
		return normal.Dot(v) + dist;
	}

	float Dot(const Vector4& v)
	{
		return normal.Dot(v.xyz()) + dist*v.w;
	}

	float Dot(const Vector4& v)	const
	{
		return normal.Dot(v.xyz()) + dist*v.w;
	}

	void	SetNormalize()
	{
		normal.SetNormalize();
	}
	Plane Normalize()
	{
		float nSqrt = normal.Magnitude();
		if(nSqrt)
		{
			return Plane(normal.x/nSqrt,normal.y/nSqrt,normal.z/nSqrt,dist/nSqrt);
		}
		else
			return Plane(0.f,0.f,0.f,0.f); 
	}

	Plane Normalize()	const
	{
		return Plane(normal.Normalize(),dist);
	}

	float SolveForX(float Y, float Z)const;
	float SolveForY(float X, float Z)const;
	float SolveForZ(float X, float Y)const;
};

//:	SolveForX
//----------------------------------------------------------------------------------------
//
//	Given Z and Y, Solve for X on the plane 
//
//-------------------------------------------------------------------------------------://
inline float Plane::SolveForX(float Y, float Z)const
{
	//Ax + By + Cz + D = 0
	// Ax = -(By + Cz + D)
	// x = -(By + Cz + D)/A

	if (normal[0])
	{
		return ( -(normal[1]*Y + normal[2]*Z + dist) / normal[0] );
	}

	return (0.0f);
}

//:	SolveForY
//----------------------------------------------------------------------------------------
//
//	Given X and Z, Solve for Y on the plane 
//
//-------------------------------------------------------------------------------------://
inline float Plane::SolveForY(float X, float Z)const
{
	//Ax + By + Cz + D = 0
	// By = -(Ax + Cz + D)
	// y = -(Ax + Cz + D)/B

	if (normal[1])
	{
		return ( -(normal[0]*X + normal[2]*Z + dist) / normal[1] );
	}

	return (0.0f);
}

//:	SolveForZ
//----------------------------------------------------------------------------------------
//
//	Given X and Y, Solve for Z on the plane 
//
//-------------------------------------------------------------------------------------://
inline float Plane::SolveForZ(float X, float Y)const
{
	//Ax + By + Cz + D = 0
	// Cz = -(Ax + By + D)
	// z = -(Ax + By + D)/C

	if (normal[2])
	{
		return ( -(normal[0]*X + normal[1]*Y + dist) / normal[2] );
	}

	return (0.0f);
}

 
inline unsigned int F2DW( float f ) { return *((unsigned int*)&f); }

//---------------------RectI-------------------
struct  RectI
{
	long    left;
	long    top;
	long    right;
	long    bottom;

	RectI()
	{
		left = top = right = bottom = 0;
	}

	RectI(int inLeft,int inTop,int inRight,int inBottom)
	{
		left	= inLeft;
		top		= inTop;
		right	= inRight;
		bottom	= inBottom;
		
	}

	RectI(VectorI2 pt,VectorI2 extent)
	{
		left = pt.x;
		top  = pt.y;
		right = pt.x + extent.x;
		bottom = pt.y + extent.y;
		
	}

	void	SetAll(long value)
	{
		left = top = right = bottom = value;
	}

	void Set(long in_left, long in_right, long in_top, long in_bottom)
	{
		top = in_top;
		bottom = in_bottom;
		left = in_left;
		right = in_right;
	}

	int		GetWidth() const
	{
		return right - left;
	}

	int		GetHeight() const
	{
		return bottom - top;
	}

	VectorI2	GetExtent()
	{
		return VectorI2(right - left,bottom - top);
	}

	bool	IsValidRect() const 
	{
		return (right - left> 0 && bottom - top > 0); 
	}

	void CollapseRect( RectI &rectRef )
	{
		// Inset by padding
		left	+= rectRef.left;
		top		+= rectRef.top;
		right	-= rectRef.right;
		bottom	-= rectRef.bottom;
	}
	void ExpandRect( RectI &rectRef )
	{
		// Inset by padding
		left	-= rectRef.left;
		top		-= rectRef.top	;
		right	+= rectRef.right;
		bottom	+= rectRef.bottom	;
	}

	void	InflateRect(int x,int y)
	{
		left -= x;
		right += x;
		top -= y;
		bottom += y;
	}

	void	OffsetRect(int x,int y)
	{
		left += x;
		right += x;
		top += y;
		bottom += y;
	}
#ifdef WIN32
	operator RECT&()
	{
		return *this;
	}
#endif

#ifdef ANDROID
	operator ARect&()
	{
		return *this;
	}
#endif


	bool	operator ==(RectI& rc)
	{
		return left == rc.left && right == rc.right && top == rc.top && bottom == rc.bottom; 
	}

	bool	operator !=(RectI& rc)
	{
		return left != rc.left || right != rc.right || top != rc.top || bottom != rc.bottom; 
	}

	bool	operator ==(const RectI& rc)
	{
		return left == rc.left && right == rc.right && top == rc.top && bottom == rc.bottom; 
	}

	bool	operator !=(const RectI& rc)
	{
		return left != rc.left || right != rc.right || top != rc.top || bottom != rc.bottom; 
	}

	bool	operator ==(const RectI& rc)	const
	{
		return left == rc.left && right == rc.right && top == rc.top && bottom == rc.bottom; 
	}

	bool	operator !=(const RectI& rc)  const
	{
		return left != rc.left || right != rc.right || top != rc.top || bottom != rc.bottom; 
	}


	bool	Intersect(const RectI& rc)	const
	{
		int tmpleft =  std::max(left,rc.left);
		int tmpright = std::min(right,rc.right);
		int tmptop = std::max(top,rc.top);
		int tmpbottom = std::min(bottom,rc.bottom);

		if(tmpright < tmpleft || tmptop > tmpbottom)
			return false;
		return true;		
	}

	bool	PointInRect(int x,int y)
	{
		return left <= x && x<= right && y>= top && y<= bottom;
	}
	 
};


//--------------------------RectF---------------------------


struct  RectF
{
	float left;
	float top;
	float right;
	float bottom;

	RectF()
	{
		left = top = right = bottom = 0.0f;
	}

	RectF(float inLeft,float inTop,float inRight,float inBottom)
	{
		left	= inLeft;
		top		= inTop;
		right	= inRight;
		bottom	= inBottom;
		
	}

	RectF(Vector2 pt,Vector2 extent)
	{
		left = pt.x;
		top  = pt.y;
		right = pt.x + extent.x;
		bottom = pt.y + extent.y;
		
	}

	void	SetAll(float value)
	{
		left = top = right = bottom = value;
	}

	void Set(  float in_left, float in_right, float in_top, float in_bottom)
	{
		top = in_top;
		bottom = in_bottom;
		left = in_left;
		right = in_right;
	}

	float		GetWidth() const
	{
		return right - left;
	}

	float		GetHeight() const
	{
		return bottom - top;
	}

	Vector2	GetExtent()
	{
		return Vector2(right - left,bottom - top);
	}

	bool	IsValidRect() const 
	{
		return (right - left> FLT_EPSILON && bottom - top > FLT_EPSILON); 
	}

	void CollapseRect( RectF &rectRef )
	{
		// Inset by padding
		left	+= rectRef.left;
		top		+= rectRef.top;
		right	-= rectRef.right;
		bottom	-= rectRef.bottom;
	}
	void ExpandRect( RectF &rectRef )
	{
		// Inset by padding
		left	-= rectRef.left;
		top		-= rectRef.top	;
		right	+= rectRef.right;
		bottom	+= rectRef.bottom	;
	}

	void	InflateRect(float x,float y)
	{
		left += x;
		right -= x;
		top += y;
		bottom -= y;
	}

	void	OffsetRect(float x,float y)
	{
		left += x;
		right += x;
		top += y;
		bottom += y;
	}

	bool	operator ==(RectF& rc)
	{
		return (left - rc.left)<FLT_EPSILON && (right - rc.right)<FLT_EPSILON  && (top - rc.top)<FLT_EPSILON  && (bottom - rc.bottom)<FLT_EPSILON ; 
	}

	bool	operator !=(RectF& rc)
	{
		return (left - rc.left)>FLT_EPSILON ||(right - rc.right)>FLT_EPSILON  || (top - rc.top)>FLT_EPSILON  || (bottom - rc.bottom)>FLT_EPSILON ; 
	}

	bool	operator ==(const RectI& rc)
	{
		return (left - rc.left)<FLT_EPSILON && (right - rc.right)<FLT_EPSILON  && (top - rc.top)<FLT_EPSILON  && (bottom - rc.bottom)<FLT_EPSILON ; 
	}

	bool	operator !=(const RectI& rc)
	{
		return (left - rc.left)>FLT_EPSILON ||(right - rc.right)>FLT_EPSILON  || (top - rc.top)>FLT_EPSILON  || (bottom - rc.bottom)>FLT_EPSILON ; 
	}

	bool	operator ==(const RectI& rc)	const
	{
		return (left - rc.left)<FLT_EPSILON && (right - rc.right)<FLT_EPSILON  && (top - rc.top)<FLT_EPSILON  && (bottom - rc.bottom)<FLT_EPSILON ; 
	}

	bool	operator !=(const RectI& rc)  const
	{
		return (left - rc.left)>FLT_EPSILON ||(right - rc.right)>FLT_EPSILON  || (top - rc.top)>FLT_EPSILON  || (bottom - rc.bottom)>FLT_EPSILON ; 
	}


	bool	Intersect(const RectF& rc)	const
	{
		float tmpleft = std::max (left,rc.left);
		float tmpright = std::min(right,rc.right);
		float tmptop = std::max(top,rc.top);
		float tmpbottom = std::min(bottom,rc.bottom);

		if(tmpright < tmpleft || tmptop > tmpbottom)
			return false;
		return true;		
	}

	bool	PointInRect(float x,float y)
	{
		return left <= x && x<= right && y>= top && y<= bottom;
	}

	bool	PointInRect(int x,int y)
	{
		return PointInRect((float)x,(float)y);
	}
	 
};


//--------- Utility Functions ------
//const var
static const Matrix44 Cst_MatIndetity(1.0f,0.0f,0.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,
	0.0f,0.0f,1.0f,0.0f,
	0.0f,0.0f,0.0f,1.0f);
 
 
static Vector3 Cst_X_AXIS(1.f,0.f,0.f);
static Vector3 Cst_Y_AXIS(0.f,1.f,0.f);
static Vector3 Cst_Z_AXIS(0.f,0.f,1.f);

static Vector2 Cst_ZERO_VEC2(0.f,0.f);
static Vector3 Cst_ZERO_VEC3(0.f,0.f,0.f);
static Vector3 Cst_ONE_VEC3(1.f,1.f,1.f);
static Vector3 Cst_HALFONE_VEC3(0.5,0.5,0.5);
//function
// Quaternion	QuaternionInterpolate( Quaternion a, const Quaternion& b, float interp );
 Quaternion	QuaternionSlerp(const Quaternion &q0,const Quaternion &q1,float alpha) ;

 Quaternion	QuaternionInterpolate(const Quaternion &q0,const Quaternion &q1,float alpha) ;

 Quaternion	QuaternionExtrapolate(const Quaternion &q0,const Quaternion &q1,float t) ;

 bool			QuaternionIsNan(const Quaternion &q0);

 Quaternion	RotationArc(Vector3 v0, Vector3 v1 ); // returns quat q where q*v0=v1 

 Matrix44		MatrixFromQuatVec(const Quaternion &q, const Vector3 &v);

 Matrix44		MatrixRigidInverse(const Matrix44 &m);

 Matrix44*		MatrixPerspectiveLH(   Matrix44 *pOut,
	     float w,
	 float h,
	 float zn,
	 float zf );
 Matrix44*		MatrixPerspectiveFovLH(  Matrix44 *pOut,
	     float fovy, 
	     float aspect, 
	     float zn, 
	     float zf );
 Matrix44*		MatrixPerspectiveFovRH(  Matrix44 *pOut,
	     float fovy, 
	     float aspect, 
	     float zn, 
	     float zf );

 Matrix44*		MatrixLookAtLH(  Matrix44 *pOut,
	     const Vector3* eye, 
	     const Vector3* at, 
	     const Vector3* up);

 Matrix44*		MatrixLookAtRH(  Matrix44 *pOut,
	     const Vector3* eye, 
	     const Vector3* at, 
	     const Vector3* up);

 Matrix44*		MatrixOrthoLH(
	  Matrix44 *pOut,
	 float w,
	 float h,
	 float zn,
	 float zf
	);

 Matrix44*		MatrixOrthoRH(
	  Matrix44 *pOut,
	 float w,
	 float h,
	 float zn,
	 float zf
	);
 Matrix44*		MatrixLookViewDir(
	  Matrix44 *pOut,
	     Vector3* pVecEye,
	     Vector3* pVecDir,
			Vector3*	pVecUp);
	

 Matrix33*		Matrix33RotateAngle( Matrix33* pOutMat,const float angle_radians);

 Matrix44*		MatrixTranslation( Matrix44* pOutMat,     const Vector3 &t);
 Matrix44*		MatrixTranslation( Matrix44* pOutMat,     float x,     float y,     float z);
 Matrix44*		MatrixScaling( Matrix44* pOutMat,      float x,     float y,     float z);

 Matrix44*		MatrixRotationX( Matrix44* pOutMat,     const float angle_radians);
 Matrix44*		MatrixRotationY( Matrix44* pOutMat,     const float angle_radians);
 Matrix44*		MatrixRotationZ( Matrix44* pOutMat,     const float angle_radians);

 Matrix44*		MatrixMultiply( Matrix44* pOutMat,     const Matrix44* pInMat0,     const Matrix44* pInMat1);

 Matrix44*		MatrixRotationYawPitchRoll(  Matrix44 *pOut,
	     float Yaw,
	     float Pitch,
	     float Roll);

 Matrix44*		MatrixRotationQuaternion(  Matrix44 *pOut,
	     const Quaternion* pRot);
// Build a quaternion from a rotation matrix.
 Quaternion* QuaternionRotationMatrix
	( Quaternion *pOut, const Matrix44 *pM);
 Quaternion* QuaternionRotationMatrix
	( Quaternion *pOut, const Matrix33 *pM);

 Matrix44 * MatrixRotationAxis(
	  Matrix44 *pOut,
	     const Vector3 *pV,
	     float Angle
	);

  void MatrixDecompose(
	  Vector3 *pOutScale,
	  Quaternion *pOutRotation,
	  Vector3 *pOutTranslation,
	  const Matrix44 *pM
	);

	Matrix44 * MatrixReflect(
	  Matrix44 *pOut,
	     const Plane *pPlane
	);

	Plane * PlaneTransform(
	  Plane *pOut,
	     const Plane *pP,
	     const Matrix44 *pM
	); //可能有问题
 
 Vector3		PlaneLineIntersection(const Plane &plane, const Vector3 &p0, const Vector3 &p1);
 Vector3		PlaneProject(const Plane &plane, const Vector3 &point);
 

 Vector3		LineProject(const Vector3 &p0, const Vector3 &p1, const Vector3 &a);  // projects a onto infinite line p0p1
 float			LineProjectTime(const Vector3 &p0, const Vector3 &p1, const Vector3 &a);
 Vector3		ThreePlaneIntersection(const Plane &p0,const Plane &p1, const Plane &p2);
 int			PolyHit(const Vector3 *vert,const int n,const Vector3 &v0, const Vector3 &v1, Vector3 *impact=NULL, Vector3 *normal=NULL);
 int			BoxInside(const Vector3 &p,const Vector3 &bmin, const Vector3 &bmax) ;
 int			BoxIntersect(const Vector3 &v0, const Vector3 &v1, const Vector3 &bmin, const Vector3 &bmax, Vector3 *impact);
 float			DistanceBetweenLines(const Vector3 &ustart, const Vector3 &udir, const Vector3 &vstart, const Vector3 &vdir, Vector3 *upoint=NULL, Vector3 *vpoint=NULL);
 Vector3		TriNormal(const Vector3 &v0, const Vector3 &v1, const Vector3 &v2);
 Vector3		NormalOf(const Vector3 *vert, const int n);
 Quaternion	VirtualTrackBall(const Vector3 &cop, const Vector3 &cor, const Vector3 &dir0, const Vector3 &dir1);
 
 float		GetSignedAngleBetweenVectors(Vector3* vec0, Vector3* vec1,const Vector3* vecRightOfVecs);
 
 Vector3  GetClosestPointToRefPt(Vector3* ptVec,unsigned int ptNum,const Vector3& refPt);  
 
 Matrix44	GetTransformMatrixOfObj(Vector3 ObjXAixe,
	Vector3 ObjYAixe,
	Vector3 ObjZAixe);
  
 inline bool PtIsOnPlane(Vector3 p,Plane& plane)
{
	float dist = plane.Dot(p);
	return dist < 0.1 && dist > -0.1;
}
 
 inline Vector3	Vec3DXToOGL(float x,float y,float z)  
{
	return Vector3(x,y,-z);
}
 inline Vector3	Vec3OGLToDX(float x,float y,float z)
{
	return Vector3(x,y,-z);
}
 inline Quaternion QuaternionDXToOGL(float x,float y,float z,float w)
{
	return Quaternion(x,y,-z,w);
}
 inline Quaternion QuaternionOGLToDX(float x,float y,float z,float w)
{
	return Quaternion(x,y,-z,w);
}


 Vector3*	Vec3GetEulerFromQuaternion(Vector3* out,const Quaternion* q); //y pitch,x yaw,z roll
 
 Quaternion*	QuaternionDecompress(Quaternion* out,short x,short y,short z);

 Vector3*	CalcPtByUVInFace(Vector3* pOut,Vector3* p0,Vector3* p1,Vector3* p2,float U,float V);

// function D3DXLineIntersectTri find line and triangle intersect state
// pa face-vertex0 ,pb face-vertex1, pc face-vertex2
// p1 line-start, p2 line-end
// p intersect point to return
 bool  LineIntersectTri(const Vector3 *pa,const Vector3 *pb,
	const Vector3 *pc,const Vector3 *p1,
	const Vector3 *p2,Vector3 *p,float *pU,float* pV);
 bool  RayIntersectTri(const Vector3 *pa,const Vector3 *pb,
	const Vector3 *pc,const Vector3 *RayOrigin,
	const Vector3 *Raydir,float *pU,float* pV,float* pDist );

 float PlaneIntersectRay( Vector3* ptOut,
	 Plane* plane, 
	 const Vector3* rayDir, 
	 const Vector3*  rayOrigin );

 bool	IsPointInsideTri(const Vector3 *pa,const Vector3 *pb,
	const Vector3 *pc,const Vector3 *p,float *pU,float* pV);

/// Simple reflection equation - pass in a vector and a normal to reflect off of
 inline Vector3* Vec3ReflectNormal( Vector3* outVec, Vector3* inVec, const Vector3* normal )
{
	*outVec = *inVec - (*normal) * ( inVec->Dot(*normal )) * 2.0f ;
	return outVec;
}

 inline Vector3* Vec3ReflectNormal( Vector3* outVec, Vector3* inVec, const Plane* pPlane )
{
	*outVec = *inVec - (pPlane->normal) * ( inVec->Dot(pPlane->normal )) * 2.0f ;
	return outVec;
}

 inline Vector3* Vec3ReflectCoord( Vector3* outPt, Vector3* inPt, const Plane* pPlane )
{
	float dist =  pPlane->Dot(*inPt );
	*outPt = *inPt - pPlane->normal * dist * 2.0;
	return outPt;
}

 inline void  Vec3OrthoNormalize(Vector3& normal, Vector3& tangent)
 {
	 normal = normal.Safenormalize();
	 Vector3 proj = normal*tangent.Dot(normal);
	 tangent =  tangent - proj;
	 tangent = tangent.Safenormalize();
 }

 void Vec3GetYawPitchAngles( const Vector3 &vec, float &yawAng, float &pitchAng );
 Vector3* Vec3GetVectorOfYawPitch( Vector3*  vec, float yawAng, float pitchAng );

 void	CreateObliqueNearPlaneClippingMatrix(const Plane* clipPlane,const Matrix44* matrix, Matrix44* outClipMatrix);


//calc reflect plane no normalize 
 Matrix44 * MatrixReflectDirect(
	  Matrix44 *pOut,
	     const Plane *pPlane
	);

	Matrix44*	MatrixClipProj(	Matrix44* pOutProj,	Matrix44* pViewMat,	Matrix44* pOldProjMat,     const Plane *pOldClipPlane);

 int	LineInterscetPoint(const Vector2& start,const Vector2& end,Vector2& pt);

 int	LineInterscetPoint(Vector2& start,Vector2& end,Vector2& pt);

 bool	LineIntersectLine( Vector2 lineStart0, Vector2 lineEnd0, Vector2 lineStart1, Vector2 lineEnd1,float& t,Vector2& outPt);
 
 bool	IsPointInsideTri(const Vector2 *pa,const Vector2 *pb,
	const Vector2 *pc,const Vector2 *p);


 Vector3*  Vec3TransformCoordByQuaternion(Vector3* out, const Vector3* v,const Quaternion* q);
 

 inline  bool IsEqualF( float a, float b, const float epsilon = FLT_EPSILON )
{
	float diff = a - b;
	return diff > -epsilon && diff < epsilon; 
}

 inline bool IsZeroF(const float val, const float epsilon = FLT_EPSILON )
{
	return (val > -epsilon) && (val < epsilon);
}

 inline float ClampFToZero(float& input)
{
	if (input < FLT_EPSILON && input > -FLT_EPSILON)
		input = 0.0f;

	return input;
}

 inline float ClampF(float val, float low, float high)
{
	return (float)std::max(std::min(val, high), low);
}
 
 inline float  Sgn(float a)
 {
	 if (a > 0.0F) return (1.0F);
	 if (a < 0.0F) return (-1.0F);
	 return (0.0F);
 }


#endif