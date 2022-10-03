#include "CommonMath.h"
 
 
float  Sqr(float a) {return a*a;}
//float   ClampF(float a) {return Min(1.0f,Max(0.0f,a));}



float Roundf(float a,float precision)
{
	return std::floor(0.5f+a/precision)*precision;
}


float Interpolate(const float &f0,const float &f1,float alpha)
{
	return f0*(1-alpha) + f1*alpha;
}




int     Argmin(float a[],int n)
{
	int r=0;
	for(int i=1;i<n;i++) 
	{
		if(a[i]<a[r]) 
		{
			r = i;
		}
	}
	return r;
}
 

 Matrix33 operator*( const Matrix33& a, const Matrix33& b )
{ 
	return Matrix33(a.x*b,a.y*b,a.z*b); 
}


 Vector4   operator*( const Vector4&   v, const Matrix44& m )
{
	return Vector4(v.Dot(Vector4(m.x.x,m.y.x,m.z.x,m.w.x)),
					v.Dot(Vector4(m.x.y,m.y.y,m.z.y,m.w.y)),
					v.Dot(Vector4(m.x.z,m.y.z,m.z.z,m.w.z)),
					v.Dot(Vector4(m.x.w,m.y.w,m.z.w,m.w.w)));
}

//Vector4   operator*( const Vector4&   v, const Matrix44* m )
//{
//	return Vector4(v.Dot(Vector4(m->x.x,m->y.x,m->z.x,m->w.x)),
//		v.Dot(Vector4(m->x.y,m->y.y,m->z.y,m->w.y)),
//		v.Dot(Vector4(m->x.z,m->y.z,m->z.z,m->w.z)),
//		v.Dot(Vector4(m->x.w,m->y.w,m->z.w,m->w.w)));
//}


 Matrix44 operator*( const Matrix44& a, const Matrix44& b )
{
	return Matrix44(a.x*b,a.y*b,a.z*b,a.w*b);
}


 Vector3*	Vec3TransformCoord(Vector3* vecOut,const Vector3* vecIn,const Matrix44* matIn)
{
	Vector4	tmp(vecIn->x,vecIn->y,vecIn->z,1);
	tmp = Vector4(tmp.Dot(Vector4(matIn->x.x,matIn->y.x,matIn->z.x,matIn->w.x)),
		tmp.Dot(Vector4(matIn->x.y,matIn->y.y,matIn->z.y,matIn->w.y)),
		tmp.Dot(Vector4(matIn->x.z,matIn->y.z,matIn->z.z,matIn->w.z)),
		tmp.Dot(Vector4(matIn->x.w,matIn->y.w,matIn->z.w,matIn->w.w)));

	if(tmp.w)
		*vecOut = tmp.xyz()/tmp.w;
	else
	{
		vecOut->x = 0;
		vecOut->y = 0;
		vecOut->z = 0;
	}
	return vecOut;
}

 Vector3*	Vec3TransformNormal(Vector3* vecOut,const Vector3* vecIn,const Matrix44* matIn)
{
	Vector4	tmp(vecIn->x,vecIn->y,vecIn->z,0);
	tmp = Vector4(tmp.Dot(Vector4(matIn->x.x,matIn->y.x,matIn->z.x,matIn->w.x)),
		tmp.Dot(Vector4(matIn->x.y,matIn->y.y,matIn->z.y,matIn->w.y)),
		tmp.Dot(Vector4(matIn->x.z,matIn->y.z,matIn->z.z,matIn->w.z)),
		tmp.Dot(Vector4(matIn->x.w,matIn->y.w,matIn->z.w,matIn->w.w)));
 
	*vecOut = tmp.xyz();
	
	return vecOut;
}


 Vector2*	Vec2TransformCoord(Vector2* vecOut,const Vector2* vecIn,const Matrix33* matIn)
{
	Vector3	tmp(vecIn->x,vecIn->y,1);
	tmp = Vector3(tmp.Dot(Vector3(matIn->x.x,matIn->y.x,matIn->z.x)),
		tmp.Dot(Vector3(matIn->x.y,matIn->y.y,matIn->z.y)),
		tmp.Dot(Vector3(matIn->x.z,matIn->y.z,matIn->z.z)));
	if(tmp.z)
	{
		vecOut->x = tmp.x/tmp.z;
		vecOut->y = tmp.y/tmp.z;
	}
	else
	{
		vecOut->x = 0;
		vecOut->y = 0;
	}
	return vecOut;
}

 Vector2*	Vec2TransformNormal(Vector2* vecOut,const Vector2* vecIn,const Matrix33* matIn)
{
	Vector3	tmp(vecIn->x,vecIn->y,0);
	tmp = Vector3(tmp.Dot(Vector3(matIn->x.x,matIn->y.x,matIn->z.x)),
		tmp.Dot(Vector3(matIn->x.y,matIn->y.y,matIn->z.y)),
		tmp.Dot(Vector3(matIn->x.z,matIn->y.z,matIn->z.z)));

	(*vecOut).x = tmp.x/tmp.z;
	(*vecOut).y = tmp.y/tmp.z;
	
	return vecOut;
}

// CAABBBox*	AABBBoxTransform(CAABBBox* vecOut,const CAABBBox* vecIn,const Matrix44* mat)
//{
//	assert(vecOut != vecIn);
//	vecOut->Reset();
//	Vector3 vec;
//	for(int i= 0;i<8;i++)
//	{ 
//		vecIn->GetVecBounds(i,&vec);
//		Vec3TransformCoord(&vec,&vec,mat);
//		vecOut->IncreaseBox(vec);
//	}
//	return vecOut;
//}

 Vector3 Vec3Interpolate(const Vector3 &v0,const Vector3 &v1,float alpha) 
{
	return v0*(1-alpha) + v1*alpha;
}

 Vector3 Vec3VectorMin(const Vector3 &a,const Vector3 &b)
{
	return Vector3(Min(a.x,b.x),Min(a.y,b.y),Min(a.z,b.z));
}
 Vector3 Vec3VectorMax(const Vector3 &a,const Vector3 &b)
{
	return Vector3(Max(a.x,b.x),Max(a.y,b.y),Max(a.z,b.z));
}
 
Vector3	Matrix33::GetColumn(int index) const
{
	const float* pFloatArray0 = &x.x;
	const float* pFloatArray1 = &y.x;
	const float* pFloatArray2 = &z.x;
	return Vector3(pFloatArray0[index],pFloatArray1[index],pFloatArray2[index]);
}

Vector3	Matrix33::GetRow(int index) const
{
	const Vector3* pVector = &x;
	return pVector[index];
}
 

Vector4	Matrix44::GetColumn(int index)  const
{
	const float* pFloatArray0 = &x.x;
	const float* pFloatArray1 = &y.x;
	const float* pFloatArray2 = &z.x;
	const float* pFloatArray3 = &w.x;

	return Vector4(pFloatArray0[index],pFloatArray1[index],pFloatArray2[index],pFloatArray3[index]);
}

Vector4	Matrix44::GetRow(int index)  const
{
	//float* pFloatArray = &x.x;
	//return *(Vector4*)(pFloatArray + sizeof(Vector4)*index);	
	const Vector4* pVector = &x;
	return pVector[index];
}

Matrix33		Matrix44::GetRotationPart()	const
{
	return Matrix33(x.xyz(),y.xyz(),z.xyz());
}

void	Matrix44::SetRotationPart(Matrix33& mat33)
{
	x.x = mat33.x.x;
	x.y = mat33.x.y;
	x.z = mat33.x.z;

	y.x = mat33.y.x;
	y.y = mat33.y.y;
	y.z = mat33.y.z;

	z.x = mat33.z.x;
	z.y = mat33.z.y;
	z.z = mat33.z.z;
}

Vector3		Matrix44::GetPositionPart()
{
	return Vector3(w.x,w.y,w.z);
}

const Vector3		Matrix44::GetPositionPart()	const
{
	return Vector3(w.x,w.y,w.z);
}

void	Matrix44::SetPositionPart(Vector3& vecPos)
{
	w.x = vecPos.x;
	w.y = vecPos.y;
	w.z = vecPos.z;
}

//------------- Plane --------------
Plane::Plane(const Vector3 &ptInPlane,const Vector3 &planeNormal)
{
	normal = planeNormal.Normalize();
	dist = -ptInPlane.Dot(normal); 
}

Plane::Plane(const Vector3 &pt0,const Vector3& pt1,Vector3& pt2)
{
	normal = (pt1 - pt0).Cross(pt2 - pt0);
	normal.SetNormalize();
	dist = - pt0.Dot(normal);
}

void Plane::Transform(const Vector3 &position, const Quaternion &orientation) {
	//   Transforms the plane to the space defined by the 
	//   given position/orientation.
	Vector3 newnormal;
	Vector3 origin;

	newnormal = (orientation.Inverse()*normal).xyz();
	origin = (orientation.Inverse()*(-normal*dist - position)).xyz();

	normal = newnormal;
	dist = -newnormal.Dot(origin);
}

void	Plane::Transform(const Matrix44& transMat)
{
	//   Transforms the plane to the space defined by the 
	//   given transMat
	Matrix44 tmpMatTrans = transMat;
	tmpMatTrans = tmpMatTrans.Inverse();
	tmpMatTrans = tmpMatTrans.Transpose();

	Vector4 v(normal.x,normal.y,normal.z,dist);
	v = v*tmpMatTrans;
	normal.x = v.x;
	normal.y = v.y;
	normal.z = v.z;
	dist = v.w;
}

 

//--------- utility functions -------------


 Quaternion QuaternionSlerp( const Quaternion& q0, const Quaternion& q1, float interp )
{
	Quaternion a = q0;
	if(a.Dot(q1) <0.0) 
	{
		a.w=-a.w;
		a.x=-a.x;
		a.y=-a.y;
		a.z=-a.z;
	}
	float d = a.Dot(q1);
	if(d>=1.0) {
		return a;
	}
	float theta = std::acos(d);
	if(theta==0.0f) { return(a);}
	return a*(std::sin(theta-interp*theta)/std::sin(theta)) + q1*(std::sin(interp*theta)/std::sin(theta));
	
	//float magnitude = sqrt(q0.MagnitudeSq() * q1.MagnitudeSq()); 
	//assert(magnitude > float(0));

	//float product = q0.Dot(q1) / magnitude;
	//if (abs(product) != float(1))
	//{
	//	// Take care of long angle case see http://en.wikipedia.org/wiki/Slerp
	//	const float sign = (product < 0) ? float(-1) : float(1);

	//	const float theta = acos(sign * product);
	//	const float s1 = sin(sign * interp * theta);   
	//	const float d = float(1.0) / sin(theta);
	//	const float s0 = sin((float(1.0) - interp) * theta);

	//	return Quaternion(
	//		(q0.x * s0 + q1.x * s1) * d,
	//		(q0.y * s0 + q1.y * s1) * d,
	//		(q0.z * s0 + q1.z * s1) * d,
	//		(q0.w * s0 + q1.w * s1) * d);
	//}
	//else
	//{
	//	return q0;
	//}
}

 Quaternion QuaternionInterpolate(const Quaternion &q0,const Quaternion &q1,float t) 
{
	//-----------------------------------
	// Calculate the cosine of the angle:
	Quaternion res;
	double cosOmega = q0.Dot( q1 );

	//-----------------------------------
	// adjust signs if necessary:

	float sign2;
	if ( cosOmega < 0.0 )
	{
		cosOmega = -cosOmega;
		sign2 = -1.0f;
	}
	else
		sign2 = 1.0f;

	//-----------------------------------
	// calculate interpolating coeffs:

	double scale1, scale2;
	if ( (1.0 - cosOmega) > 0.00001 )
	{
		// standard case
		double omega = std::acos(cosOmega);
		double sinOmega = std::sin(omega);
		scale1 = std::sin((1.0 - t) * omega) / sinOmega;
		scale2 = sign2 * std::sin(t * omega) / sinOmega;
	}
	else
	{
		// if quats are very close, just do linear interpolation
		scale1 = 1.0 - t;
		scale2 = sign2 * t;
	}


	//-----------------------------------
	// actually do the interpolation:

	res.x = float(scale1 * q0.x + scale2 * q1.x);
	res.y = float(scale1 * q0.y + scale2 * q1.y);
	res.z = float(scale1 * q0.z + scale2 * q1.z);
	res.w = float(scale1 * q0.w + scale2 * q1.w);
	return res;
}

 Quaternion	 QuaternionExtrapolate(const Quaternion &q0,const Quaternion &q1,float t) 
{
	Quaternion res;
	// assert t >= 0 && t <= 1
	// q0 is value at time = 0
	// q1 is value at time = t
	// Computes quaternion at time = 1
	double flip,cos = q0.x * q1.x + q0.y * q1.y + q0.z * q1.z + q0.w * q1.w;
	if (cos < 0.0)
	{
		cos = -cos;
		flip = -1.0;
	}
	else
		flip = 1.0;

	double s1,s2;
	if ((1.0 - cos) > 0.00001)
	{
		double om = std::acos(cos) / t;
		double sd = 1.0 / std::sin(t * om);
		s1 = flip * std::sin(om) * sd;
		s2 = std::sin((1.0 - t) * om) * sd;
	}
	else
	{
		// If quats are very close, do linear interpolation
		s1 = flip / t;
		s2 = (1.0 - t) / t;
	}

	res.x = float(s1 * q1.x - s2 * q0.x);
	res.y = float(s1 * q1.y - s2 * q0.y);
	res.z = float(s1 * q1.z - s2 * q0.z);
	res.w = float(s1 * q1.w - s2 * q0.w);

	return res;
}

	bool			QuaternionIsNan(const Quaternion &q0)
{
	return IsNan(q0.x) ||IsNan(q0.y)|| IsNan(q0.z) ||IsNan(q0.w) ;
}

 Matrix44 MatrixRigidInverse(const Matrix44 &m)
{
	Matrix44 trans_inverse;
	MatrixTranslation(&trans_inverse,-m.w.xyz());
	Matrix44 rot   = m;
	rot.w = Vector4(0,0,0,1);
	return trans_inverse * rot.Transpose();
}

 Matrix44*		MatrixPerspectiveLH(  Matrix44 *pOut,
													        float w,
													        float h,
													        float zn,
													        float zf )
{
	pOut->x.x = 2*zn/w;
	pOut->y.y = 2*zn/h;
	pOut->z.z = zf/(zf- zn);
	pOut->z.w = 1;
	pOut->w.z = zn*zf/(zn - zf);
	pOut->w.w = 0;
	return pOut;	
}

 Matrix44*		MatrixPerspectiveFovLH(  Matrix44 *pOut,
											     float fovy, 
											     float aspect, 
											     float zn, 
											     float zf )
{
	float h = 1.0f/std::tan(fovy/2.0f); // view space height
	float w = h / aspect ;  // view space width
	pOut->SetIndetity();
	pOut->x.x = w;
	pOut->y.y = h;
	pOut->z.z = zf/(zf- zn);
	pOut->z.w = 1;
	pOut->w.z = zn*zf/(zn - zf);
	pOut->w.w = 0;
	return pOut;
	//return Matrix44(
	//	w, 0, 0             ,   0,
	//	0, h, 0             ,   0,
	//	0, 0, zf/(zf-zn)    ,  -1,
	//	0, 0, zn*zf/(zn-zf) ,   0 );
}

 Matrix44*		MatrixPerspectiveFovRH(  Matrix44 *pOut,
											       float fovy, 
											       float aspect, 
											       float zn, 
											       float zf )
{
	float h = 1.0f/ std::tan(fovy/2.0f); // view space height
	float w = h / aspect ;  // view space width
	pOut->SetIndetity();
	pOut->x.x = w;
	pOut->y.y = h;
	pOut->z.z = zf/(zn- zf);
	pOut->z.w = -1;
	pOut->w.z = zn*zf/(zn - zf);
	pOut->w.w = 0;
	return pOut;
	//return Matrix44(
	//	w, 0, 0             ,   0,
	//	0, h, 0             ,   0,
	//	0, 0, zf/(zn-zf)    ,  -1,
	//	0, 0, zn*zf/(zn-zf) ,   0 );
}

 Matrix44*		MatrixLookAtLH(  Matrix44 *pOut,
									     const Vector3* eye, 
									     const Vector3* at, 
									     const Vector3* up)
{

	//zaxis = normal(At - Eye)
	//	xaxis = normal(cross(Up, zaxis))
	//	yaxis = cross(zaxis, xaxis)

	//	xaxis.x           yaxis.x           zaxis.x          0
	//	xaxis.y           yaxis.y           zaxis.y          0
	//	xaxis.z           yaxis.z           zaxis.z          0
	//	-dot(xaxis, eye)  -dot(yaxis, eye)  -dot(zaxis, eye)  1

	pOut->SetIndetity();
	Vector3 zaxis = (*at - *eye).Normalize();
	Vector3 xaxis = up->Cross(zaxis);
	
	//if(xaxis.MagnitudeSq()< FLT_EPSILON)
	//{
	//	assert(false);
	//}
	xaxis.Safenormalize();
	Vector3 yaxis = zaxis.Cross(xaxis);


	pOut->x.x = xaxis.x;
	pOut->x.y = yaxis.x;
	pOut->x.z = zaxis.x;

	pOut->y.x = xaxis.y;
	pOut->y.y = yaxis.y;
	pOut->y.z = zaxis.y;

	pOut->z.x = xaxis.z;
	pOut->z.y = yaxis.z;
	pOut->z.z = zaxis.z;
	
	pOut->w.x = -xaxis.Dot(*eye);
	pOut->w.y = -yaxis.Dot(*eye);
	pOut->w.z = -zaxis.Dot(*eye);
	pOut->w.w = 1.0f;

	return pOut;
}

 Matrix44*		MatrixLookAtRH(  Matrix44 *pOut,
									       const Vector3* eye, 
									       const Vector3* at, 
									       const Vector3* up)
{
	//zaxis = normal(Eye - At)
	//	xaxis = normal(cross(Up, zaxis))
	//	yaxis = cross(zaxis, xaxis)

	//xaxis.x           yaxis.x           zaxis.x          0
	//xaxis.y           yaxis.y           zaxis.y          0
	//xaxis.z           yaxis.z           zaxis.z          0
	//-dot(xaxis, eye)  -dot(yaxis, eye)  -dot(zaxis, eye)  1


	pOut->SetIndetity();
	Vector3 zaxis = (*eye - *at).Normalize();
	Vector3 xaxis = up->Cross(zaxis);

	if(xaxis.MagnitudeSq()< FLT_EPSILON)
	{
		xaxis = Vector3(1,0,0);
	}
	xaxis.SetNormalize();
	Vector3 yaxis = zaxis.Cross(xaxis);

	pOut->w.w = 1.0f;
	pOut->x.x = xaxis.x;
	pOut->x.y = yaxis.x;
	pOut->x.z = zaxis.x;

	pOut->y.x = xaxis.y;
	pOut->y.y = yaxis.y;
	pOut->y.z = zaxis.y;

	pOut->z.x = xaxis.z;
	pOut->z.y = yaxis.z;
	pOut->z.z = zaxis.z;

	pOut->w.x = -xaxis.Dot(*eye);
	pOut->w.y = -yaxis.Dot(*eye);
	pOut->w.z = -zaxis.Dot(*eye);


	return pOut;

}

 Matrix44*		MatrixOrthoLH(
									   Matrix44 *pOut,
									      float w,
									      float h,
									      float zn,
									      float zf
									 )
{
	//2/w  0    0           0
	//	0    2/h  0           0
	//	0    0    1/(zf-zn)   0
	//	0    0   -zn/(zf-zn)  1
	pOut->SetIndetity();
	pOut->x.x = 2.0f/w;
	pOut->y.y = 2.0f/h;
	pOut->z.z = 1/(zf - zn);
	pOut->w.z = -zn/(zf - zn);
	return pOut;
}

 
 Matrix44*		MatrixOrthoRH(
									   Matrix44 *pOut,
									      float w,
									      float h,
									      float zn,
									      float zf
									 )
{
	// 2/w  0    0           0
	//	0    2/h  0           0
	//	0    0    1/(zn-zf)   0
	//	0    0    zn/(zn-zf)  1

	pOut->SetIndetity();
	pOut->x.x = 2.0f/w;
	pOut->y.y = 2.0f/h;
	pOut->z.z = 1/(zn - zf);
	pOut->w.z =  zn/(zn - zf);
	return pOut;
}

 Matrix44*		MatrixLookViewDir(
									    Matrix44 *pOut,
									       Vector3* pVecEye,
									       Vector3* pVecDir,
									  		Vector3*	pVecUp)
{ 
	Vector3 xaxis = pVecDir->Cross(*pVecUp);
	xaxis.SetNormalize();
	Vector3 yaxis = xaxis.Cross(*pVecDir);
	yaxis.SetNormalize();

	pOut->w.w = 1.0f;
	pOut->x.x = xaxis.x;
	pOut->x.y = yaxis.x;
	pOut->x.z = pVecDir->x;

	pOut->y.x = xaxis.y;
	pOut->y.y = yaxis.y;
	pOut->y.z = pVecDir->y;

	pOut->z.x = xaxis.z;
	pOut->z.y = yaxis.z;
	pOut->z.z = pVecDir->z;

	pOut->w.x = -xaxis.Dot(*pVecEye);
	pOut->w.y = -yaxis.Dot(*pVecEye);
	pOut->w.z = -pVecDir->Dot(*pVecEye);

	return pOut;

}
 
 Matrix33*		Matrix33RotateAngle( Matrix33* pOutMat,const float angle_radians)
{
	pOutMat->x.x = cos(angle_radians);
	pOutMat->x.y = sin(angle_radians);
	pOutMat->x.z = 0;
	pOutMat->y.x = -sin(angle_radians);
	pOutMat->y.y = cos(angle_radians);
	pOutMat->y.z = 0;
	pOutMat->z.x = 0;
	pOutMat->z.y = 0;
	pOutMat->z.z = 1;
	return pOutMat;
}

 Matrix44*		MatrixTranslation( Matrix44* pOutMat, const Vector3 &t)
{
	pOutMat->SetIndetity();
	pOutMat->w.x = t.x;
	pOutMat->w.y = t.y;
	pOutMat->w.z = t.z;
	return pOutMat;
}
 Matrix44*		MatrixTranslation( Matrix44* pOutMat, float x, float y, float z)
{
	pOutMat->SetIndetity();
	pOutMat->w.x = x;
	pOutMat->w.y = y;
	pOutMat->w.z = z;
	return pOutMat;
}
 Matrix44*		MatrixScaling( Matrix44* pOutMat,  float x, float y, float z)
{
	pOutMat->SetIndetity();
	pOutMat->x.x = x;
	pOutMat->y.y = y;
	pOutMat->z.z = z;
	return pOutMat;
}

 Matrix44*		MatrixRotationX( Matrix44* pOutMat, const float angle_radians)
{
	float s =  std::sin(angle_radians);
	float c =  std::cos(angle_radians);

	pOutMat->SetIndetity();
	pOutMat->y.y = c;
	pOutMat->y.z = s;
	pOutMat->z.y = -s;
	pOutMat->z.z = c;
	return pOutMat;
}
 Matrix44*		MatrixRotationY( Matrix44* pOutMat, const float angle_radians)
{
	float s =  std::sin(angle_radians);
	float c =  std::cos(angle_radians);
	pOutMat->SetIndetity();
	pOutMat->x.x = c;
	pOutMat->z.z = c;
	pOutMat->x.z = -s;
	pOutMat->z.x = s;
	return pOutMat;

}
 Matrix44*		MatrixRotationZ( Matrix44* pOutMat, const float angle_radians)
{
	float s =  std::sin(angle_radians);
	float c =  std::cos(angle_radians);
	pOutMat->SetIndetity();
	pOutMat->x.x = c;
	pOutMat->y.y = c;
	pOutMat->x.y = s;
	pOutMat->y.x = -s;
	return pOutMat;
}

 Matrix44*		MatrixMultiply( Matrix44* pOutMat, const Matrix44* pInMat0, const Matrix44* pInMat1)
{
	*pOutMat = *pInMat0* (*pInMat1);
	return pOutMat;
}
 
 Matrix44*		MatrixRotationYawPitchRoll(  Matrix44 *pOut,
												       float Yaw,
												       float Pitch,
												       float Roll)
{
	Quaternion q;
	q.SetYawPitchRollInRadian(Yaw,Pitch,Roll);
	pOut->SetIndetity();
	Matrix33 rot = q.Getmatrix();
	pOut->SetRotationPart(rot) ;
	return pOut;
}

 Matrix44*		MatrixRotationQuaternion(  Matrix44 *pOut,
															     const Quaternion* pRot)
{
	Quaternion qtNormalRot = pRot->Normalize();
	pOut->SetIndetity();
	Matrix33 rot = qtNormalRot.Getmatrix();
	pOut->SetRotationPart(rot) ;
	return pOut;
}

 Quaternion* QuaternionRotationMatrix( Quaternion *pOut, const Matrix44 *pM)
{
	Matrix33 rotPart = pM->GetRotationPart();
	return QuaternionRotationMatrix(pOut,&rotPart);
}

 Quaternion* QuaternionRotationMatrix( Quaternion *pOut, const Matrix33 *pM)
{
	float trace = pM->x.x + pM->y.y + pM->z.z;

	float temp[4];

	if (trace > float(0.0)) 
	{
		float s = sqrt(trace + float(1.0));
		temp[3]=(s * float(0.5));
		s = float(0.5) / s;

		temp[0]=((pM->y.z - pM->z.y) * s);
		temp[1]=((pM->z.x - pM->x.z) * s);
		temp[2]=((pM->x.y - pM->y.x) * s);
	} 
	else 
	{
		int i = pM->x.x < pM->y.y ? 
			(pM->y.y< pM->z.z ? 2 : 1) :
			(pM->x.x < pM->z.z ? 2 : 0); 
		int j = (i + 1) % 3;  
		int k = (i + 2) % 3;

		float s = sqrt((*pM)(i,i) - (*pM)(j,j) - (*pM)(k,k) + float(1.0));
		temp[i] = s * float(0.5);
		s = float(0.5) / s;

		temp[3] = ((*pM)(j,k) - (*pM)(k,j)) * s;
		temp[j] = ((*pM)(i,j) + (*pM)(j,i)) * s;
		temp[k] = ((*pM)(i,k) + (*pM)(k,i)) * s;
	}
	pOut->x = temp[0];
	pOut->y = temp[1];
	pOut->z = temp[2];
	pOut->w = temp[3];
	
	return pOut;
}

 Matrix44 * MatrixRotationAxis(
	  Matrix44 *pOut,
	     const Vector3 *pV,
	     float Angle
	)
{
	Quaternion rot(*pV,Angle);
	pOut->SetIndetity();
	Matrix33 rotMat = rot.Getmatrix();
	pOut->SetRotationPart(rotMat);
	return pOut;
}

  void MatrixDecompose(
									    Vector3 *pOut,
									    Quaternion *pOutRotation,
									    Vector3 *pOutTranslation,
									    const Matrix44 *pM
									  )
{
	pOut->x = pM->x.xyz().Magnitude();
	pOut->y = pM->y.xyz().Magnitude();
	pOut->z = pM->z.xyz().Magnitude();

	*pOutTranslation = pM->w.xyz();
	Matrix33 tmpMat33 = Matrix33(pM->x.xyz()*(1.0f/pOut->x),
		pM->y.xyz()*(1.0f/pOut->y),
		pM->z.xyz()*(1.0f/pOut->z));
	QuaternionRotationMatrix(pOutRotation,&tmpMat33);
	
}

	Matrix44 * MatrixReflect(
	  Matrix44 *pOut,
	     const Plane *pPlane
	)
{

	//P = normalize(Plane);

	//-2 * P.a * P.a + 1  -2 * P.b * P.a      -2 * P.c * P.a        0
	//	-2 * P.a * P.b      -2 * P.b * P.b + 1  -2 * P.c * P.b        0
	//	-2 * P.a * P.c      -2 * P.b * P.c      -2 * P.c * P.c + 1    0
	//	-2 * P.a * P.d      -2 * P.b * P.d      -2 * P.c * P.d        1
	Plane P = pPlane->Normalize();
	
	pOut->x.x = -2 * P.normal.x * P.normal.x + 1;
	pOut->x.y = -2 * P.normal.y * P.normal.x;
	pOut->x.z = -2 * P.normal.z * P.normal.x;
	pOut->x.w = 0;
	pOut->y.x = -2 * P.normal.x * P.normal.y;
	pOut->y.y = -2 * P.normal.y * P.normal.y + 1;
	pOut->y.z = -2 * P.normal.z * P.normal.y;
	pOut->y.w = 0;
	pOut->z.x = -2 * P.normal.x * P.normal.z;
	pOut->z.y = -2 * P.normal.y * P.normal.z;
	pOut->z.z = -2 * P.normal.z * P.normal.z + 1;
	pOut->z.w = 0;
	pOut->w.x = -2 * P.normal.x * P.dist;
	pOut->w.y = -2 * P.normal.y * P.dist;
	pOut->w.z = -2 * P.normal.z * P.dist;
	pOut->w.w = 1;
	return pOut;
}


	Plane * PlaneTransform(
	  Plane *pOut,
	     const Plane *pP,
	     const Matrix44 *pM
	)
{
	pOut->normal.x = (*pM)(0,0) * pP->normal.x + (*pM)(1,0) * pP->normal.y + (*pM)(2,0) * pP->normal.z + (*pM)(3,0) * pP->dist;
	pOut->normal.y = (*pM)(0,1) * pP->normal.x + (*pM)(1,1) * pP->normal.y + (*pM)(2,1) * pP->normal.z + (*pM)(3,1) * pP->dist;
	pOut->normal.z = (*pM)(0,2) * pP->normal.x + (*pM)(1,2) * pP->normal.y + (*pM)(2,2) * pP->normal.z + (*pM)(3,2) * pP->dist;
	pOut->dist = (*pM)(0,3) * pP->normal.x + (*pM)(1,3) * pP->normal.y + (*pM)(2,3) * pP->normal.z + (*pM)(3,3) * pP->dist;
	return pOut;

	//Vector4 vecPlane(pP->normal.x,pP->normal.y,pP->normal.z,pP->dist);
	//vecPlane = vecPlane*(*pM);
	//pOut->normal.x = vecPlane.x;
	//pOut->normal.y = vecPlane.y;
	//pOut->normal.z = vecPlane.z;
	//pOut->dist = vecPlane.w;
	//return pOut;
}

//        RotationArc()
// Given two vectors v0 and v1 this function
// returns quaternion q where q*v0==v1.
// Routine taken from game programming gems.
 Quaternion RotationArc(Vector3 v0,Vector3 v1){
	Quaternion q;
	v0 = v0.Normalize();  // Comment these two lines out if you know its not needed.
	v1 = v1.Normalize();  // If vector is already unit length then why do it again?
	Vector3  c = v0.Cross(v1);
	float   d = v0.Dot(v1);
	if(d<=-1.0f) { return Quaternion(1,0,0,0);} // 180 about x axis
	float   s = std::sqrt((1+d)*2);
	q.x = c.x / s;
	q.y = c.y / s;
	q.z = c.z / s;
	q.w = s /2.0f;
	return q;
}


 Matrix44 MatrixFromQuatVec(const Quaternion &q, const Vector3 &v) 
{
	// builds a 4x4 transformation matrix based on orientation q and translation v 
	float qx2 = q.x*q.x;
	float qy2 = q.y*q.y;
	float qz2 = q.z*q.z;

	float qxqy = q.x*q.y;
	float qxqz = q.x*q.z;
	float qxqw = q.x*q.w;
	float qyqz = q.y*q.z;
	float qyqw = q.y*q.w;
	float qzqw = q.z*q.w;

	return Matrix44(
		1-2*(qy2+qz2),  
		2*(qxqy+qzqw),
		2*(qxqz-qyqw),  
		0            ,  
		2*(qxqy-qzqw),  
		1-2*(qx2+qz2),
		2*(qyqz+qxqw),  
		0            ,  
		2*(qxqz+qyqw),  
		2*(qyqz-qxqw),  
		1-2*(qx2+qy2),  
		0    , 
		v.x ,
		v.y ,
		v.z ,
		1.0f );
}

//可能有问题，要测试
#define PLANEDIREPSILON 0.00000001f
 float PlaneIntersectRay( Vector3* ptOut,
													 Plane* plane, 
													 const Vector3* rayDir, 
													 const Vector3*  rayOrigin )
{
	//const float t = (-plane->d - D3DXVec3Dot(&normal, rayOrigin)) / D3DXVec3Dot(&normal, rayDir);
	float fDotdir = plane->normal.Dot(*rayDir); 
 
	float fDis = plane->Dot(*rayOrigin); 

	float t= -fDis/fDotdir;

	if(t < 0)
		return t;	
	*ptOut = *rayOrigin + t*(*rayDir);
	return t;
 

}



 Vector3 PlaneLineIntersection(const Plane &plane, const Vector3 &p0, const Vector3 &p1)
{
	// returns the point where the line p0-p1 intersects the plane n&d
	Vector3 dif;
	dif = p1-p0;
	float dn= plane.normal.Dot(dif);
	float t = -(plane.dist+plane.normal.Dot(p0) )/dn;
	return p0 + (dif*t);
}

 Vector3 PlaneProject(const Plane &plane, const Vector3 &point)
{
	return point - plane.normal * (point.Dot(plane.normal)+plane.dist);
}

 Vector3 LineProject(const Vector3 &p0, const Vector3 &p1, const Vector3 &a)
{
	Vector3 w;
	w = p1-p0;
	float t= w.Dot((a-p0)) / (Sqr(w.x)+Sqr(w.y)+Sqr(w.z));
	return p0+ w*t;
}


 float LineProjectTime(const Vector3 &p0, const Vector3 &p1, const Vector3 &a)
{
	Vector3 w;
	w = p1-p0;
	float t= w.Dot((a-p0)) / (Sqr(w.x)+Sqr(w.y)+Sqr(w.z));
	return t;
}



 Vector3 TriNormal(const Vector3 &v0, const Vector3 &v1, const Vector3 &v2)
{
	// return the normal of the triangle
	// inscribed by v0, v1, and v2
	Vector3 cp= v1-v0.Cross(v2-v1);
	float m= cp.Magnitude();
	if(m==0) return Vector3(1,0,0);
	return cp*(1.0f/m);
}



 int BoxInside(const Vector3 &p, const Vector3 &bmin, const Vector3 &bmax) 
{
	return (p.x >= bmin.x && p.x <=bmax.x && 
		p.y >= bmin.y && p.y <=bmax.y && 
		p.z >= bmin.z && p.z <=bmax.z );
}


 int BoxIntersect(const Vector3 &v0, const Vector3 &v1, const Vector3 &bmin, const Vector3 &bmax,Vector3 *impact)
{
	if(BoxInside(v0,bmin,bmax))
	{
		*impact=v0;
		return 1;
	}
	if(v0.x<=bmin.x && v1.x>=bmin.x) 
	{
		float a = (bmin.x-v0.x)/(v1.x-v0.x);
		//v.x = bmin.x;
		float vy =  (1-a) *v0.y + a*v1.y;
		float vz =  (1-a) *v0.z + a*v1.z;
		if(vy>=bmin.y && vy<=bmax.y && vz>=bmin.z && vz<=bmax.z) 
		{
			impact->x = bmin.x;
			impact->y = vy;
			impact->z = vz;
			return 1;
		}
	}
	else if(v0.x >= bmax.x  &&  v1.x <= bmax.x) 
	{
		float a = (bmax.x-v0.x)/(v1.x-v0.x);
		//v.x = bmax.x;
		float vy =  (1-a) *v0.y + a*v1.y;
		float vz =  (1-a) *v0.z + a*v1.z;
		if(vy>=bmin.y && vy<=bmax.y && vz>=bmin.z && vz<=bmax.z) 
		{
			impact->x = bmax.x;
			impact->y = vy;
			impact->z = vz;
			return 1;
		}
	}
	if(v0.y<=bmin.y && v1.y>=bmin.y) 
	{
		float a = (bmin.y-v0.y)/(v1.y-v0.y);
		float vx =  (1-a) *v0.x + a*v1.x;
		//v.y = bmin.y;
		float vz =  (1-a) *v0.z + a*v1.z;
		if(vx>=bmin.x && vx<=bmax.x && vz>=bmin.z && vz<=bmax.z) 
		{
			impact->x = vx;
			impact->y = bmin.y;
			impact->z = vz;
			return 1;
		}
	}
	else if(v0.y >= bmax.y  &&  v1.y <= bmax.y) 
	{
		float a = (bmax.y-v0.y)/(v1.y-v0.y);
		float vx =  (1-a) *v0.x + a*v1.x;
		// vy = bmax.y;
		float vz =  (1-a) *v0.z + a*v1.z;
		if(vx>=bmin.x && vx<=bmax.x && vz>=bmin.z && vz<=bmax.z)
		{
			impact->x = vx;
			impact->y = bmax.y;
			impact->z = vz;
			return 1;
		}
	}
	if(v0.z<=bmin.z && v1.z>=bmin.z) 
	{
		float a = (bmin.z-v0.z)/(v1.z-v0.z);
		float vx =  (1-a) *v0.x + a*v1.x;
		float vy =  (1-a) *v0.y + a*v1.y;
		// v.z = bmin.z;
		if(vy>=bmin.y && vy<=bmax.y && vx>=bmin.x && vx<=bmax.x) 
		{
			impact->x = vx;
			impact->y = vy;
			impact->z = bmin.z;
			return 1;
		}
	}
	else if(v0.z >= bmax.z  &&  v1.z <= bmax.z) 
	{
		float a = (bmax.z-v0.z)/(v1.z-v0.z);
		float vx =  (1-a) *v0.x + a*v1.x;
		float vy =  (1-a) *v0.y + a*v1.y;
		// v.z = bmax.z;
		if(vy>=bmin.y && vy<=bmax.y && vx>=bmin.x && vx<=bmax.x) 
		{
			impact->x = vx;
			impact->y = vy;
			impact->z = bmax.z;
			return 1;
		}
	}
	return 0;
}


 float DistanceBetweenLines(const Vector3 &ustart, const Vector3 &udir, const Vector3 &vstart, const Vector3 &vdir, Vector3 *upoint, Vector3 *vpoint)
{
	Vector3 cp;
	cp = udir.Cross(vdir).Normalize();

	float distu = -cp.Dot(ustart);
	float distv = -cp.Dot(vstart);
	float dist = (float)fabs(distu-distv);
	if(upoint) 
	{
		Plane plane;
		plane.normal = vdir.Cross(cp).Normalize();
		plane.dist = -plane.normal.Dot(vstart);
		*upoint = PlaneLineIntersection(plane,ustart,ustart+udir);
	}
	if(vpoint) 
	{
		Plane plane;
		plane.normal = udir.Cross(cp).Normalize();
		plane.dist = -plane.normal.Dot(ustart);
		*vpoint = PlaneLineIntersection(plane,vstart,vstart+vdir);
	}
	return dist;
}


 Quaternion VirtualTrackBall(const Vector3 &cop, const Vector3 &cor, const Vector3 &dir1, const Vector3 &dir2) 
{
	// routine taken from game programming gems.
	// Implement track ball functionality to spin stuf on the screen
	//  cop   center of projection
	//  cor   center of rotation
	//  dir1  old mouse direction 
	//  dir2  new mouse direction
	// pretend there is a sphere around cor.  Then find the points
	// where dir1 and dir2 intersect that sphere.  Find the
	// rotation that takes the first point to the second.
	float m;
	// compute plane 
	Vector3 nrml = cor - cop;
	float fudgefactor = 1.0f/(nrml.Magnitude() * 0.25f); // since trackball proportional to distance from cop
	nrml = nrml.Normalize();
	float dist = -nrml.Dot(cor);
	Vector3 u= PlaneLineIntersection(Plane(nrml,dist),cop,cop+dir1);
	u=u-cor;
	u=u*fudgefactor;
	m= u.Magnitude();
	if(m>1)
	{
		u/=m;
	}
	else 
	{
		u=u - (nrml * std::sqrt(1-m*m));
	}
	Vector3 v= PlaneLineIntersection(Plane(nrml,dist),cop,cop+dir2);
	v=v-cor;
	v=v*fudgefactor;
	m= v.Magnitude();
	if(m>1) 
	{
		v/=m;
	}
	else 
	{
		v=v - (nrml * std::sqrt(1-m*m));
	}
	return RotationArc(u,v);
}


int countpolyhit=0;
 int PolyHit(const Vector3 *vert, const int n, const Vector3 &v0, const Vector3 &v1, Vector3 *impact, Vector3 *normal)
{
	countpolyhit++;
	int i;
	Vector3 nrml(0,0,0);
	for(i=0;i<n;i++) 
	{
		int i1=(i+1)%n;
		int i2=(i+2)%n;
		nrml = nrml + (vert[i1]-vert[i]).Cross(vert[i2]-vert[i1]);
	}

	float m = nrml.Magnitude();
	if(m==0.0)
	{
		return 0;
	}
	nrml = nrml * (1.0f/m);
	float dist = -nrml.Dot(vert[0]);
	float d0,d1;
	if((d0=v0.Dot(nrml)+dist) <0  ||  (d1=v1.Dot(nrml)+dist) >0) 
	{        
		return 0;
	}

	Vector3 the_point; 
	// By using the cached plane distances d0 and d1
	// we can optimize the following:
	//     the_point = planelineintersection(nrml,dist,v0,v1);
	float a = d0/(d0-d1);
	the_point = v0*(1-a) + v1*a;


	int inside=1;
	for(int j=0;inside && j<n;j++) 
	{
		// let inside = 0 if outside
		Vector3 pp1,pp2,side;
		pp1 = vert[j] ;
		pp2 = vert[(j+1)%n];
		side = (pp2-pp1).Cross((the_point-pp1));
		inside = (nrml.Dot(side) >= 0.0);
	}
	if(inside) 
	{
		if(normal){*normal=nrml;}
		if(impact){*impact=the_point;}
	}
	return inside;
}


 float GetSignedAngleBetweenVectors(Vector3* vec0, Vector3* vec1,const Vector3* vecRightOfVecs)
{
	Vector3 tmpvec0 = vec0->Safenormalize();
	Vector3 tmpvec1 = vec1->Safenormalize();

	float fDot = tmpvec0.Dot(tmpvec1);
	Vector3 vecCross = tmpvec0.Cross(tmpvec1); 

	float fRightDot = vecCross.Dot(*vecRightOfVecs); 

	if(fDot > 1.0 - FLT_EPSILON) 
		return 0.0f;
	if(fDot< -1.0 + FLT_EPSILON)
		return FLT_PI;

	float angleBetween = std::acos(fDot);
	if(fDot < -1 || fDot > 1)
		angleBetween = 0.0f;

	if (fRightDot < 0.0f)
		angleBetween  = -angleBetween;

	return angleBetween;
}

 
 Vector3  GetClosestPointToRefPt(Vector3* ptVec,unsigned int ptNum,const Vector3& refPt)  
{
	Vector3 closest;
	float fMinLengthToRefPt = (ptVec[0].x - refPt.x)*(ptVec[0].x - refPt.x) +
		(ptVec[0].y - refPt.y)*(ptVec[0].y - refPt.y) +
		(ptVec[0].z - refPt.z)*(ptVec[0].z - refPt.z);
	float fLengthToRefPt = fMinLengthToRefPt;
	closest = ptVec[0];
	for(unsigned int i= 1;i<ptNum;i++)
	{
		fLengthToRefPt = (ptVec[0].x - refPt.x)*(ptVec[0].x - refPt.x) +
			(ptVec[0].y - refPt.y)*(ptVec[0].y - refPt.y) +
			(ptVec[0].z - refPt.z)*(ptVec[0].z - refPt.z);
		if(fLengthToRefPt < fMinLengthToRefPt)
		{
			closest = ptVec[i];
			fMinLengthToRefPt = fLengthToRefPt;
		}
	}

	return closest;
}

 Matrix44	GetTransformMatrixOfObj(Vector3 ObjXAixe,
													   Vector3 ObjYAixe,
													   Vector3 ObjZAixe)
{
	Matrix44 mat;
	mat.SetIndetity();
	mat.x.x = ObjXAixe.x ;
	mat.x.y = ObjXAixe.y ;
	mat.x.z = ObjXAixe.z ;

	mat.y.x = ObjYAixe.x ;
	mat.y.y = ObjYAixe.y ;
	mat.y.z = ObjYAixe.z ;

	mat.z.x = ObjZAixe.x ;
	mat.z.y = ObjZAixe.y ;
	mat.z.z = ObjZAixe.z ;

	return mat;
}

 Vector3*	Vec3GetEulerFromQuaternion(Vector3* out,const Quaternion* q)
{
	//// Extract sin(pitch)

	float sp = -2.0f * (q->y*q->z - q->w*q->x);

	// Check for Gimbel lock, giving slight tolerance for numerical imprecision

	if (fabs(sp) > 0.9999f) {

		// Looking straight up or down

		out->y = FLT_PI*0.5f * sp;

		// Compute heading, slam bank to zero

		out->x  = std::atan2(-q->x*q->z + q->w*q->y, 0.5f - q->y*q->y - q->z*q->z);
		out->z  = 0.0f;

	} else {

		// Compute angles.  We don't have to use the "safe" asin
		// function because we already checked for range errors when
		// checking for Gimbel lock

		out->y	= std::asin(sp);
		out->x	= std::atan2(q->x*q->z + q->w*q->y, 0.5f - q->x*q->x - q->y*q->y);
		out->z 	= std::atan2(q->x*q->y + q->w*q->z, 0.5f - q->x*q->x - q->z*q->z);
	}
 
	return out;
}



 Quaternion*	QuaternionDecompress(Quaternion* out,short s0,short s1,short s2)
{
	int which = ((s1 & 1) << 1) | (s2 & 1);
	s1 &= 0xfffe;
	s2 &= 0xfffe;

	static const float scale = 1.0f / 32767.0f / 1.41421f;

	if (which == 3) {
		out->x = s0 * scale;
		out->y = s1 * scale;
		out->z = s2 * scale;

		out->w = 1 - (out->x*out->x) - (out->y*out->y) - (out->z*out->z);
		if (out->w > FLT_EPSILON)
			out->w = sqrt(out->w);
	}
	else if (which == 2) {
		out->x = s0 * scale;
		out->y = s1 * scale;
		out->w = s2 * scale;

		out->z = 1 - (out->x*out->x) - (out->y*out->y) - (out->w*out->w);
		if (out->z > FLT_EPSILON)
			out->z = sqrt(out->z);
	}
	else if (which == 1) {
		out->x = s0 * scale;
		out->z = s1 * scale;
		out->w = s2 * scale;

		out->y = 1 - (out->x*out->x) - (out->z*out->z) - (out->w*out->w);
		if (out->y > FLT_EPSILON)
			out->y = sqrt(out->y);
	}
	else {
		out->y = s0 * scale;
		out->z = s1 * scale;
		out->w = s2 * scale;

		out->x = 1 - (out->y*out->y) - (out->z*out->z) - (out->w*out->w);
		if (out->x > FLT_EPSILON)
			out->x = sqrt(out->x);
	}
	return  out;
}


 Vector3*	CalcPtByUVInFace(Vector3* pOut,Vector3* p0,Vector3* p1,Vector3* p2,float U,float V)
{
	*pOut = *p0 + U*(*p1 - *p0) + V*(*p2 - *p0);
	return pOut;
}



// bool	IsPointInsideTri(const Vector3 *pa,const Vector3 *pb,
//											  const Vector3 *pc,const Vector3 *p,float *pU,float* pV)
//{
//	Vector3 f(pb->x -pa->x ,pb->y - pa->y ,pb->z - pa->z);
//	Vector3 g(pc->x -pa->x ,pc->y - pa->y ,pc->z - pa->z);
//
//	const float a = f.Dot(f);
//	const float b = f.Dot(g);
//	const float c = g.Dot(g);
//
//	const Vector3 vp(p->x - pa->x,p->y - pa->y ,p->z - pa->z);
//	const float d = vp.Dot(f);
//	const float e = vp.Dot(g);
//
//	float x = (d*c) - (e*b);
//	float y = (e*a) - (d*b);
//	const float ac_bb  = (a*c) - (b*b);
//	assert(abs(ac_bb) > FLT_EPSILON);   //all point in the same line
//
//	float inv_ac_bb = 1.0f/ac_bb;
//	*pU= x*inv_ac_bb;
//	*pV= y*inv_ac_bb;
//
//	return ( *pU>= 0) && (*pV >= 0) && (*pU+*pV <= 1);
//}

 bool  LineIntersectTri(const Vector3 *pa,const Vector3 *pb,
											 const Vector3 *pc,const Vector3 *p1,
											 const Vector3 *p2,Vector3 *p,float *pU,float* pV)
{
	double d;  //face plane dist
	double denom,mu;
	Vector3 n;

	/*calc the parameters for the plane*/
	/* Calculate the parameters for the plane */
	n.x = (pb->y - pa->y)*(pc->z - pa->z) - (pb->z - pa->z)*(pc->y - pa->y);
	n.y = (pb->z - pa->z)*(pc->x - pa->x) - (pb->x - pa->x)*(pc->z - pa->z);
	n.z = (pb->x - pa->x)*(pc->y - pa->y) - (pb->y - pa->y)*(pc->x - pa->x);
	n = n.Normalize();
	d = - n.x * pa->x - n.y * pa->y - n.z * pa->z;

	/* Calculate the position on the line that intersects the plane */
	denom = n.x * (p2->x - p1->x) + n.y * (p2->y - p1->y) + n.z * (p2->z - p1->z);
	if (std::abs(denom) < FLT_EPSILON)         /* Line and plane don't intersect */
	{
		//line Parallel the tri,or on the plane of tri

		return false;
	}
	mu = - (d + n.x * p1->x + n.y * p1->y + n.z * p1->z) / denom;	
	if (mu <= FLT_EPSILON || mu > 1)   /* Intersection not along line segment */
		return false;
	p->x = p1->x + mu * (p2->x - p1->x);
	p->y = p1->y + mu * (p2->y - p1->y);
	p->z = p1->z + mu * (p2->z - p1->z);

	return IsPointInsideTri(pa,pb,pc,p,pU,pV);
}

 bool  RayIntersectTri(const Vector3 *pa,const Vector3 *pb,
									  const Vector3 *pc,const Vector3 *RayOrigin,
									  const Vector3 *Raydir,float *pU,float* pV,float* pDist)
{
	Vector3 e1 = *pb - *pa;
	Vector3 e2 = *pc - *pa;
	Vector3 h = Raydir->Cross(e2);
	float a = e1.Dot(h);
 
	if (a > -0.00001 && a < 0.00001)
		return(false);

	float f = 1/a;
	Vector3 s = *RayOrigin - *pa;
	float u = f*s.Dot(h);
 
	if (u < 0.0 || u > 1.0)
		return(false);
	
	Vector3 q = s.Cross(e1);
	float v = f * (*Raydir).Dot(q);
	if (v < 0.0 || u + v > 1.0)
		return(false);
	// at this stage we can compute t to find out where
	// the intersection point is on the line
	*pDist = f * e2.Dot(q);
	if (*pDist > 0) // ray intersection
	{
		Vector3 intersectPt = *RayOrigin + *Raydir*(*pDist);
		return IsPointInsideTri(pa,pb,pc,&intersectPt,pU,pV);
	}
	else // this means that there is a line intersection
		// but not a ray intersection
		return (false); 


}

 bool	IsPointInsideTri(const Vector3 *pa,const Vector3 *pb,
											  const Vector3 *pc,const Vector3 *p,float *pU,float* pV)
{
	Vector3 f(pb->x -pa->x ,pb->y - pa->y ,pb->z - pa->z);
	Vector3 g(pc->x -pa->x ,pc->y - pa->y ,pc->z - pa->z);

	const float a = f.Dot(f);
	const float b = f.Dot(g);
	const float c = g.Dot(g);

	const Vector3 vp(p->x - pa->x,p->y - pa->y ,p->z - pa->z);
	const float d = vp.Dot(f);
	const float e = vp.Dot(g);

	float x = (d*c) - (e*b);
	float y = (e*a) - (d*b);
	const float ac_bb  = (a*c) - (b*b);
	assert(std::abs(ac_bb) > FLT_EPSILON);   //all point in the same line

	float inv_ac_bb = 1.0f/ac_bb;
	*pU= x*inv_ac_bb;
	*pV= y*inv_ac_bb;

	return ( *pU>= 0) && (*pV >= 0) && (*pU+*pV <= 1);
}

 void Vec3GetYawPitchAngles( const Vector3 &vec, float &yawAng, float &pitchAng )
{
	yawAng = std::atan2( vec.x, vec.z );
	if( yawAng < 0.0f )
		yawAng += FLT_2PI;

	if(std::abs(vec.x) > std::abs(vec.z) )
		pitchAng = std::atan2(std::abs(vec.y), std::abs(vec.x) );
	else
		pitchAng = std::atan2(std::abs(vec.y), std::abs(vec.z) );
	if( vec.y < 0.0f )
		pitchAng = -pitchAng;
}

Vector3* Vec3GetVectorOfYawPitch( Vector3*  vec, float yawAng, float pitchAng )
{
	Vector3  pnt( 0.0f, 0.0f, 1.0f );

	Matrix44 matRot;
	MatrixRotationYawPitchRoll(&matRot,yawAng,pitchAng,0);
	Vec3TransformCoord(vec,&pnt,&matRot);
 
	return vec;
}


 void CreateObliqueNearPlaneClippingMatrix(const Plane* pClipPlane,const Matrix44* matrix, Matrix44* outClipMatrix)
{
	Vector4	q;

	// Calculate the clip-space corner point opposite the clipping plane
	// as (sgn(clipPlane.x), sgn(clipPlane.y), 1, 1) and
	// transform it into camera space by multiplying it
	// by the inverse of the projection matrix

	q.x = Sgn(pClipPlane->normal.x) / (*matrix)(0,0);
	q.y = Sgn(pClipPlane->normal.y) / (*matrix)(0,0);
	q.z = 1.0F;
	q.w = (1.0F -(*matrix)(1,1)) / (*matrix)(2,1);

	// Calculate the scaled plane vector
	//	dot = (a*x + b*y + c*z + d*w)
	Plane c = (*pClipPlane) * (1.0F / pClipPlane->Dot(q)); //pClipPlane,&q));

	// Replace the third column of the projection matrix
	*outClipMatrix = *matrix;
	(*outClipMatrix)(0,2) = c.normal.x;
	(*outClipMatrix)(0,1) = c.normal.y;
	(*outClipMatrix)(1,1) = c.normal.z;
	(*outClipMatrix)(2,1) = c.dist;
}


 Matrix44 * MatrixReflectDirect(
	  Matrix44 *pOut,
	     const Plane *pPlane
	)
{
	(*pOut)(0,0) = -2 * pPlane->normal.x * pPlane->normal.x + 1; (*pOut)(0,1)  = -2 * pPlane->normal.y * pPlane->normal.x ;		(*pOut)(0,2) = -2 * pPlane->normal.z * pPlane->normal.x ;	(*pOut)(0,3) =  0;
	(*pOut)(1,0) = -2 * pPlane->normal.x * pPlane->normal.y  ;	 (*pOut)(1,1)  = -2 * pPlane->normal.y * pPlane->normal.y  + 1;	(*pOut)(1,2) = -2 * pPlane->normal.z * pPlane->normal.y  ;	(*pOut)(1,3) =  0;
	(*pOut)(2,0) = -2 * pPlane->normal.x * pPlane->normal.z ;	 (*pOut)(2,1)  = -2 * pPlane->normal.y * pPlane->normal.z ;		(*pOut)(2,2) = -2 * pPlane->normal.z * pPlane->normal.z + 1;(*pOut)(2,3) =  0;
	(*pOut)(3,0) = -2 * pPlane->normal.x * pPlane->dist ;	     (*pOut)(3,1)  = -2 * pPlane->normal.y * pPlane->dist ;			(*pOut)(3,2) = -2 * pPlane->normal.z * pPlane->dist  ;		(*pOut)(3,3) =  1;
	return pOut;
}

inline float sgn(float a)
{
	if (a > 0.0F) return (1.0F);
	if (a < 0.0F) return (-1.0F);
	return (0.0F);
}

	Matrix44*	MatrixClipProj(	Matrix44* pOutProj,	Matrix44* pViewMat,	Matrix44* pOldProjMat,     const Plane *pOldClipPlane)
{ 
	//<<Oblique Frustum Clipping>> where:http://www.nvidia.com/object/oblique_frustum_clipping.html
	Matrix44 WorldToProjection;
	MatrixMultiply(&WorldToProjection,pViewMat,pOldProjMat);
 
	Plane projClipPlane;
	Matrix44 matInvTranspos = WorldToProjection.Inverse().Transpose();
	PlaneTransform(&projClipPlane,pOldClipPlane,&matInvTranspos);
	// transform clip plane into projection space
 
	if (projClipPlane.dist == 0)  // or less than a really small value
	{
		// plane is perpendicular to the near plane
		*pOutProj = *pOldProjMat;
		return pOutProj;
	}


	if (projClipPlane.dist > 0)
	{
		// flip plane to point away from eye
		Plane clipPlane(-pOldClipPlane->normal.x,-pOldClipPlane->normal.y,-pOldClipPlane->normal.z,-pOldClipPlane->dist);
		// transform clip plane into projection space
		PlaneTransform(&projClipPlane, &clipPlane, &WorldToProjection);

	}

	Matrix44 matClipProj;
	matClipProj.SetIndetity();
	// put projection space clip plane in Z column
	matClipProj(0, 2) = projClipPlane.normal.x;
	matClipProj(1, 2) = projClipPlane.normal.y;
	matClipProj(2, 2) = projClipPlane.normal.z;
	matClipProj(3, 2) = projClipPlane.dist;

	// multiply into projection matrix
	*pOutProj = *pOldProjMat * matClipProj; 

	return pOutProj;
}

 int	LineInterscetPoint(const Vector2& start,const Vector2& end,Vector2& pt)
{
	float res = (end.y-start.y)*pt.x + (start.x-end.x)*pt.y + (end.x*start.y-start.x*end.y);
	if(res< FLT_EPSILON && res> -FLT_EPSILON)
		return 0;
	if(res > FLT_EPSILON)
		return 1;
	else
		return -1;
}

 int	LineInterscetPoint(Vector2& start,Vector2& end,Vector2& pt)
{
	float res = (end.y-start.y)*pt.x + (start.x-end.x)*pt.y + (end.x*start.y-start.x*end.y);
	if(res< FLT_EPSILON && res> -FLT_EPSILON)
		return 0;
	if(res > FLT_EPSILON)
		return 1;
	else
		return -1;
}

 bool LineIntersectLine( Vector2 lineStart0, Vector2 lineEnd0, Vector2 lineStart1, Vector2 lineEnd1,float& t,Vector2& outPt)
{
	float denom = ((lineEnd1.y - lineStart1.y) * (lineEnd0.x - lineStart0.x)) - ((lineEnd1.x - lineStart1.x) * (lineEnd0.y - lineStart0.y));
	float numerator = ((lineEnd1.x - lineStart1.x) * (lineStart0.y - lineStart1.y)) - ((lineEnd1.y - lineStart1.y) * (lineStart0.x - lineStart1.x));

	float numerator2 = ((lineEnd0.x - lineStart0.x) * (lineStart0.y - lineStart1.y)) - ((lineEnd0.y - lineStart0.y) * (lineStart0.x - lineStart1.x));

	if ( denom == 0.0f )
	{
		if ( numerator == 0.0f && numerator2 == 0.0f )
		{
			return false;//COINCIDENT;
		}
		return false;// PARALLEL;
	}
	float ua = numerator / denom;
	float ub = numerator2/ denom;

	if(ua >= 0.0f && ua <= 1.0f && ub >= 0.0f && ub <= 1.0f)
	{
		t = ua;
		outPt = (lineEnd0 - lineStart0)*ua + lineStart0;
		return true;
	}
	else
		return false;
}

float Sign(const Vector2* p1, const Vector2* p2, const Vector2* p3)
{
	return (p1->x - p3->x) * (p2->y - p3->y) - (p2->x - p3->x) * (p1->y - p3->y);
}

 
 bool	IsPointInsideTri(const Vector2 *pa,const Vector2 *pb,
										const Vector2 *pc,const Vector2 *p)
{
	bool b1, b2, b3;

	b1 = Sign(p, pa, pb) < 0.0f;
	b2 = Sign(p, pb, pc) < 0.0f;
	b3 = Sign(p, pc, pa) < 0.0f;

	return ((b1 == b2) && (b2 == b3));
}

 Vector3*  Vec3TransformCoordByQuaternion(Vector3* out, const Vector3* v,const Quaternion* q)
{

	Quaternion p(v->x,v->y,v->z,0);
	Quaternion temp(-q->x, -q->y, -q->z, q->w);

	temp = temp*p*(*q);

	out->x = temp.x;
	out->y = temp.y;
	out->z = temp.z;

	return out;
}