/*
I have no idea how to seperate Geometry and Body.
こちらでメモリを確保するか、それとも参照のみを操るか。
ひとまず、こちらで確保する。


UnionFind
GJK
Voronoi

固定配列ベクターを作ってみる。


*/

/*
M a = f - J^T lambda
J v = - beta C
->
J v(n+1) = - e J v(n) - beta C
m (v(n+1) - v(n)) / dt= f - J^T lambda
->
v(n+1) = v(n) + M^-1 (f dt - J^T lambda dt) 
->
J v(n+1) = J v(n) + J M^-1 (f dt - J^T lambda dt) 
->
- e J v(n) - beta C = J v(n) + J M^-1 (f dt - J^T lambda dt) 
->
0 = J ( (1 + e) v(n) + M^-1 f dt) + beta C - J M^-1 J^T lambda dt 
->
J M^-1 J^T lambda dt = J ( (1 + e) v(n) + M^-1 f dt) + beta C
->
K lambda dt = rhs
K := J M^-1 J^T
rhs(v) := J ( (1 + e) v(n) + M^-1 f dt) + beta C 
->
lambda dt = K^-1 rhs(v)

J vの大きさの和をとり、十分に小さくなるまで繰り返し計算する。

K = Jl1 m1^-1 Jl1^T + Ja1 I1^-1 Ja1^T + Jl2 m2^-1 Jl2^T + Ja2 I2^-1 Ja2^T

K12 =
[Ja1 Jl1] * diag(m^-1,I^-1) * [Ja2 Jl2]^T
=
Ja1 m^-1 Ja2^T + Jl1 I^-1 Jl2^T
行列サイズは、n_c1 x n_c2となる。

0 < beta < 2/dt

J M^-1 J^T
Mi: M^-1
Jxy, Myi
x: joint id
y: body id	

J =
J00   0 
J10 J11
M^-1 = 
M0i   0
    0 M1i   
->
J M^-1 =
J00*M0i       0
J10*M0i J11*M1i
->
J M^-1 J^T =
J00*M0i*J00^T  J00*M0i*J10^T
J10*M0i*J00^T  J10*M0i*J10^T + J11*M1i*J11^T
	
J =
J00   0 J02
	0 J11 J12
J20 J21   0
M^-1 = 
M0i   0   0
    0 M1i   0
    0   0 M2i
->
J M^-1 =
J00*M0i       0 J02*M2i
	    0 J11*M1i J12*M2i
J20*M0i J21*M1i       0
->
J M^-1 J^T =
J00*M0i*J00^T + J02*M2i*J02^T                  J02*M2i*J12^T  J00*M0i*J20^T
	            J12*M2i*J02^T  J11*M1i*J11^T + J12*M2i*J12^T                  J11*M1i*J21^T
J20*M0i*J00^T                  J21*M1i*J11^T                  J20*M0i*J20^T + J21*M1i*J21^T

J =
J00   0 J02
J10 J11   0
	0 J21 J22
J30   0 J32
M^-1 = 
M0i   0   0
    0 M1i   0
    0   0 M2i
->
J M^-1 =
J00*M0i       0 J02*M2i
J10*M0i J11*M1i       0
	    0 J21*M1i J22*M2i
J30*M0i       0 J32*M2i
->
->
J M^-1 J^T =
J00*M0i*J00^T + J02*M2i*J02^T  J00*M0i*J10^T                                  J02*M2i*J22^T  J00*M0i*J30^T + J02*M2i*J32^T
J10*M0i*J00^T                  J10*M0i*J11^T + J11*M1i*J11^T  J11*M1i*J21^T                  J10*M0i*J30^T
	            J22*M2i*J02^T                  J21*M1i*J11^T  J21*M1i*J21^T + J22*M2i*J22^T                  J22*M2i*J32^T 
J30*M0i*J00^T + J32*M2i*J02^T  J30*M0i*J10^T                                  J32*M2i*J22^T  J30*M0i*J30^T + J32*M2i*J32^T


*/


#ifndef NSDEV_CLASSES_MULTIBODY4_HPP
#define NSDEV_CLASSES_MULTIBODY4_HPP

#include <vector>
#include <list>
#include <map>
#include <fstream>
#include <iostream>

#include <nsdev/common3/macros.h>
#include <nsdev/classes/space2.hpp>
#include <nsdev/classes/dynamic_matrix.hpp>
#include <nsdev/functions/matrix.h>
#include <nsdev/classes/function.hpp>
#include <nsdev/classes/time3.hpp>
#include <nsdev/classes/containers/stack.hpp>
#include <nsdev/win32api3/time4.hpp>

#include "containers/list.hpp"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4100) // 引数は関数の本体部で 1 度も参照されません。
#pragma warning(disable:4189) // ローカル変数が初期化されましたが、参照されていません
#endif

namespace nsdev { namespace classes { namespace multibody4 {

using classes::containers::TaskList;
using classes::containers::TaskVector;
using classes::space2::Matrix;
using classes::space2::Vector;
using classes::space2::Quaternion;
using classes::dynamic_matrix::DVector;
using classes::dynamic_matrix::DMatrix;
using classes::function::Function;
using win32api3::time4::PerformanceCounterStopWatch;
using classes::containers::Stack;
using classes::space2::orthogonal_axes;

namespace details
{
	using namespace functions::matrix;

}// namespace details


template <class Real>
class MBUpdateArgs
{
public:
	Real dt;
	Real erp;
};

template <class Real>
class MBWorldUpdateByImpulseArgs : public MBUpdateArgs<Real>
{
public:
	Real iteration_times;
	Real eps_rhs_sum;
};

template <class Real>
class MBWorkspace
{
public:
	Real rhs_sum;
	Real dt_inv;
	Real betadtinv;
};


template <class Real>
class MBTypedef
{
public:
	typedef Matrix<Real,3,3> matrix_type;
	typedef Vector<Real,3> vector_type;
	typedef Quaternion<Real> quaternion_type;

protected:
	typedef matrix_type mat_t;
	typedef vector_type vec_t;
};

struct MBBodyTypeID { enum {
	Body,
	Tire,
};};

template <class Real>
class MBBody : public MBTypedef<Real>
{
	typedef MBBody this_t; 
public:
	virtual ~MBBody(){} 

	MBBody()
	{
		r = (Real)0;
		v = (Real)0;
		w = (Real)0;
		make_identity(R);
		make_identity(q);
	}


	void set_mass(Real m) {this->m = m;}
	void set_inertia(const matrix_type& I_b) {this->I_b = vec_t(I_b._00, I_b._11, I_b._22);}
	void set_inertia(const vector_type& I_b) {this->I_b = I_b;}
	void set_linear_velocity(const vector_type& v) {this->v = v;}
	void set_angular_velocity(const vector_type& w) {this->w = w;}
	void set_position(const vector_type& r) {this->r = r;}

	void set_orientation(const matrix_type& R)
	{
		this->R = R;
		convert(q,R);
	}

	void set_orientation(const quaternion_type& q)
	{
		this->q = q;
		convert(R,q);
	}


	Real m;
	
	vector_type I_b;
	matrix_type I;

	Real m_inv;
	matrix_type I_inv;

	vector_type r; // position
	matrix_type R;
	quaternion_type q;
	vector_type v;
	vector_type w;
	vector_type a;
	vector_type alpha;

	// Solver
	size_t id;


	vector_type f;
	vector_type tau;

	vector_type f_e;
	vector_type tau_e;
};

template <class Real>
class MBTire : public MBBody<Real>
{
public:
	virtual ~MBTire(){} 

	void set_radius(Real r){radius = r;}
	
	Real radius;
};

struct MBJointTypeID { enum {
	Hinge,
	Slider,
	BallPlane,
	Ball,
	Universal,
	Clinder,
	Fixed,
	Angle,
	Distance,
	Orientation,
	Position,
};};


template <class Real>
class MBJoint : public MBTypedef<Real>
{
	typedef MBJoint this_t;
public:
	typedef MBBody<Real> body_type;

	virtual ~MBJoint(){} 
	MBJoint(){}
	MBJoint(const this_t&){}
	this_t& operator=(const this_t&)
	{
		return *this;
	}

	virtual void _precomp() = 0;
	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args) = 0;
	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args) = 0;
	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args) = 0;


	virtual size_t constraint_number() const = 0;
	virtual size_t body_number() const = 0;
	virtual const body_type& body(size_t ibody) const = 0;
	virtual Real* Jl_(size_t ibody) = 0;
	virtual Real* Ja_(size_t ibody) = 0;
	virtual Real* K_() = 0;
	virtual Real* rhs_() = 0;


};

template <class Real, size_t _N_B, size_t _N_C>
class MBJointBase : public MBJoint<Real>
{
public:
	static const size_t N_B = _N_B;
	static const size_t N_C = _N_C;

	body_type* b[N_B];
	Matrix<Real,N_C,3> Jl[N_B];
	Matrix<Real,N_C,3> Ja[N_B];
	Matrix<Real,N_C,N_C> K;
	Vector<Real,N_C> rhs;
	Vector<Real,N_C> C;
	Vector<Real,N_C> lambdadt;

	virtual size_t constraint_number() const {return N_C;}
	virtual size_t body_number() const {return N_B;}
	virtual const body_type& body(size_t ibody) const {return *b[ibody];}
	virtual Real* Jl_(size_t ibody) {return Jl[ibody].data();}
	virtual Real* Ja_(size_t ibody) {return Ja[ibody].data();}
	virtual Real* K_() {return K.data();}
	virtual Real* rhs_() {return rhs.data();}
};

// 未完
// 式から見直し
template <class Real>
class MBDistanceJoint : public MBJointBase<Real,2,1>
{
public:
	MBDistanceJoint() : distance0(0)
	{}

	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	void set_connecting_point(const vector_type& _r0, const vector_type& _r1)
	{
		r0 = _r0;
		r1 = _r1;
	}

	void connecting_point(vector_type& _r0, vector_type& _r1) const 
	{
		_r0 = r0;
		_r1 = r1;
	}

	void set_distance(Real dis){distance0 = dis;}

	virtual void _precomp()
	{
		distance = squared_length(r0 - r1);
		r0 = tmul(b[0]->R, r0 - b[0]->r);
		r1 = tmul(b[1]->R, r1 - b[1]->r);
	}

	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;

		vec_t r01 = b[0]->r - b[1]->r;
		vec_t b0Rr0(b[0]->R * r0);
		vec_t b1Rr1(b[1]->R * r1);

		vec_t tmp(r01 + b0Rr0 - b1Rr1);

		(vec_t&)Jl[0] =  (Real)2 * tmp;
		(vec_t&)Jl[1] = -(const vec_t&)Jl[0];
		(vec_t&)Ja[0] =  cross((const vec_t&)Jl[0], b0Rr0);
		(vec_t&)Ja[1] =  cross((const vec_t&)Jl[1], b1Rr1);

		C[0] = squared_length(tmp) - distance - distance0;

		rhs = (Real)0;
		rhs[0] += dot((const vec_t&)Jl[0], b[0]->v + b[0]->m_inv*b[0]->f*dt);
		rhs[0] += dot((const vec_t&)Ja[0], b[0]->w + b[0]->I_inv*b[0]->tau*dt);
		rhs[0] += dot((const vec_t&)Jl[1], b[1]->v + b[1]->m_inv*b[1]->f*dt);
		rhs[0] += dot((const vec_t&)Ja[1], b[1]->w + b[1]->I_inv*b[1]->tau*dt);

		rhs += betadtinv * C;

		K[0] = (Real)0;
		K[0] += dot((const vec_t&)Jl[0], b[0]->m_inv * (const vec_t&)Jl[0]);
		K[0] += dot((const vec_t&)Jl[1], b[1]->m_inv * (const vec_t&)Jl[1]);
		K[0] += dot((const vec_t&)Ja[0], b[0]->I_inv * (const vec_t&)Ja[0]);
		K[0] += dot((const vec_t&)Ja[1], b[1]->I_inv * (const vec_t&)Ja[1]);
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};
	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};

public:
	vector_type r0, r1;
	Real distance;
	Real distance0;
};



template <class Real>
class MBAngleJoint : public MBJointBase<Real,2,1>
{
public:
	MBAngleJoint() : theta0(0), dtheta(0)
	{
		(vec_t&)Jl[0] = vec_t((Real)0);
		(vec_t&)Jl[1] = vec_t((Real)0);
	}

	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	void set_axis(const vector_type& axis)
	{
		axis0 = axis;
		orth0 = orthogonal_axis(axis);
	}

	const vector_type& axis() const {return axis0;}

	void set_angle(Real ang){theta0 = ang;}

	Real angle() const{return theta0;}

	virtual void _precomp()
	{
		axis1 = tmul(b[1]->R, axis0);
		axis0 = tmul(b[0]->R, axis0);

		orth1 = tmul(b[1]->R, orth0);
		orth0 = tmul(b[0]->R, orth0);
	}

	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;

		vec_t b0Raxis0 = b[0]->R * axis0;
		vec_t b1Raxis1 = b[1]->R * axis1;
		vec_t b0Rorth0 = b[0]->R * orth0;
		vec_t b1Rorth1 = b[1]->R * orth1;


		(vec_t&)Ja[0] =  b0Raxis0;
		(vec_t&)Ja[1] = -b1Raxis1;

		dtheta = ::angle(b1Rorth1, b0Rorth0, b0Raxis0);

		C[0] = dtheta - theta0;

		rhs = (Real)0;
		//rhs[0] += dot((const vec_t&)Jl[0], b[0]->v + b[0]->m_inv*b[0]->f*dt);
		rhs[0] += dot((const vec_t&)Ja[0], b[0]->w + b[0]->I_inv*b[0]->tau*dt);
		//rhs[0] += dot((const vec_t&)Jl[1], b[1]->v + b[1]->m_inv*b[1]->f*dt);
		rhs[0] += dot((const vec_t&)Ja[1], b[1]->w + b[1]->I_inv*b[1]->tau*dt);

		rhs += betadtinv * C;

		K[0] = (Real)0;
		//K[0] += dot((const vec_t&)Jl[0], b[0]->m_inv * (const vec_t&)Jl[0]);
		//K[0] += dot((const vec_t&)Jl[1], b[1]->m_inv * (const vec_t&)Jl[1]);
		K[0] += dot((const vec_t&)Ja[0], b[0]->I_inv * (const vec_t&)Ja[0]);
		K[0] += dot((const vec_t&)Ja[1], b[1]->I_inv * (const vec_t&)Ja[1]);
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};
	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};

public:
	vector_type axis0, orth0;
	vector_type axis1, orth1;
	Real theta0;
	Real dtheta;
};


template <class Real>
class MBBallPlaneJoint : public MBJointBase<Real,2,1>
{
public:
	MBBallPlaneJoint()
	{}

	void link_for_ball_constraint(body_type* b0)
	{
		b[0] = b0;
	}

	void link_for_plane_constraint(body_type* b1)
	{
		b[1] = b1;
	}

	void set_plane_direction(const vector_type& n)
	{
		this->n = n;
	}

	void set_point_on_plane(const vector_type& s)
	{
		this->s = s;
	}

	const vector_type& plane_direction() const {return n;}
	const vector_type& point_on_plane() const {return s;}

	virtual void _precomp()
	{
		d = dot(n, b[0]->r + s - b[1]->r);
		//offset = tmul(b[1]->R, b[0]->r + s - b[1]->r);

		s = tmul(b[0]->R, s);
		n = tmul(b[1]->R, n);
	}

	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;

		vec_t rij_(b[0]->r - b[1]->r);
		vec_t b0Rs(b[0]->R * s);
		vec_t b1Rn(b[1]->R * n);
		//vec_t offset_(b[1]->R * offset);

		(vec_t&)Jl[0] =  b1Rn;
		(vec_t&)Jl[1] = -b1Rn;
		(vec_t&)Ja[0] = cross(-b1Rn, b0Rs);
		(vec_t&)Ja[1] = cross( b1Rn, rij_ + b0Rs);

		C = dot(b1Rn, rij_ + b0Rs) - d;
		//C = dot(b1Rn, rij_ + b0Rs - offset_);

		rhs = (Real)0;
		rhs += dot((const vec_t&)Jl[0], b[0]->v + b[0]->m_inv*b[0]->f*dt);
		rhs += dot((const vec_t&)Ja[0], b[0]->w + b[0]->I_inv*b[0]->tau*dt);
		rhs += dot((const vec_t&)Jl[1], b[1]->v + b[1]->m_inv*b[1]->f*dt);
		rhs += dot((const vec_t&)Ja[1], b[1]->w + b[1]->I_inv*b[1]->tau*dt);

		rhs += betadtinv * C;

		K = (Real)0;
		K[0] += dot((const vec_t&)Jl[0], b[0]->m_inv * (const vec_t&)Jl[0]);
		K[0] += dot((const vec_t&)Jl[1], b[1]->m_inv * (const vec_t&)Jl[1]);
		K[0] += dot((const vec_t&)Ja[0], b[0]->I_inv * (const vec_t&)Ja[0]);
		K[0] += dot((const vec_t&)Ja[1], b[1]->I_inv * (const vec_t&)Ja[1]);
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};
	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};

public:
	vector_type s;
	vector_type n;
	//vector_type offset;
	Real d;
};


template <class Real>
class MBSkyBallJoint : public MBJointBase<Real,1,3>
{
public:
	MBSkyBallJoint()
	{
		make_identity(Jl[0]);
	}

	void set_rotation_point(const vector_type& rp)
	{
		r0 = rp;
	}

	void set_sky_point(const vector_type& sp)
	{
		s = sp;
	}

	void link(body_type* b0)
	{
		b[0] = b0;
	}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r0 - b[0]->r);
	}


	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;
		vec_t b0Rr0 = b[0]->R * r0;
		
		Ja[0] =  - skew(b[0]->R * r0);
		C = b[0]->r + b[0]->R * r0 - s; 

		rhs = betadtinv * C;
		//rhs += Jl[0] * (b[0]->v + b[0]->m_inv*b[0]->f*dt);
		rhs += b[0]->v + b[0]->m_inv*b[0]->f*dt;
		//rhs += Ja[0] * (b[0]->w + b[0]->I_inv*b[0]->tau*dt);
		rhs += cross(-b0Rr0, (b[0]->w + b[0]->I_inv*b[0]->tau*dt));	
		
		K = (Real)0;
		//K += Jl[0] * b[0]->m_inv * transposed(Jl[0]);
		K += diag(vec_t(b[0]->m_inv));
		K += Ja[0] * b[0]->I_inv * transposed(Ja[0]);
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		vec_t b0Rr0 = b[0]->R * r0;
		Ja[0] =  skew(-b0Rr0);
		C = b[0]->r + b0Rr0 - s; 
	}

	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		Real betadtinv = ws.betadtinv;
		vec_t b0Rr0 = b[0]->R * r0;

		rhs = (Real)0;
		//rhs += Jl[0] * (b[0]->v);
		rhs += b[0]->v;
		//rhs += Ja[0] * (b[0]->w);
		rhs += cross(-b0Rr0, b[0]->w);
		
		ws.rhs_sum += squared_length(rhs);
		
		rhs += betadtinv * C;
		

		K = (Real)0;
		K += diag(vec_t(b[0]->m_inv));
		K += Ja[0] * b[0]->I_inv * transposed(Ja[0]);
		lambdadt = inverse(K) * rhs;
		
		b[0]->v -= b[0]->m_inv * transposed(Jl[0]) * lambdadt;
		b[0]->w -= b[0]->I_inv * transposed(Ja[0]) * lambdadt;
	}

private:
	vector_type r0;
	vector_type s ;
};

template <class Real>
class MBBallJoint : public MBJointBase<Real,2,3>
{
public:
	MBBallJoint()
	{
		make_identity(Jl[0]);
		make_identity(Jl[1]);
		Jl[1] = -Jl[1];
	}

	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	void set_rotation_point(const vector_type& rp)
	{
		r1 = rp;
	}
	
	const vector_type& rotation_point() const {return r1;}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r1 - b[0]->r);
		r1 = tmul(b[1]->R, r1 - b[1]->r);
	}

	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		//static int counter = 0;
		
		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;
		vec_t b0Rr0 = b[0]->R * r0;
		vec_t b1Rr1 = b[1]->R * r1;


		Ja[0] =  skew(-b0Rr0);
		Ja[1] =  skew( b1Rr1);

		C = b[0]->r + b0Rr0 - (b[1]->r + b1Rr1); 

		rhs = (Real)0;
		//rhs += Jl[0] * (b[0]->v + b[0]->m_inv*b[0]->f*dt);		
		rhs += b[0]->v + b[0]->m_inv*b[0]->f*dt;
		//rhs += Ja[0] * (b[0]->w + b[0]->I_inv*b[0]->tau*dt);	
		rhs -= cross(b0Rr0, b[0]->w + b[0]->I_inv*b[0]->tau*dt);	
		//rhs += Jl[1] * (b[1]->v + b[1]->m_inv*b[1]->f*dt);
		rhs -= b[1]->v + b[1]->m_inv*b[1]->f*dt;
		//rhs += Ja[1] * (b[1]->w + b[1]->I_inv*b[1]->tau*dt);
		rhs += cross(b1Rr1, b[1]->w + b[1]->I_inv*b[1]->tau*dt);	

		rhs += betadtinv * C;
		
#if 0
		if(counter == 0)
		{
			using namespace std;
			cout << "C: " << C << endl;
			cout << "b0Rr0: " << b0Rr0 << endl;
			cout << "b1Rr1: " << b1Rr1 << endl;
			cout << "b[0]->I_inv: " << b[0]->I_inv << endl;
			cout << "b[1]->I_inv: " << b[1]->I_inv << endl;
			cout << "rhs: " << rhs << endl;
			cout << "betadtinv: " << betadtinv << endl;
			cout << "b[0]->v: " << b[0]->v << endl;
			cout << "b[1]->v: " << b[1]->v << endl;
			cout << "b[0]->w: " << b[0]->w << endl;
			cout << "b[1]->w: " << b[1]->w << endl;
			cout << "b[0]->f: " << b[0]->f << endl;
			cout << "b[1]->f: " << b[1]->f << endl;
			cout << "b[0]->tau: " << b[0]->tau << endl;
			cout << "b[1]->tau: " << b[1]->tau << endl;
		}
#endif

		K = (Real)0;
		//K += Jl[0] * b[0]->m_inv * transposed(Jl[0]);
		//K += Jl[1] * b[1]->m_inv * transposed(Jl[1]);
		K += diag(vec_t(b[0]->m_inv));
		K += diag(vec_t(b[1]->m_inv));
		K += Ja[0] * b[0]->I_inv * transposed(Ja[0]);
		K += Ja[1] * b[1]->I_inv * transposed(Ja[1]);

		//++counter;
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		vec_t b0Rr0 = b[0]->R * r0;
		vec_t b1Rr1 = b[1]->R * r1;
		Ja[0] = skew(-b0Rr0);
		Ja[1] = skew( b1Rr1);
		C = b[0]->r + b0Rr0 - (b[1]->r + b1Rr1); 
	}

	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;
		vec_t b0Rr0 = b[0]->R * r0;
		vec_t b1Rr1 = b[1]->R * r1;

		rhs = (Real)0;
		//rhs += Jl[0] * (b[0]->v);
		rhs += b[0]->v;
		//rhs += Jl[1] * (b[1]->v);
		rhs -= b[1]->v;
		//rhs += Ja[0] * (b[0]->w);
		rhs += cross(-b0Rr0, b[0]->w);
		//rhs += Ja[1] * (b[1]->w);
		rhs += cross( b1Rr1, b[1]->w);

		ws.rhs_sum += squared_length(rhs);

		rhs += betadtinv * C;


		K = (Real)0;
		K += diag(vec_t(b[0]->m_inv));
		K += diag(vec_t(b[1]->m_inv));
		K += Ja[0] * b[0]->I_inv * transposed(Ja[0]);
		K += Ja[1] * b[1]->I_inv * transposed(Ja[1]);
		lambdadt = inverse(K) * rhs;

		b[0]->v -= b[0]->m_inv * transposed(Jl[0]) * lambdadt;
		b[0]->w -= b[0]->I_inv * transposed(Ja[0]) * lambdadt;
		b[1]->v -= b[1]->m_inv * transposed(Jl[1]) * lambdadt;
		b[1]->w -= b[1]->I_inv * transposed(Ja[1]) * lambdadt;
	}

private:
	vector_type r0, r1;
};



template <class Real>
class MBUniversalJoint : public MBJointBase<Real,2,4>
{
public:
	MBUniversalJoint()
	{
		make_identity((mat_t&)Jl[0]);
		make_identity((mat_t&)Jl[1]);
		(mat_t&)Jl[1] = -(const mat_t&)Jl[1];

		(vec_t&)Jl[0](3) = vec_t((Real)0);
		(vec_t&)Jl[1](3) = vec_t((Real)0);
	}

	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	void set_rotation_point(const vector_type& rp)
	{
		r1 = rp;
	}

	void set_axis(const vector_type& axis1, const vector_type& axis2)
	{
		axis1i = axis1;
		axis2j = axis2;
		axdot = dot(axis2j,axis1i);
	}

	void axis(vector_type& ax1, vector_type& ax2) const
	{
		ax1 = axis1i;
		ax2 = axis2j;
	}
	
	const vector_type& rotation_point() const {return r1;}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r1 - b[0]->r);
		r1 = tmul(b[1]->R, r1 - b[1]->r);
		axis1i = tmul(b[0]->R, axis1i);
		axis2j = tmul(b[1]->R, axis2j);
	}

	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		typedef Vector<Real,3> vec3_t;
		typedef Matrix<Real,1,3> mat13_t;
		typedef Matrix<Real,3,3> mat33_t;

		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;
		vec_t b0Rr0 = b[0]->R * r0;
		vec_t b1Rr1 = b[1]->R * r1;

		vec_t axis1i_ = b[0]->R * axis1i;
		vec_t axis2j_ = b[1]->R * axis2j;

		(mat_t&)Ja[0] =  skew(-b0Rr0);
		(mat_t&)Ja[1] =  skew( b1Rr1);

		(vec_t&)C = b[0]->r + b0Rr0 - (b[1]->r + b1Rr1); 

		rhs = (Real)0;
		//(vec_t&)rhs += Jl[0] * (b[0]->v + b[0]->m_inv*b[0]->f*dt);		
		(vec_t&)rhs += b[0]->v + b[0]->m_inv*b[0]->f*dt;
		//(vec_t&)rhs += Ja[0] * (b[0]->w + b[0]->I_inv*b[0]->tau*dt);	
		(vec_t&)rhs += cross(-b0Rr0,(b[0]->w + b[0]->I_inv*b[0]->tau*dt));	
		//(vec_t&)rhs += Jl[1] * (b[1]->v + b[1]->m_inv*b[1]->f*dt);
		(vec_t&)rhs -= b[1]->v + b[1]->m_inv*b[1]->f*dt;
		//(vec_t&)rhs += Ja[1] * (b[1]->w + b[1]->I_inv*b[1]->tau*dt);
		(vec_t&)rhs += cross(b1Rr1,(b[1]->w + b[1]->I_inv*b[1]->tau*dt));	

		vec_t Ja03 = cross(axis1i_, axis2j_);
		(vec_t&)Ja[0](3) =   Ja03; 
		(vec_t&)Ja[1](3) = - Ja03; 

		C[3] = dot(axis2j_, axis1i_) - axdot;

		//rhs[3] += dot((const vec_t&)Jl[0](3), b[0]->v + b[0]->m_inv*b[0]->f*dt); 
		//rhs[3] += dot((const vec_t&)Jl[1](3), b[1]->v + b[1]->m_inv*b[1]->f*dt); 
		rhs[3] += dot((const vec_t&)Ja[0](3), b[0]->w + b[0]->I_inv*b[0]->tau*dt); 
		rhs[3] += dot((const vec_t&)Ja[1](3), b[1]->w + b[1]->I_inv*b[1]->tau*dt); 

		rhs += betadtinv * C;


		K = (Real)0;

		//add(K, 0, 0, (const mat33_t&)Jl[0](0) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		//add(K, 0, 0, (const mat33_t&)Jl[1](0) * b[1]->m_inv * transposed((const mat33_t&)Jl[1](0)));
		add(K, 0, 0, mat33_t(b[0]->m_inv, b[0]->m_inv, b[0]->m_inv));
		add(K, 0, 0, mat33_t(b[1]->m_inv, b[1]->m_inv, b[1]->m_inv));
		add(K, 0, 0, (const mat33_t&)Ja[0](0) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));
		add(K, 0, 0, (const mat33_t&)Ja[1](0) * b[1]->I_inv * transposed((const mat33_t&)Ja[1](0)));

		//add(K, 3, 0, (const mat13_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		//add(K, 3, 0, (const mat13_t&)Jl[1](3) * b[1]->m_inv * transposed((const mat33_t&)Jl[1](0)));
		add(K, 3, 0, (const mat13_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));
		add(K, 3, 0, (const mat13_t&)Ja[1](3) * b[1]->I_inv * transposed((const mat33_t&)Ja[1](0)));

		//add(K, 3, 3, (const mat13_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat13_t&)Jl[0](3)));
		//add(K, 3, 3, (const mat13_t&)Jl[1](3) * b[1]->m_inv * transposed((const mat13_t&)Jl[1](3)));
		add(K, 3, 3, (const mat13_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat13_t&)Ja[0](3)));
		add(K, 3, 3, (const mat13_t&)Ja[1](3) * b[1]->I_inv * transposed((const mat13_t&)Ja[1](3)));

	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){}

	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){}

private:
	vector_type r0, r1;
	vector_type axis1i, axis2j;
	Real axdot;
	vector_type orth1, orth2;
};


template <class Real>
class MBSkySliderJoint : public MBJointBase<Real,1,5>
{
public:
	MBSkySliderJoint()
	{
		(mat_t&)Jl[0](0) = mat_t((Real)0);

		(vec_t&)Ja[0](3) = vec_t((Real)0);
		(vec_t&)Ja[0](4) = vec_t((Real)0);
	}

	void link(body_type* b0)
	{
		b[0] = b0;
	}

	void set_axis(const vector_type& axis_)
	{
		orth0i = axis_;
		orthogonal_axes(orth1j, orth2j, orth0i);
	}

	virtual void _precomp()
	{
		rij = tmul(b[0]->R, b[0]->r);
		orth0i = tmul(b[0]->R, orth0i);
		orth1i = tmul(b[0]->R, orth1j);
	}
	
	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		typedef Vector<Real,3> vec3_t;
		typedef Matrix<Real,2,3> mat23_t;
		typedef Matrix<Real,3,3> mat33_t;

		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;

		vec_t rij_(b[0]->r);
		vec_t b0Rrij(b[0]->R * rij);
		vec_t orth0i_(b[0]->R * orth0i);
		vec_t orth1i_(b[0]->R * orth1i);
		vec_t orth1j_(orth1j);
		vec_t orth2j_(orth2j);

		(vec_t&)Ja[0](0) = cross(orth1i_,orth2j_);
		(vec_t&)Ja[0](1) = cross(orth0i_,orth1j_);
		(vec_t&)Ja[0](2) = cross(orth0i_,orth2j_);

		C[0] = dot(orth2j_,orth1i_);
		C[1] = dot(orth1j_,orth0i_);
		C[2] = dot(orth2j_,orth0i_);
		
		rhs = (Real)0;
		//(vec_t&)rhs[0] += (mat_t&)Jl[0](0) * (b[0]->v + b[0]->m_inv*b[0]->f*dt);		
		(vec_t&)rhs[0] += (mat_t&)Ja[0](0) * (b[0]->w + b[0]->I_inv*b[0]->tau*dt);	


		(vec_t&)Jl[0](3) =  orth1j_;
		(vec_t&)Jl[0](4) =  orth2j_;

		vec_t rij_b0Rrij(rij_ - b0Rrij);
		C[3] = dot(orth1j_, rij_b0Rrij);
		C[4] = dot(orth2j_, rij_b0Rrij);

		rhs[3] += dot((const vec_t&)Jl[0](3), b[0]->v + b[0]->m_inv*b[0]->f*dt); 
		//rhs[3] += dot((const vec_t&)Ja[0](3), b[0]->w + b[0]->I_inv*b[0]->tau*dt); 
		rhs[4] += dot((const vec_t&)Jl[0](4), b[0]->v + b[0]->m_inv*b[0]->f*dt); 
		//rhs[4] += dot((const vec_t&)Ja[0](4), b[0]->w + b[0]->I_inv*b[0]->tau*dt); 

		rhs += betadtinv * C;

		K = (Real)0;

		//add(K, 0, 0, (const mat33_t&)Jl[0](0) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		add(K, 0, 0, (const mat33_t&)Ja[0](0) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));

		//add(K, 3, 0, (const mat23_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		//add(K, 3, 0, (const mat23_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));

		add(K, 3, 3, (const mat23_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat23_t&)Jl[0](3)));
		//add(K, 3, 3, (const mat23_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat23_t&)Ja[0](3)));
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};
	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};


	vector_type rij;
	vector_type orth0i, orth1i, orth1j, orth2j;
};

template <class Real>
class MBSliderJoint : public MBJointBase<Real,2,5>
{
public:
	MBSliderJoint()
	{
		(mat_t&)Jl[0](0) = mat_t((Real)0);
		(mat_t&)Jl[1](0) = mat_t((Real)0);

		(vec_t&)Ja[0](3) = vec_t((Real)0);
		(vec_t&)Ja[0](4) = vec_t((Real)0);
	}

	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	void set_axis(const vector_type& axis_)
	{
		orth0i = axis_;
		orthogonal_axes(orth1j, orth2j, orth0i);
	}
	
	const vector_type& axis() const {return orth0i;}
	const vector_type& c_axis() const {return b[0]->R * orth0i;}


	virtual void _precomp()
	{
		rij = tmul(b[0]->R, b[0]->r - b[1]->r);
		orth0i = tmul(b[0]->R, orth0i);
		orth1i = tmul(b[0]->R, orth1j);
		orth1j = tmul(b[1]->R, orth1j);
		orth2j = tmul(b[1]->R, orth2j);
	}
	
	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		typedef Vector<Real,3> vec3_t;
		typedef Matrix<Real,2,3> mat23_t;
		typedef Matrix<Real,3,3> mat33_t;

		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;

		vec_t rij_(b[0]->r - b[1]->r);
		vec_t b0Rrij(b[0]->R * rij);
		vec_t orth0i_(b[0]->R * orth0i);
		vec_t orth1i_(b[0]->R * orth1i);
		vec_t orth1j_(b[1]->R * orth1j);
		vec_t orth2j_(b[1]->R * orth2j);

		(vec_t&)Ja[0](0) = cross(orth1i_,orth2j_);
		(vec_t&)Ja[0](1) = cross(orth0i_,orth1j_);
		(vec_t&)Ja[0](2) = cross(orth0i_,orth2j_);
		(vec_t&)Ja[1](0) = -(const vec_t&)Ja[0](0);
		(vec_t&)Ja[1](1) = -(const vec_t&)Ja[0](1);
		(vec_t&)Ja[1](2) = -(const vec_t&)Ja[0](2);

		C[0] = dot(orth2j_,orth1i_);
		C[1] = dot(orth1j_,orth0i_);
		C[2] = dot(orth2j_,orth0i_);
		
		rhs = (Real)0;
		//(vec_t&)rhs[0] += (mat_t&)Jl[0](0) * (b[0]->v + b[0]->m_inv*b[0]->f*dt);		
		(vec_t&)rhs[0] += (mat_t&)Ja[0](0) * (b[0]->w + b[0]->I_inv*b[0]->tau*dt);	
		//(vec_t&)rhs[0] += (mat_t&)Jl[1](0) * (b[1]->v + b[1]->m_inv*b[1]->f*dt);
		(vec_t&)rhs[0] += (mat_t&)Ja[1](0) * (b[1]->w + b[1]->I_inv*b[1]->tau*dt);


		(vec_t&)Jl[0](3) =  orth1j_;
		(vec_t&)Jl[1](3) = -orth1j_;
		(vec_t&)Ja[1](3) = cross(orth1j_, rij_);
		(vec_t&)Jl[0](4) =  orth2j_;
		(vec_t&)Jl[1](4) = -orth2j_;
		(vec_t&)Ja[1](4) = cross(orth2j_, rij_);

		vec_t rij_b0Rrij(rij_ - b0Rrij);
		C[3] = dot(orth1j_, rij_b0Rrij);
		C[4] = dot(orth2j_, rij_b0Rrij);

		rhs[3] += dot((const vec_t&)Jl[0](3), b[0]->v + b[0]->m_inv*b[0]->f*dt); 
		rhs[3] += dot((const vec_t&)Jl[1](3), b[1]->v + b[1]->m_inv*b[1]->f*dt); 
		//rhs[3] += dot((const vec_t&)Ja[0](3), b[0]->w + b[0]->I_inv*b[0]->tau*dt); 
		rhs[3] += dot((const vec_t&)Ja[1](3), b[1]->w + b[1]->I_inv*b[1]->tau*dt); 
		rhs[4] += dot((const vec_t&)Jl[0](4), b[0]->v + b[0]->m_inv*b[0]->f*dt); 
		rhs[4] += dot((const vec_t&)Jl[1](4), b[1]->v + b[1]->m_inv*b[1]->f*dt); 
		//rhs[4] += dot((const vec_t&)Ja[0](4), b[0]->w + b[0]->I_inv*b[0]->tau*dt); 
		rhs[4] += dot((const vec_t&)Ja[1](4), b[1]->w + b[1]->I_inv*b[1]->tau*dt); 
		
		rhs += betadtinv * C;

		K = (Real)0;

		//add(K, 0, 0, (const mat33_t&)Jl[0](0) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		//add(K, 0, 0, (const mat33_t&)Jl[1](0) * b[1]->m_inv * transposed((const mat33_t&)Jl[1](0)));
		add(K, 0, 0, (const mat33_t&)Ja[0](0) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));
		add(K, 0, 0, (const mat33_t&)Ja[1](0) * b[1]->I_inv * transposed((const mat33_t&)Ja[1](0)));

		//add(K, 3, 0, (const mat23_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		//add(K, 3, 0, (const mat23_t&)Jl[1](3) * b[1]->m_inv * transposed((const mat33_t&)Jl[1](0)));
		//add(K, 3, 0, (const mat23_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));
		add(K, 3, 0, (const mat23_t&)Ja[1](3) * b[1]->I_inv * transposed((const mat33_t&)Ja[1](0)));

		add(K, 3, 3, (const mat23_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat23_t&)Jl[0](3)));
		add(K, 3, 3, (const mat23_t&)Jl[1](3) * b[1]->m_inv * transposed((const mat23_t&)Jl[1](3)));
		//add(K, 3, 3, (const mat23_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat23_t&)Ja[0](3)));
		add(K, 3, 3, (const mat23_t&)Ja[1](3) * b[1]->I_inv * transposed((const mat23_t&)Ja[1](3)));
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};
	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};


	vector_type rij;
	vector_type orth0i, orth1i, orth1j, orth2j;
};


template <class Real>
class MBSkyHingeJoint : public MBJointBase<Real,1,5>
{
public:
	MBSkyHingeJoint()
	{
		make_identity((mat_t&)Jl[0](0));
		(vec_t&)Jl[0](3) = vec_t((Real)0);
		(vec_t&)Jl[0](4) = vec_t((Real)0);
	}

	void link(body_type* b0)
	{
		b[0] = b0;
	}

	void set_sky_point(const vector_type& sp)
	{
		s = sp;
	}

	void set_rotation_point(const vector_type& rp)
	{
		r0 = rp;
	}

	void set_axis(const vector_type& axis)
	{
		parai = axis;
		orthogonal_axes(orth1j, orth2j, parai);
	}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r0 - b[0]->r);

		parai  = tmul(b[0]->R, parai );
	}

	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		typedef Vector<Real,3> vec3_t;
		typedef Matrix<Real,2,3> mat23_t;
		typedef Matrix<Real,3,3> mat33_t;

		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;

		vec_t b0Rr0 = b[0]->R * r0;
		vec_t parai_(b[0]->R * parai);
		vec_t orth1j_(orth1j);
		vec_t orth2j_(orth2j);

		(mat_t&)Ja[0](0) = skew(-b0Rr0);

		(vec_t&)C[0] = b[0]->r + b0Rr0 - s; 
		
		rhs = (Real)0;
		//(vec_t&)rhs[0] += (mat_t&)Jl[0](0) * (b[0]->v + b[0]->m_inv*b[0]->f*dt);		
		(vec_t&)rhs[0] += b[0]->v + b[0]->m_inv*b[0]->f*dt;
		//(vec_t&)rhs[0] += (mat_t&)Ja[0](0) * (b[0]->w + b[0]->I_inv*b[0]->tau*dt);	
		(vec_t&)rhs[0] += cross(-b0Rr0,(b[0]->w + b[0]->I_inv*b[0]->tau*dt));	

		vec_t Ja03 = cross(parai, orth1j_);
		vec_t Ja04 = cross(parai, orth2j_);

		(vec_t&)Ja[0](3) =  Ja03;
		(vec_t&)Ja[0](4) =  Ja04;

		C[3] = dot(orth1j_, parai_); 
		C[4] = dot(orth2j_, parai_); 

		rhs[3] += dot((const vec_t&)Ja[0](3), b[0]->w + b[0]->I_inv*b[0]->tau*dt); 
		rhs[4] += dot((const vec_t&)Ja[0](4), b[0]->w + b[0]->I_inv*b[0]->tau*dt); 

		rhs += betadtinv * C;


		K = (Real)0;

		//add(K, 0, 0, (const mat33_t&)Jl[0](0) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		add(K, 0, 0, mat33_t(b[0]->m_inv, b[0]->m_inv, b[0]->m_inv));
		add(K, 0, 0, (const mat33_t&)Ja[0](0) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));

		//add(K, 3, 0, (const mat23_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		add(K, 3, 0, (const mat23_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));

		//add(K, 3, 3, (const mat23_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat23_t&)Jl[0](3)));
		add(K, 3, 3, (const mat23_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat23_t&)Ja[0](3)));
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};
	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};

	vector_type s;
	vector_type r0;
	
	vector_type parai;
	vector_type orth1j, orth2j;
};

template <class Real>
class MBHingeJoint : public MBJointBase<Real,2,5>
{
public:
	MBHingeJoint()
	{
		make_identity((mat_t&)Jl[0](0));
		make_identity((mat_t&)Jl[1](0));
		(mat_t&)Jl[1](0) = -(const mat_t&)Jl[1](0);
		(vec_t&)Jl[0](3) = vec_t((Real)0);
		(vec_t&)Jl[1](3) = vec_t((Real)0);
		(vec_t&)Jl[0](4) = vec_t((Real)0);
		(vec_t&)Jl[1](4) = vec_t((Real)0);
	}

	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	void set_rotation_point(const vector_type& rp)
	{
		r1 = rp;
	}

	void set_axis(const vector_type& axis)
	{
		parai = axis;
		orthogonal_axes(orth1j, orth2j, parai);
	}

	const vector_type& axis() const {return parai;}
	const vector_type& rotation_point() const {return r1;}

	const vector_type& c_axis() const {return b[0]->R * parai;}
	const vector_type& c_rotation_point() const {return b[1]->R * r1;}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r1 - b[0]->r);
		r1 = tmul(b[1]->R, r1 - b[1]->r);

		parai  = tmul(b[0]->R, parai );
		orth1j = tmul(b[1]->R, orth1j);
		orth2j = tmul(b[1]->R, orth2j);
	}

	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		typedef Vector<Real,3> vec3_t;
		typedef Matrix<Real,2,3> mat23_t;
		typedef Matrix<Real,3,3> mat33_t;

		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;

		vec_t b0Rr0 = b[0]->R * r0;
		vec_t b1Rr1 = b[1]->R * r1;
		vec_t parai_(b[0]->R * parai);
		vec_t orth1j_(b[1]->R * orth1j);
		vec_t orth2j_(b[1]->R * orth2j);

		(mat_t&)Ja[0](0) = skew(-b0Rr0);
		(mat_t&)Ja[1](0) = skew( b1Rr1);

		(vec_t&)C[0] = b[0]->r + b0Rr0 - (b[1]->r + b1Rr1);
		
		rhs = (Real)0;
		//(vec_t&)rhs[0] += (mat_t&)Jl[0](0) * (b[0]->v + b[0]->m_inv*b[0]->f*dt);
		(vec_t&)rhs[0] += b[0]->v + b[0]->m_inv*b[0]->f*dt;
		//(vec_t&)rhs[0] += (mat_t&)Ja[0](0) * (b[0]->w + b[0]->I_inv*b[0]->tau*dt);
		(vec_t&)rhs[0] += cross(-b0Rr0,(b[0]->w + b[0]->I_inv*b[0]->tau*dt));
		//(vec_t&)rhs[0] += (mat_t&)Jl[1](0) * (b[1]->v + b[1]->m_inv*b[1]->f*dt);
		(vec_t&)rhs[0] -= b[1]->v + b[1]->m_inv*b[1]->f*dt;
		//(vec_t&)rhs[0] += (mat_t&)Ja[1](0) * (b[1]->w + b[1]->I_inv*b[1]->tau*dt);
		(vec_t&)rhs[0] += cross(b1Rr1,(b[1]->w + b[1]->I_inv*b[1]->tau*dt));

		vec_t Ja03 = cross(parai_, orth1j_);
		vec_t Ja04 = cross(parai_, orth2j_);

		(vec_t&)Ja[0](3) =  Ja03;
		(vec_t&)Ja[1](3) = -Ja03;
		(vec_t&)Ja[0](4) =  Ja04;
		(vec_t&)Ja[1](4) = -Ja04;

		C[3] = dot(orth1j_, parai_); 
		C[4] = dot(orth2j_, parai_); 

		//rhs[3] += dot((const vec_t&)Jl[0](3), b[0]->v + b[0]->m_inv*b[0]->f*dt); 
		//rhs[3] += dot((const vec_t&)Jl[1](3), b[1]->v + b[1]->m_inv*b[1]->f*dt); 
		rhs[3] += dot((const vec_t&)Ja[0](3), b[0]->w + b[0]->I_inv*b[0]->tau*dt); 
		rhs[3] += dot((const vec_t&)Ja[1](3), b[1]->w + b[1]->I_inv*b[1]->tau*dt);
		//rhs[4] += dot((const vec_t&)Jl[0](4), b[0]->v + b[0]->m_inv*b[0]->f*dt); 
		//rhs[4] += dot((const vec_t&)Jl[1](4), b[1]->v + b[1]->m_inv*b[1]->f*dt); 
		rhs[4] += dot((const vec_t&)Ja[0](4), b[0]->w + b[0]->I_inv*b[0]->tau*dt); 
		rhs[4] += dot((const vec_t&)Ja[1](4), b[1]->w + b[1]->I_inv*b[1]->tau*dt); 

		rhs += betadtinv * C;


		K = (Real)0;

		//add(K, 0, 0, (const mat33_t&)Jl[0](0) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		//add(K, 0, 0, (const mat33_t&)Jl[1](0) * b[1]->m_inv * transposed((const mat33_t&)Jl[1](0)));
		add(K, 0, 0, mat33_t(b[0]->m_inv, b[0]->m_inv, b[0]->m_inv));
		add(K, 0, 0, mat33_t(b[1]->m_inv, b[1]->m_inv, b[1]->m_inv));
		add(K, 0, 0, (const mat33_t&)Ja[0](0) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));
		add(K, 0, 0, (const mat33_t&)Ja[1](0) * b[1]->I_inv * transposed((const mat33_t&)Ja[1](0)));

		//add(K, 3, 0, (const mat23_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		//add(K, 3, 0, (const mat23_t&)Jl[1](3) * b[1]->m_inv * transposed((const mat33_t&)Jl[1](0)));
		add(K, 3, 0, (const mat23_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));
		add(K, 3, 0, (const mat23_t&)Ja[1](3) * b[1]->I_inv * transposed((const mat33_t&)Ja[1](0)));

		//add(K, 3, 3, (const mat23_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat23_t&)Jl[0](3)));
		//add(K, 3, 3, (const mat23_t&)Jl[1](3) * b[1]->m_inv * transposed((const mat23_t&)Jl[1](3)));
		add(K, 3, 3, (const mat23_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat23_t&)Ja[0](3)));
		add(K, 3, 3, (const mat23_t&)Ja[1](3) * b[1]->I_inv * transposed((const mat23_t&)Ja[1](3)));
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};
	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};

public:
	vector_type r0, r1;
	vector_type parai;
	vector_type orth1j, orth2j;
};

template <class Real>
class MBFixedJoint : public MBJointBase<Real,2,6>
{
public:
	MBFixedJoint()
	{
		(mat_t&)Jl[0] = mat_t((Real)0);
		(mat_t&)Jl[1] = mat_t((Real)0);

		make_identity((mat_t&)Jl[0](3));
		make_identity((mat_t&)Jl[1](3));
		(mat_t&)Jl[1](3) = -(const mat_t&)Jl[1](3);

		(mat_t&)Ja[1](3) = mat_t((Real)0);
	}

	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}


	virtual void _precomp()
	{
		rij = tmul(b[0]->R, b[0]->r - b[1]->r);
		orth0i = vec_t(1, 0, 0);
		orth1j = cross_y(orth0i,(Real)1);
		orth2j = cross(orth0i, orth1j);

		orth0i = tmul(b[0]->R, orth0i);
		orth1i = tmul(b[0]->R, orth1j);
		orth1j = tmul(b[1]->R, orth1j);
		orth2j = tmul(b[1]->R, orth2j);
	}
	
	virtual void _update(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args)
	{
		typedef Vector<Real,3> vec3_t;
		typedef Matrix<Real,2,3> mat23_t;
		typedef Matrix<Real,3,3> mat33_t;

		Real dt = args.dt;
		Real betadtinv = ws.betadtinv;

		vec_t rij_(b[0]->r - b[1]->r);
		vec_t b0Rrij(b[0]->R * rij);
		vec_t orth0i_(b[0]->R * orth0i);
		vec_t orth1i_(b[0]->R * orth1i);
		vec_t orth1j_(b[1]->R * orth1j);
		vec_t orth2j_(b[1]->R * orth2j);

		(vec_t&)Ja[0](0) = cross(orth1i_,orth2j_);
		(vec_t&)Ja[0](1) = cross(orth0i_,orth1j_);
		(vec_t&)Ja[0](2) = cross(orth0i_,orth2j_);
		(vec_t&)Ja[1](0) = -(const vec_t&)Ja[0](0);
		(vec_t&)Ja[1](1) = -(const vec_t&)Ja[0](1);
		(vec_t&)Ja[1](2) = -(const vec_t&)Ja[0](2);

		C[0] = dot(orth2j_,orth1i_);
		C[1] = dot(orth1j_,orth0i_);
		C[2] = dot(orth2j_,orth0i_);
		
		rhs = (Real)0;
		//(vec_t&)rhs[0] += (mat_t&)Jl[0](0) * (b[0]->v + b[0]->m_inv*b[0]->f*dt);		
		(vec_t&)rhs[0] += (mat_t&)Ja[0](0) * (b[0]->w + b[0]->I_inv*b[0]->tau*dt);	
		//(vec_t&)rhs[0] += (mat_t&)Jl[1](0) * (b[1]->v + b[1]->m_inv*b[1]->f*dt);
		(vec_t&)rhs[0] += (mat_t&)Ja[1](0) * (b[1]->w + b[1]->I_inv*b[1]->tau*dt);


		(mat_t&)Ja[0](3) = skew(b0Rrij);

		vec_t rij_b0Rrij(rij_ - b0Rrij);
		(vec_t&)C[3] = rij_ - b0Rrij;

		//(vec_t&)rhs[3] += (mat_t&)Jl[0](3) * (b[0]->v + b[0]->m_inv*b[0]->f*dt);		
		(vec_t&)rhs[3] += b[0]->v + b[0]->m_inv*b[0]->f*dt;		
		//(vec_t&)rhs[3] += (mat_t&)Ja[0](3) * (b[0]->w + b[0]->I_inv*b[0]->w*dt);		
		(vec_t&)rhs[3] += cross(b0Rrij, b[0]->w + b[0]->I_inv*b[0]->tau*dt);	
		//(vec_t&)rhs[3] += (mat_t&)Jl[1](3) * (b[1]->v + b[1]->m_inv*b[1]->f*dt);
		(vec_t&)rhs[3] -= b[1]->v + b[1]->m_inv*b[1]->f*dt;
		//(vec_t&)rhs[3] += (mat_t&)Ja[1](3) * (b[1]->w + b[1]->I_inv*b[1]->tau*dt);
		
		rhs += betadtinv * C;

		K = (Real)0;

		//add(K, 0, 0, (const mat33_t&)Jl[0](0) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		//add(K, 0, 0, (const mat33_t&)Jl[1](0) * b[1]->m_inv * transposed((const mat33_t&)Jl[1](0)));
		add(K, 0, 0, (const mat33_t&)Ja[0](0) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));
		add(K, 0, 0, (const mat33_t&)Ja[1](0) * b[1]->I_inv * transposed((const mat33_t&)Ja[1](0)));

		//add(K, 3, 0, (const mat33_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](0)));
		//add(K, 3, 0, (const mat33_t&)Jl[1](3) * b[1]->m_inv * transposed((const mat33_t&)Jl[1](0)));
		add(K, 3, 0, (const mat33_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](0)));
		//add(K, 3, 0, (const mat33_t&)Ja[1](3) * b[1]->I_inv * transposed((const mat33_t&)Ja[1](0)));

		//add(K, 3, 3, (const mat33_t&)Jl[0](3) * b[0]->m_inv * transposed((const mat33_t&)Jl[0](3)));
		add(K, 3, 3, mat33_t(b[0]->m_inv, b[0]->m_inv, b[0]->m_inv));
		//add(K, 3, 3, (const mat23_t&)Jl[1](3) * b[1]->m_inv * transposed((const mat33_t&)Jl[1](3)));
		add(K, 3, 3, mat33_t(b[1]->m_inv, b[1]->m_inv, b[1]->m_inv));
		add(K, 3, 3, (const mat33_t&)Ja[0](3) * b[0]->I_inv * transposed((const mat33_t&)Ja[0](3)));
		//add(K, 3, 3, (const mat33_t&)Ja[1](3) * b[1]->I_inv * transposed((const mat33_t&)Ja[1](3)));
	}

	virtual void _preupdate_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};
	virtual void _update_by_impulse(MBWorkspace<Real>& ws, const MBUpdateArgs<Real>& args){};


	vector_type rij;
	vector_type orth0i, orth1i, orth1j, orth2j;
};

// spring
// damper
// 

struct MBForceElementTypeID { enum {
	Spring,
	Damper,
	TorsionalSpringDamper,
	SpringDamper,
};};

template <class Real>
class MBForceElement : public MBTypedef<Real>
{
	typedef MBForceElement this_t;
public:
	typedef MBBody<Real> body_type;

	virtual ~MBForceElement(){} 
	MBForceElement(){}
	MBForceElement(const this_t&){}
	this_t& operator=(const this_t&)
	{
		return *this;
	}

	virtual size_t body_number() const = 0;
	virtual const body_type& body(size_t ibody) const = 0;

	virtual void _precomp() = 0;
	virtual void _update() = 0;

};

template <class Real, size_t _N_B>
class MBForceElementBase : public MBForceElement<Real>
{
public:
	typedef MBBody<Real> body_type;

	static const size_t N_B = _N_B;

	virtual size_t body_number() const {return N_B;}
	virtual const body_type& body(size_t ibody) const {return *b[ibody];}

	body_type* b[N_B];
};


template <class Real>
class MBSkySpringDamper : public MBForceElementBase<Real,1>
{
public:
	MBSkySpringDamper() : k(0), c(0)
	{}

	void set_connecting_point(const vector_type& rp){ r0 = rp; }
	void set_sky_point(const vector_type& sp){ s = sp;}
	void set_natural_length(Real l){ l0 = l;}
	void set_spring_factor(Real _k){ k = _k;}
	void set_damper_factor(Real _c){ c = _c;}

	void link(body_type* b0)
	{
		b[0] = b0;
	}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r0 - b[0]->r);
	}

	virtual void _update()
	{
		vector_type b0Rr0 = b[0]->R * r0;
		vector_type dr = b[0]->r + b0Rr0 - s;
		Real l = length(dr);
		stroke = l - l0;

		vector_type dv = b[0]->v + cross(b[0]->w, b0Rr0); 

		vector_type dr_unit = unit_s(dr);
		vector_type f0 = - (k * stroke + c * dot(dr_unit,dv)) * dr_unit;
		vector_type tau0 = cross(b0Rr0, f0);
	
		b[0]->f += f0;
		b[0]->tau += tau0;
	}

public:
	vector_type r0;
	vector_type s;
	Real k;
	Real c;
	Real l0;
	Real stroke;
};

template <class Real>
class MBSpringDamper : public MBForceElementBase<Real,2>
{
public:
	MBSpringDamper() : k(0), c(0)
	{}

	void set_connecting_point(const vector_type& _r0, const vector_type& _r1)
	{
		r0 = _r0;
		r1 = _r1;
	}

	void set_natural_length(Real l) { l0 = l;}
	void set_spring_factor(Real _k){ k = _k;}
	void set_damper_factor(Real _c){ c = _c;}

	Real natural_length() const {return l0;}
	Real spring_factor() const {return k;}
	Real damper_factor() const {return c;}
	void connecting_point(vector_type& _r0, vector_type& _r1) const 
	{
		_r0 = r0;
		_r1 = r1;
	}

	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r0 - b[0]->r);
		r1 = tmul(b[1]->R, r1 - b[1]->r);
	}

	virtual void _update()
	{
		vector_type b0Rr0 = b[0]->R * r0;
		vector_type b1Rr1 = b[1]->R * r1;
		vector_type dr = b[0]->r + b0Rr0 - (b[1]->r + b1Rr1);
		Real l = length(dr);
		stroke = l - l0;

		vector_type dv = b[0]->v + cross(b[0]->w, b0Rr0) - (b[1]->v + cross(b[1]->w, b1Rr1)); 

		vector_type dr_unit = unit_s(dr);
		vector_type f0 = - (k * stroke + c * dot(dr_unit,dv)) * dr_unit;
		vector_type f1 = - f0;
		vector_type tau0 = cross(b0Rr0, f0);
		vector_type tau1 = cross(b1Rr1, f1);
	
		b[0]->f += f0;
		b[0]->tau += tau0;
		b[1]->f += f1;
		b[1]->tau += tau1;
	}

public:
	vector_type r0, r1;
	Real k;
	Real c;
	Real l0;
	Real stroke;
};

template <class Real>
class MBSkySpring : public MBForceElementBase<Real,1>
{
public:
	MBSkySpring() : k(0)
	{}

	void set_connecting_point(const vector_type& rp){ r0 = rp; }
	void set_sky_point(const vector_type& sp){ s = sp;}
	void set_natural_length(Real l){ l0 = l;}
	void set_factor(Real _k){ k = _k;}

	void link(body_type* b0)
	{
		b[0] = b0;
	}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r0 - b[0]->r);
	}

	virtual void _update()
	{
		vector_type b0Rr0 = b[0]->R * r0;
		vector_type dr = b[0]->r + b0Rr0 - s;
		Real l = length(dr);
		stroke = l - l0;

		vector_type dr_unit = unit_s(dr);
		vector_type f0 = - (k * stroke) * dr_unit;
		vector_type tau0 = cross(b0Rr0, f0);
	
		b[0]->f += f0;
		b[0]->tau += tau0;
	}

public:
	vector_type r0;
	vector_type s;
	Real k;
	Real l0;
	Real stroke;
};

template <class Real>
class MBSpring : public MBForceElementBase<Real,2>
{
public:
	MBSpring() : k(0)
	{}

	void set_connecting_point(const vector_type& _r0, const vector_type& _r1)
	{
		r0 = _r0;
		r1 = _r1;
	}

	void set_natural_length(Real l) { l0 = l;}
	void set_factor(Real _k){ k = _k;}

	Real natural_length() const {return l0;}
	Real factor() const {return k;}
	void connecting_point(vector_type& _r0, vector_type& _r1) const 
	{
		_r0 = r0;
		_r1 = r1;
	}

	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r0 - b[0]->r);
		r1 = tmul(b[1]->R, r1 - b[1]->r);
	}

	virtual void _update()
	{
		vector_type b0Rr0 = b[0]->R * r0;
		vector_type b1Rr1 = b[1]->R * r1;
		vector_type dr = b[0]->r + b0Rr0 - (b[1]->r + b1Rr1);
		Real l = length(dr);
		stroke = l - l0;

		vector_type dr_unit = unit_s(dr);
		vector_type f0 = - (k * stroke) * dr_unit;
		vector_type f1 = - f0;
		vector_type tau0 = cross(b0Rr0, f0);
		vector_type tau1 = cross(b1Rr1, f1);
	
		b[0]->f += f0;
		b[0]->tau += tau0;
		b[1]->f += f1;
		b[1]->tau += tau1;
	}

public:
	vector_type r0, r1;
	Real k;
	Real l0;
	Real stroke;
};

template <class Real>
class MBSkyDamper : public MBForceElementBase<Real,1>
{
public:
	MBSkyDamper() : c(0)
	{}

	void set_connecting_point(const vector_type& rp){ r0 = rp; }
	void set_sky_point(const vector_type& sp){ s = sp;}
	void set_factor(Real _c){ c = _c;}

	void link(body_type* b0)
	{
		b[0] = b0;
	}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r0 - b[0]->r);
	}

	virtual void _update()
	{
		vector_type b0Rr0 = b[0]->R * r0;
		vector_type dr = b[0]->r + b0Rr0 - s;

		vector_type dv = b[0]->v + cross(b[0]->w, b0Rr0); 

		vector_type dr_unit = unit_s(dr);
		vector_type f0 = - (c * dot(dr_unit,dv)) * dr_unit;
		vector_type tau0 = cross(b0Rr0, f0);
	
		b[0]->f += f0;
		b[0]->tau += tau0;
	}

public:
	vector_type r0;
	vector_type s;
	Real c;
};

template <class Real>
class MBDamper : public MBForceElementBase<Real,2>
{
public:
	MBDamper() : c(0)
	{}

	void set_connecting_point(const vector_type& _r0, const vector_type& _r1)
	{
		r0 = _r0;
		r1 = _r1;
	}

	void set_factor(Real _c){ c = _c;}

	Real factor() const {return c;}
	void connecting_point(vector_type& _r0, vector_type& _r1) const 
	{
		_r0 = r0;
		_r1 = r1;
	}


	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	virtual void _precomp()
	{
		r0 = tmul(b[0]->R, r0 - b[0]->r);
		r1 = tmul(b[1]->R, r1 - b[1]->r);
	}

	virtual void _update()
	{
		vector_type b0Rr0 = b[0]->R * r0;
		vector_type b1Rr1 = b[1]->R * r1;
		vector_type dr = b[0]->r + b0Rr0 - (b[1]->r + b1Rr1);

		vector_type dv = b[0]->v + cross(b[0]->w, b0Rr0) - (b[1]->v + cross(b[1]->w, b1Rr1)); 

		vector_type dr_unit = unit_s(dr);
		vector_type f0 = - (c * dot(dr_unit,dv)) * dr_unit;
		vector_type f1 = - f0;
		vector_type tau0 = cross(b0Rr0, f0);
		vector_type tau1 = cross(b1Rr1, f1);
	
		b[0]->f += f0;
		b[0]->tau += tau0;
		b[1]->f += f1;
		b[1]->tau += tau1;
	}

public:
	vector_type r0, r1;
	Real c;
};


template <class Real>
class MBSkyTorsionalSpringDamper : public MBForceElementBase<Real,1>
{
	MBSkyTorsionalSpringDamper()
	{}
};

template <class Real>
class MBTorsionalSpringDamper : public MBForceElementBase<Real,2>
{
public:
	MBTorsionalSpringDamper() : theta0(0), k(0), c(0) 
	{
	}

	void set_axis(const vector_type& axis)
	{
		axis0 = axis;
		orth0 = orthogonal_axis(axis);
		//vec_t d;
		//orthogonal_axes(d, orth0, axis);
	}

	void set_spring_factor(Real _k){ k = _k;}
	void set_damper_factor(Real _c){ c = _c;}

	void set_natural_angle(Real theta0)
	{
		this->theta0 = theta0;
	}

	Real spring_factor() const {return k;}
	Real damper_factor() const {return c;}
	Real natural_angle() const {return theta0;}
	const vector_type& axis() const {return axis0;}


	void link(body_type* b0, body_type* b1)
	{
		b[0] = b0;
		b[1] = b1;
	}

	virtual void _precomp()
	{
		axis1 = tmul(b[1]->R, axis0);
		axis0 = tmul(b[0]->R, axis0);

		orth1 = tmul(b[1]->R, orth0);
		orth0 = tmul(b[0]->R, orth0);
	}

	virtual void _update()
	{
		vec_t b0Raxis0 = b[0]->R * axis0;
		vec_t b1Raxis1 = b[1]->R * axis1;
		vec_t b0Rorth0 = b[0]->R * orth0;
		vec_t b1Rorth1 = b[1]->R * orth1;


		Real w0 = dot(b0Raxis0, b[0]->w);
		Real w1 = dot(b1Raxis1, b[1]->w);

		dtheta = angle(b1Rorth1, b0Rorth0, b0Raxis0);
		//dtheta = angle(b0Rorth0, b1Rorth1, b0Raxis0);

		Real t = - (k * (dtheta - theta0) + c * (w0 - w1)); 
		vec_t tau0 =   t * b0Raxis0;
		vec_t tau1 = - t * b1Raxis1;

		b[0]->tau += tau0;
		b[1]->tau += tau1;
	}

public:
	vector_type axis0, orth0;
	vector_type axis1, orth1;
	Real k;
	Real c;
	Real dtheta;
	Real theta0;
};

// 基準のbody3とbody1にtorsion、body2にforceに
template <class Real>
class MBTorsionToForceSpringDamper : public MBForceElementBase<Real,2>
{
	MBTorsionToForceSpringDamper()
	{}
};

template <class Real>
class MBTorsionTo2ForceSpringDamper : public MBForceElementBase<Real,2>
{
	MBTorsionTo2ForceSpringDamper()
	{}
};

template <class Real>
class MBWorld : public MBTypedef<Real>
{
public:
	typedef MBUpdateArgs<Real> update_args_type;

	typedef MBBody<Real> body_type;
	typedef MBJoint<Real> joint_type;
	typedef MBForceElement<Real> force_element_type;
	//typedef TaskList<body_type, MBTire<Real> > bodies_type;
	//typedef TaskList<joint_type, MBFixedJoint<Real> > joints_type;
	//typedef TaskList<force_element_type, MBTorsionalSpringDamper<Real> > force_elements_type;
	typedef TaskVector<body_type, MBTire<Real> > bodies_type;
	typedef TaskVector<joint_type, MBFixedJoint<Real> > joints_type;
	typedef TaskVector<force_element_type, MBTorsionalSpringDamper<Real> > force_elements_type;

	typedef Function<void(const update_args_type&)> on_add_force_type;
	typedef Function<void(const update_args_type&)> on_rewrite_velocity_type;


	MBWorld()
	{
		g = (Real)0;

		stop_watch.init();
	}

	~MBWorld()
	{
		release();
	}

	void init(size_t n_body, size_t n_joint, size_t n_fe)
	{
		bodies.assign(n_body);
		joints.assign(n_joint);
		force_elements.assign(n_fe);
	}
	
#define CREATE_BODY(name) \
	void create(name<Real>*& body_ptr)\
	{\
		name<Real> body; bodies.push_back(body); body_ptr = (name<Real>*)&bodies.back();\
	}\

	CREATE_BODY(MBBody)
	CREATE_BODY(MBTire)

#undef CREATE_BODY

#define CREATE_JOINT(name) \
	void create(name<Real>*& joint_ptr)\
	{\
		name<Real> joint; joints.push_back(joint); joint_ptr = (name<Real>*)&joints.back();\
	}\

	CREATE_JOINT(MBDistanceJoint)
	CREATE_JOINT(MBAngleJoint)
	CREATE_JOINT(MBBallPlaneJoint)
	CREATE_JOINT(MBSkyBallJoint)
	CREATE_JOINT(MBBallJoint)
	CREATE_JOINT(MBUniversalJoint)
	CREATE_JOINT(MBSkySliderJoint)
	CREATE_JOINT(MBSliderJoint)
	CREATE_JOINT(MBSkyHingeJoint)
	CREATE_JOINT(MBHingeJoint)
	CREATE_JOINT(MBFixedJoint)

#undef CREATE_JOINT


#define CREATE_FORCE_ELEMENT(name) \
	void create(name<Real>*& obj_ptr)\
	{\
		name<Real> obj; force_elements.push_back(obj); obj_ptr = (name<Real>*)&force_elements.back();\
	}\

	CREATE_FORCE_ELEMENT(MBSkySpringDamper)
	CREATE_FORCE_ELEMENT(MBSpringDamper)
	CREATE_FORCE_ELEMENT(MBSkySpring)
	CREATE_FORCE_ELEMENT(MBSpring)
	CREATE_FORCE_ELEMENT(MBSkyDamper)
	CREATE_FORCE_ELEMENT(MBDamper)
	CREATE_FORCE_ELEMENT(MBSkyTorsionalSpringDamper)
	CREATE_FORCE_ELEMENT(MBTorsionalSpringDamper)
#undef CREATE_FORCE_ELEMENT

	void precomp()
	{
		using namespace std;
		
		size_t N_b = 0;
		size_t ibody=0;
		for(typename bodies_type::iterator it=bodies.begin(); it!=bodies.end(); ++it)
		{
			body_type& b = *it;
			b.id = ibody;

			++ibody;
			N_b += 6;
		}

		size_t N_c = 0;
		size_t N_j = joints.size();
		for(typename joints_type::iterator it=joints.begin(); it!=joints.end(); ++it)
		{
			joint_type& c = *it;
			c._precomp();
			
			N_c += c.constraint_number();
		}

		//size_t N_f = force_elements.size();
		for(typename force_elements_type::iterator it=force_elements.begin(); it!=force_elements.end(); ++it)
		{
			force_element_type& f = *it;
			f._precomp();
		}

		memory.release();

		memory.Minv.assign(N_b, N_b, 0);
		memory.J.assign(N_c, N_b, 0);
		memory.JMinv.assign(N_c, N_b, 0);
		memory.W.assign(6, 6, 0);
		memory.K.assign(N_c,N_c, 0);
		memory.rhs.assign(N_c, 0);
		memory.lambda.assign(N_c, 0);
		memory.w_c.assign(N_c, 0);
		memory.a.assign(N_b, 0);
		memory.v.assign(N_b, 0);
		memory.p.assign((N_b/6)*13, 0);
		memory.f.assign(N_b, 0);
		memory.Fc.assign(N_b, 0);
		
		memory.joint_pairs.assign(N_j, N_j);

		memory.constraint_numbers.assign(N_j, 0);
		memory.K_map.assign(N_j, N_j);



		{
			size_t j0=0;
			for(typename joints_type::iterator it0=joints.begin(); it0!=joints.end(); ++it0)
			{
				const joint_type& c0 = *it0;
				size_t j1=0;
				auto it0next = it0;
				++it0next;
				for(typename joints_type::iterator it1=joints.begin(); it1!=it0next; ++it1)
				{
					const joint_type& c1 = *it1;
					auto& joint_pair = memory.joint_pairs(j0,j1);

#if 0
					if(j0 == j1)
					{
						++j1;
						continue;
					}
#endif			
					for(size_t i0=0; i0<c0.body_number(); ++i0)
					{
						for(size_t i1=0; i1<c1.body_number(); ++i1)
						{
							if(&c0.body(i0) == &c1.body(i1))
							{
								Memory::JointPair::Body body;
								body.joint0.body_index = i0;
								body.joint1.body_index = i1; 
								body.body = &c0.body(i0); 
								joint_pair.bodies.push_back(body);
							}
						}
					}

					++j1;
				}

				++j0;
			}
		}

		{
			size_t j=0;
			size_t iconstraint = 0;
			for(typename joints_type::iterator it=joints.begin(); it!=joints.end(); ++it)
			{
				joint_type& c = *it;
				memory.constraint_numbers[j] = c.constraint_number();
				iconstraint += c.constraint_number();
				++j;
			}

			for(size_t j0=0; j0<N_j; ++j0)
			{
				for(size_t j1=0; j1<=j0; ++j1)
				{
					const auto& joint_pair = memory.joint_pairs(j0,j1);
					for(size_t ib=0; ib<joint_pair.bodies.size(); ++ib)
					{
						memory.K_map(j0,j1).push(joint_pair.bodies[ib].body->id);
					}
				}
			}
		}
	}


	void update(const MBUpdateArgs<Real>& args)
	{
#define MAT(a,m,i,j) a[(i)*(m) + (j)]
#define MAT0(a,m,i) a[(i)*(m)]

		using namespace details;

		//static int debug_counter=0;
		//static std::ofstream fout("debug__multibody4.txt");
		//if(debug_counter == 0) write_memory(fout, *this);

		Real dt = args.dt;
		ws.dt_inv = (Real)1/dt; 
		ws.betadtinv = args.erp*ws.dt_inv; 

		stop_watch.start();
		
		for(typename bodies_type::iterator it=bodies.begin(); it!=bodies.end(); ++it)
		{
			body_type& b = *it;

			b.m_inv = (Real)1 / b.m;
			b.I = b.R * mat_t(b.I_b._0, b.I_b._1, b.I_b._2) * transposed(b.R);
			inverse(b.I_inv, b.I);
			
			b.f = (Real)0;
			b.tau = (Real)0;
		}

		for(typename force_elements_type::iterator it=force_elements.begin(); it!=force_elements.end(); ++it)
		{
			force_element_type& f = *it;
			f._update();
		}

		on_add_force(args);

		size_t ibody = 0;
		size_t n_body = 0;
		for(typename bodies_type::iterator it=bodies.begin(); it!=bodies.end(); ++it)
		{
			body_type& b = *it;
			size_t ibody6 = ibody*6;
			size_t ibody13 = ibody*13;

			b.f += b.m * g;
			b.tau += - cross(b.w, b.I * b.w); 

			memory.Minv(ibody6 + 0, ibody6 + 0) = b.m_inv;
			memory.Minv(ibody6 + 1, ibody6 + 1) = b.m_inv;
			memory.Minv(ibody6 + 2, ibody6 + 2) = b.m_inv;
	
			for(size_t i=0; i<3; ++i)
			{
				for(size_t j=0; j<3; ++j)
				{
					memory.Minv(ibody6 + 3 + i, ibody6 + 3 + j) = b.I_inv(i,j);
				}
			}

			(vector_type&)		memory.f[ibody6    ] = b.f;
			(vector_type&)		memory.f[ibody6 + 3] = b.tau;

			//(vector_type&)		memory.v[ibody6  ] = b.v;
			//(vector_type&)		memory.v[ibody6+3] = b.w;

			// list -> memory
			(vector_type&)		memory.p[ibody13    ] = b.v;
			(vector_type&)		memory.p[ibody13 + 3] = b.w;
			(vector_type&)		memory.p[ibody13 + 6] = b.r;
			(quaternion_type&)	memory.p[ibody13 + 9] = b.q;

			++ibody;
		}
		n_body = ibody;


		elapsed_times[0] = stop_watch.elapsed_time();

		
		stop_watch.start();
		size_t iconstraint=0;
		size_t n_constraint=0;
		size_t n_joint=0;

		static int count=0;
		for(typename joints_type::iterator it=joints.begin(); it!=joints.end(); ++it)
		{
			joint_type& c = *it;
			c._update(ws, args);

			size_t n_b = c.body_number();
			size_t n_c = c.constraint_number();

			// Jのセット
			for(size_t ib=0; ib<n_b; ++ib)
			{
				Real* Jl = c.Jl_(ib);
				Real* Ja = c.Ja_(ib);
				size_t ibody = c.body(ib).id;
				size_t ibody6 = ibody*6;

				for(size_t ic=0; ic<n_c; ++ic)
				{
					size_t ic3 = ic*3;
					size_t icc = iconstraint + ic;
					Real* Jicc = &memory.J(icc, ibody6); 
					Jicc[ 0] = Jl[ic3    ];
					Jicc[ 1] = Jl[ic3 + 1];
					Jicc[ 2] = Jl[ic3 + 2];
					Jicc[ 3] = Ja[ic3    ];
					Jicc[ 4] = Ja[ic3 + 1];
					Jicc[ 5] = Ja[ic3 + 2];
				}
			}

			// rhsのセット
			for(size_t ic=0; ic<n_c; ++ic)
			{
				size_t icc = iconstraint + ic;
				Real* rhs = c.rhs_();
				memory.rhs[icc] = rhs[ic];
			}

			iconstraint += n_c;
			++n_joint;
		}
		n_constraint = iconstraint;


		stop_watch.start();
	
#if 1
		// jointの組み合わせで、KのL部分を作成
		{
			size_t j0 = 0;
			size_t iconstraint0 = 0;
			for(typename joints_type::iterator it0=joints.begin(); it0!=joints.end(); ++it0)
			{
				joint_type& c0 = *it0;
				size_t n_c0 = c0.constraint_number();
				
				size_t j1 = 0;
				size_t iconstraint1 = 0;
				auto it0next = it0;
				++it0next;
				for(typename joints_type::iterator it1=joints.begin(); it1!=it0next; ++it1)
				{
					joint_type& c1 = *it1;
					size_t n_c1 = c1.constraint_number();

					if(j0 == j1)
					{
						Real* K = c1.K_();
						size_t n_c = c1.constraint_number();

						for(size_t i0=0; i0<n_c; ++i0)
						{
							Real* Kic0ii0 = &memory.K(iconstraint0+i0);
							Real* Ki0n_c = &K[i0*n_c];
							for(size_t i1=0; i1<n_c; ++i1)
							{
								Kic0ii0[iconstraint1+i1] = Ki0n_c[i1];  
							}
						}

					    iconstraint1 += n_c1;
						++j1;
						continue;
					}

					for(size_t i0=0; i0<n_c0; ++i0)
					{
						Real* Kic0ii0 = &memory.K(iconstraint0+i0);
						for(size_t i1=0; i1<n_c1; ++i1)
						{
							Kic0ii0[iconstraint1+i1] = 0;
						}
					}

					auto& joint_pair = memory.joint_pairs(j0,j1);
					size_t n_b = joint_pair.bodies.size();

					for(size_t ib=0; ib<n_b; ++ib)
					{
						auto& b = joint_pair.bodies[ib];

						Real* Jl0b = c0.Jl_(b.joint0.body_index);
						Real* Ja0b = c0.Ja_(b.joint0.body_index);
						Real* Jlb1 = c1.Jl_(b.joint1.body_index);
						Real* Jab1 = c1.Ja_(b.joint1.body_index);

						// W = J0 m^-1
						for(size_t i0=0; i0<n_c0; ++i0)
						{
							Real* Wi0 = &memory.W(i0);
							Real* Jl0b3i0 = &MAT0(Jl0b,3,i0);
							Real* Ja0b3i0 = &MAT0(Ja0b,3,i0);

							Wi0[0] = Jl0b3i0[0] * b.body->m_inv;
							Wi0[1] = Jl0b3i0[1] * b.body->m_inv;
							Wi0[2] = Jl0b3i0[2] * b.body->m_inv;

							for(size_t i1=0; i1<3; ++i1)
							{
								Real sum(0);
								for(size_t i2=0; i2<3; ++i2)
								{
									sum += Ja0b3i0[i2] * b.body->I_inv(i2,i1);
								}
								Wi0[i1+3] = sum;
							}
						}

						// (J0 M^-1) J1^T
						for(size_t i0=0; i0<n_c0; ++i0)
						{
							Real* Kic0ii0 = &memory.K(iconstraint0+i0);
							Real* Wi0 = &memory.W(i0);
							for(size_t i1=0; i1<n_c1; ++i1)
							{
								Real* Jlb13i1 = &MAT0(Jlb1,3,i1);
								Real sum(0);
								for(size_t i2=0; i2<3; ++i2)
								{
									sum += Wi0[i2] * Jlb13i1[i2];
								}
								Kic0ii0[iconstraint1+i1] += sum;
							}

							for(size_t i1=0; i1<n_c1; ++i1)
							{
								Real* Jab13i1 = &MAT0(Jab1,3,i1);
								Real sum(0);
								for(size_t i2=0; i2<3; ++i2)
								{
									sum += Wi0[i2+3] * Jab13i1[i2];
								}
								Kic0ii0[iconstraint1+i1] += sum;
							}
						}
					}
		
					iconstraint1 += n_c1;
					++j1;
				}

				iconstraint0 += n_c0;
				++j0;
			}
		}
#endif

#if 0
		// JMinv
		{
			memory.JMinv = (Real)0;
			size_t j = 0;
			size_t iconstraint = 0;
			for(typename joints_type::iterator it=joints.begin(); it!=joints.end(); ++it)
			{
				joint_type& c = *it;
				size_t n_c = c.constraint_number();
				
				size_t n_b = c.body_number();

				for(size_t ib=0; ib<n_b; ++ib)
				{
					size_t ibody = c.body(ib).id;
					size_t ibody6 = ibody*6;

					for(size_t ic=0; ic<n_c; ++ic)
					{
						memory.JMinv(iconstraint+ic, ibody6  ) = memory.J(iconstraint+ic, ibody6  ) * memory.Minv(ibody6  , ibody6  );
						memory.JMinv(iconstraint+ic, ibody6+1) = memory.J(iconstraint+ic, ibody6+1) * memory.Minv(ibody6+1, ibody6+1);
						memory.JMinv(iconstraint+ic, ibody6+2) = memory.J(iconstraint+ic, ibody6+2) * memory.Minv(ibody6+2, ibody6+2);

						for(size_t k=0; k<3; ++k)
						{
							memory.JMinv(iconstraint+ic, ibody6+3) += memory.J(iconstraint+ic, ibody6+3) * memory.Minv(ibody6+3+k, ibody6+3);
							memory.JMinv(iconstraint+ic, ibody6+4) += memory.J(iconstraint+ic, ibody6+4) * memory.Minv(ibody6+3+k, ibody6+4);
							memory.JMinv(iconstraint+ic, ibody6+5) += memory.J(iconstraint+ic, ibody6+5) * memory.Minv(ibody6+3+k, ibody6+5);
						}
					}
				}

				iconstraint += n_c;
				++j;
			}
		}

		// JMinv J^T
		{
			memory.K = (Real)0;
			
			size_t iconstraint0 = 0;
			for(size_t j0=0; j0<n_joint; ++j0)
			{
				size_t iconstraint1 = 0;
				size_t n_c0 = memory.constraint_numbers[j0];
				for(size_t j1=0; j1<=j0; ++j1)
				{
					size_t n_c1 = memory.constraint_numbers[j1];
					const auto& body_indices = memory.K_map(j0,j1);
					size_t n_b = body_indices.size();

					for(size_t ib=0; ib<n_b; ++ib)
					{
						size_t ibody = body_indices[ib];
						size_t ibody6 = ibody*6;

						for(size_t ic0=0; ic0<n_c0; ++ic0)
						{
							Real* Kicc0ic0 = &memory.K(iconstraint0+ic0);
							Real* JMinvicc0ic0 = &memory.JMinv(iconstraint0+ic0);
							for(size_t ic1=0; ic1<n_c1; ++ic1)
							{
								for(size_t ibb=0; ibb<6; ++ibb)
								{
									Kicc0ic0[iconstraint1+ic1] += JMinvicc0ic0[ibody6+ibb] * memory.J(iconstraint1+ic1, ibody6+ibb);
								}
							}
						}
					}

					iconstraint1 += n_c1;
				}

				iconstraint0 += n_c0; 
			}
		}

#endif

		//memory.K = (memory.J * memory.Minv) * transposed(memory.J);

#if 0
		{
			//size_t n_tmp = n_body*n_body*n_body*6*6*6;
			size_t n_tmp = n_body*n_body*n_body*6*6*6;
			size_t sum(0);
			for(size_t i=0; i<n_tmp; ++i)
			{
				sum = i;
			}
			memory.a[0] = sum;
		}
#endif

		elapsed_times[1] = stop_watch.elapsed_time();


		//cout << " -- memory.K, memory.rhs -- " << endl;
		//cout << memory.K << endl;
		//cout << memory.rhs << endl;


		memory.lambda = memory.rhs;

		stop_watch.start();
		// K lambda = rhs
		solve_le_by_cholesky(
			memory.lambda.data(), 
			memory.K.data(),
			memory.K.col_size(),
			memory.K.col_size(),
			memory.w_c.data()
			);
#if 0
		solve_le_by_lu_crout(
			memory.lambda.data(), 
			memory.K.data(),
			memory.K.col_size(),
			memory.K.col_size()
			);
#endif
		elapsed_times[2] = stop_watch.elapsed_time();


		// v(n+1) = v(n) + M^-1 (f dt - J^T lambda dt) 
		stop_watch.start();
	
		{
			size_t j = 0;
			size_t iconstraint = 0;

			memory.Fc = (Real)0;

			for(typename joints_type::iterator it=joints.begin(); it!=joints.end(); ++it)
			{
				joint_type& c = *it;
				size_t n_c = c.constraint_number();
				size_t n_b = c.body_number();

				for(size_t ib=0; ib<n_b; ++ib)
				{
					Real* Jl = c.Jl_(ib);
					Real* Ja = c.Ja_(ib);
					size_t ibody = c.body(ib).id;
					size_t ibody6 = ibody*6;
					for(size_t ic=0; ic<n_c; ++ic)
					{
						Real* Jlic = &MAT0(Jl,3,ic);
						Real* Jaic = &MAT0(Ja,3,ic);
						const Real& lambda = memory.lambda[iconstraint+ic];
						memory.Fc[ibody6 + 0] += Jlic[0] * lambda;  
						memory.Fc[ibody6 + 1] += Jlic[1] * lambda;  
						memory.Fc[ibody6 + 2] += Jlic[2] * lambda;  
						memory.Fc[ibody6 + 3] += Jaic[0] * lambda;  
						memory.Fc[ibody6 + 4] += Jaic[1] * lambda;  
						memory.Fc[ibody6 + 5] += Jaic[2] * lambda;  
					}
				}

				iconstraint += n_c;
			}
		}

#if 0
		transposed_matrix_product(
			memory.Fc.data(),
			memory.J.data(),
			memory.J.col_size(),
			memory.J.row_size(),
			memory.J.col_size(),
			memory.lambda.data()
			);
#endif

		elapsed_times[3] = stop_watch.elapsed_time();

		memory.Fc *= ws.dt_inv;

		stop_watch.start();
		
		{
			size_t ibody = 0;
			for(typename bodies_type::iterator it=bodies.begin(); it!=bodies.end(); ++it)
			{
				body_type& b = *it;
				size_t ibody6 = ibody*6;
		
				memory.a[ibody6    ] = b.m_inv * (memory.f[ibody6    ] - memory.Fc[ibody6    ]);
				memory.a[ibody6 + 1] = b.m_inv * (memory.f[ibody6 + 1] - memory.Fc[ibody6 + 1]);
				memory.a[ibody6 + 2] = b.m_inv * (memory.f[ibody6 + 2] - memory.Fc[ibody6 + 2]);

				for(size_t i=0; i<3; ++i)
				{
					Real sum(0);
					Real* aibody63 = &memory.a[ibody6 + 3];
					Real* fibody63 = &memory.f[ibody6 + 3];
					Real* Fcibody63 = &memory.Fc[ibody6 + 3];
					for(size_t j=0; j<3; ++j)
					{
						sum += b.I_inv(i,j) * (fibody63[j] - Fcibody63[j]);
					}
					aibody63[i] = sum;
				}

				++ibody;
			}
		}
		
		//memory.a = memory.Minv * (memory.f - memory.Fc);
		
		elapsed_times[4] = stop_watch.elapsed_time();
		

		// v,w,r,qの積分
		for(size_t ibody=0; ibody<n_body; ++ibody)
		{
			size_t ibody6=ibody*6;
			size_t ibody13=ibody*13;
			vector_type& v = (vector_type&)memory.p[ibody13  ];
			vector_type& w = (vector_type&)memory.p[ibody13+3];
			v += dt * (const vector_type&)memory.a[ibody6  ];
			w += dt * (const vector_type&)memory.a[ibody6+3];
		}

		on_rewrite_velocity(args);

		for(size_t ibody=0; ibody<n_body; ++ibody)
		{
			size_t ibody13=ibody*13;
			const vector_type& v = (vector_type&)memory.p[ibody13  ];
			const vector_type& w = (vector_type&)memory.p[ibody13+3];
			vector_type& r = (vector_type&)memory.p[ibody13+6];
			quaternion_type& q = (quaternion_type&)memory.p[ibody13+9];
			r += dt * v;
			q += dt * diff(w,q);
			normalize(q);
		}

		/*
		memory -> list
		*/
		ibody = 0;
		for(typename bodies_type::iterator it=bodies.begin(); it!=bodies.end(); ++it)
		{
			body_type& b = *it;
			size_t ibody13 = ibody*13;
			size_t ibody6 = ibody*6;

			b.a		= (const vector_type&)	memory.a[ibody6  ];
			b.alpha = (const vector_type&)	memory.a[ibody6+3];
			b.v = (const vector_type&)		memory.p[ibody13  ];
			b.w = (const vector_type&)		memory.p[ibody13+3];
			b.r = (const vector_type&)		memory.p[ibody13+6];
			b.q = (const quaternion_type&)	memory.p[ibody13+9];
			
			convert(b.R, b.q);

			++ibody;
		}

		//++debug_counter;
#undef MAT
#undef MAT0
	}

	void update(Real dt)
	{
		MBUpdateArgs<Real> args;
		args.dt = dt;
		args.erp = (Real)1;
		update(args);
	}

	void update_by_impulse(const MBWorldUpdateByImpulseArgs<Real>& args)
	{
		Real dt = args.dt;
		ws.dt_inv = (Real)1/dt; 
		ws.betadtinv = args.erp*ws.dt_inv; 

		for(typename bodies_type::iterator it=bodies.begin(); it!=bodies.end(); ++it)
		{
			body_type& b = *it;

			b.m_inv = (Real)1 / b.m;
			b.I = b.R * diag(b.I_b) * transposed(b.R);
			inverse(b.I_inv, b.I);
			
			b.f = b.f_e + b.m * g;
			b.tau = b.tau_e - cross(b.w, b.I * b.w); 
			b.v += dt * b.m_inv * (b.f); 
			b.w += dt * b.I_inv * (b.tau); 
		}

		// 一度のみ計算
		for(typename joints_type::iterator it=joints.begin(); it!=joints.end(); ++it)
		{
			// K = J M J^T = J_v M' J_v^T + J_w I J_w^T
			joint_type& c = *it;
			c._preupdate_by_impulse(ws, args);
		}

		// 繰り返し計算
		size_t count = 0;
		//cout << "-- iteration start --" << endl;
		while(count++ < args.iteration_times)
		{
			ws.rhs_sum = (Real)0;
	
			for(typename joints_type::iterator it=joints.begin(); it!=joints.end(); ++it)
			{
				// v(n + 1) = v(n) + dv
				// dv := M^-1 f dt - M^-1 J^T impulse(v)
				// impulse(v) := lambda dt = K^-1 rhs(v)
				// K := J M^-1 J^T
				// rhs(v) := J (1 + e) v(n) + rhs_c
				// rhs_c := J M^-1 f dt + beta C
				// 
				// dv_0 = M^-1 f dt
				// impulse_0 = impulse(v)
				// dv_1 = dv_0 - M^-1 J^T impulse_0
				// impulse_1 = impulse_0 + impulse(dv_1) 
				// dv_2 = dv_1 - M^-1 J^T impulse_1
				// impulse_2 = impulse_0 + impulse(dv_2)
				// dv_3 = dv_2 - M^-1 J^T impulse_2
				joint_type& c = *it;
				c._update_by_impulse(ws, args);
			}

			if(ws.rhs_sum < args.eps_rhs_sum)
			{
				break;
			}
		}

		for(typename bodies_type::iterator it=bodies.begin(); it!=bodies.end(); ++it)
		{
			body_type& b = *it;
			b.r += dt*b.v;
			b.q += dt*diff(b.w,b.q);
			normalize(b.q);
			convert(b.R, b.q);
		}
	}

	void update_by_impulse(Real dt)
	{
		MBWorldUpdateByImpulseArgs<Real> args;
		args.dt = dt;
		args.erp = (Real)1;
		args.iteration_times = 10;
		args.eps_rhs_sum = 1.;
		update_by_impulse(args);
	}

	void release()
	{
		memory.release();

		bodies.release();
		joints.release();
		force_elements.release();
	}

	void set_gravity(const vector_type& g)
	{
		this->g = g;
	}

	size_t total_contraint_number() const
	{
		size_t N_c = 0;
		for(typename joints_type::const_iterator it=joints.begin(); it!=joints.end(); ++it)
		{
			const joint_type& c = *it;
			N_c += c.constraint_number();
		}

		return N_c;
	}

public:
	vector_type g;

	bodies_type bodies;
	joints_type joints;
	force_elements_type force_elements;

	struct Memory
	{
		DMatrix<Real> Minv;
		DMatrix<Real> W;
		DMatrix<Real> K;
		DVector<Real> lambda;
		DVector<Real> w_c;
		DMatrix<Real> JMinv;
		DMatrix<Real> J;
		DVector<Real> rhs;
		DVector<Real> f;
		DVector<Real> a;
		DVector<Real> v;
		DVector<Real> p;
		DVector<Real> Fc;

		DVector<size_t> constraint_numbers;
		
		DMatrix<Stack<size_t, size_t[3]> > K_map; 



		struct JointPair
		{
			struct Joint
			{
				size_t body_index;
			};

			struct Body
			{
				Joint joint0;
				Joint joint1;
				const MBBody<Real>* body;
			};

			std::vector<Body> bodies;
		};

		struct JointPairs
		{
			JointPairs(){}

			void assign(size_t nrow, size_t ncol)
			{
				_ncol = ncol;
				joint_pairs.assign(nrow*ncol, JointPair());
			}

			JointPair& operator()(size_t irow, size_t icol)
			{
				return joint_pairs[irow*_ncol+icol];
			}

			void clear()
			{
				joint_pairs.clear();
			}

		private:
			std::vector<JointPair> joint_pairs;
			size_t _ncol;
		};

		JointPairs joint_pairs;

		void release()
		{
			Minv.clear();
			K.clear();
			lambda.clear();
			w_c.clear();
			JMinv.clear();
			J.clear();
			rhs.clear();
			f.clear();

			joint_pairs.clear();
			constraint_numbers.clear();
			K_map.clear();
		}
	};


	MBWorkspace<Real> ws;


	Memory memory;


	on_add_force_type on_add_force;
	on_rewrite_velocity_type on_rewrite_velocity;


	PerformanceCounterStopWatch<double> stop_watch;

	double elapsed_times[10];


};


typedef MBSkyBallJoint<double> mb_sky_ball_joint3d_t;

template <class Char, class Traits, class Real>
void write_containers(std::basic_ostream<Char,Traits>& os, const MBWorld<Real>& w)
{
	using namespace std;
	os << "-- Bodies -- " << endl;
	size_t N_b = w.bodies.size();
	os << "N_b:" << N_b << endl;
	for(auto it=w.bodies.begin(); it!=w.bodies.end(); ++it)
	{
		os << "mass: " << it->m << endl;
		os << "inertia: " << it->I_b << endl;
		os << "position: " << it->r << endl;
		os << "orientation: " << it->q << endl;
	}
	
	size_t N_c = 0;
	size_t N_j = w.joints.size();
	os << "-- Joints -- " << endl;
	os << "N_j:" << N_j << endl;
	for(auto it=w.joints.begin(); it!=w.joints.end(); ++it)
	{
		os << "constraint_number: " << it->constraint_number() << endl;
		os << "body_number: " << it->body_number() << endl;
		N_c += it->constraint_number();
	}
	os << "N_c:" << N_c << endl;
}


template <class Char, class Traits, class Real>
void write_memory(std::basic_ostream<Char,Traits>& os, const MBWorld<Real>& w)
{
	using namespace std;
	os << "--- memory ---" << endl;
	os << "Minv:" << endl;
	os << w.memory.Minv << endl;
	os << "J:" << endl;
	os << w.memory.J << endl;
	os << "J^T:" << endl;
	os << transposed(w.memory.J) << endl;
	os << "K:" << endl;
	os << w.memory.K << endl;
	os << "rhs:" << endl;
	os << w.memory.rhs << endl;
	os << "lambda:" << endl;
	os << w.memory.lambda << endl;
	os << "f:" << endl;
	os << w.memory.f << endl;
	os << "Fc:" << endl;
	os << w.memory.Fc << endl;
	os << "a:" << endl;
	os << w.memory.a << endl;
}


}}}

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif
