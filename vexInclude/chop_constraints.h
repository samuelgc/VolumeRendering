/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *  	Side Effects Software Inc
 *  	123 Front Street West
 *  	Toronto, Ontario
 *  	Canada   M5V 3E7
 *  	416-504-9876
 *
 * NAME:    chop_constraints.h
 */
#ifndef __chop_constraints__
#define __chop_constraints__

#include <math.h>

/* Reminder for defines in math.h
trs 0    "Scale Rot Trans"
trs 1    "Scale Trans Rot"
trs 2    "Rot Scale Trans"
trs 3    "Rot Trans Scale"
trs 4    "Trans Scale Rot"
trs 5    "Trans Rot Scale"

xyz 0    "Rx Ry Rz"
xyz 1    "Rx Rz Ry"
xyz 2    "Ry Rx Rz"
xyz 3    "Ry Rz Rx"
xyz 4    "Rz Rx Ry"
xyz 5    "Rz Ry Rx"
*/


// Build a suffix for parameter/channel named lookups
// based on current transform index
#define DECLARE_C_SUFFIX \
string C_suffix = ""; \
if( C!=0 ) \
{ \
    C_suffix = itoa(C/9); \
}

// Write TRS to the Buffered Output
void
chwritebufTRS( const vector t; const  vector r; const vector s )
{
    chresizebuf(9);
    for( int i=0; i<3; ++i )
    {
        chwritebuf(i + 0, getcomp(t,i) );
        chwritebuf(i + 3, getcomp(r,i) );
        chwritebuf(i + 6, getcomp(s,i) );
    }
    chsetattr("", "__slerp__", 3, -1, 1);
}

// Read TRS channels of input node i
void
chinputTRS( const int i; vector t; vector r; vector s )
{
    t.x = chinput( i, C+0, I );
    t.y = chinput( i, C+1, I );
    t.z = chinput( i, C+2, I );
    
    r.x = chinput( i, C+3, I );
    r.y = chinput( i, C+4, I );
    r.z = chinput( i, C+5, I );
    
    s.x = chinput( i, C+6, I );
    s.y = chinput( i, C+7, I );
    s.z = chinput( i, C+8, I );
}

// Convert a matrix to TRS given a pivot and trs/xyz orders.
void
cracktransformTRS( const matrix m;
		   vector t; vector r; vector s;
		   const vector pivot; const int trs; const int xyz)
{
    t = cracktransform( trs, xyz, /*T=*/0, pivot, m );
    r = cracktransform( trs, xyz, /*R=*/1, pivot, m );
    s = cracktransform( trs, xyz, /*S=*/2, pivot, m );
}

// Convert a matrix to TRS given a pivot and trs/xyz orders.
void
cracktransformTRS( const matrix m;
		   vector t; vector r; vector s;
		   const vector pivot_translate;
		   const vector pivot_rotate;
                   const int trs; const int xyz)
{
    vector shears = 0;
    cracktransform( trs, xyz, pivot_translate, pivot_rotate, m, t,r,s, shears );
}

// Write a Matrix to the TRS Buffered Output, given a pivot and trs/xyz orders.
void
chwritebufMatrix( const matrix m; const vector pivot; const int trs; const int xyz)
{
    vector t=0;
    vector r=0;
    vector s=1;
    cracktransformTRS( m, t,r,s, pivot, trs, xyz );
    chwritebufTRS( t,r,s );
}

// Write a Matrix to the TRS Buffered Output, ( pivot at 0,0,0 and trs and xyz set to 0 )
void
chwritebufMatrix( const matrix m )
{
    vector pivot = 0;
    int trs = 0;
    int xyz = 0;
    chwritebufMatrix( m, pivot, trs, xyz );
}


// Keep Chop tx, ty, tz, rx, ry, rz, sx, sy, sz channels
struct chopTRS
{
    vector t; // tx ty tz channels
    vector r; // rx ry rz channels
    vector s; // sx ry rz channels

    void fromIdentity()
    {
	t = 0;
	r = 0;
	s = 1;
    }

    void fromMatrix(const matrix m;
		   const vector pivot; const int trs; const int xyz)
    {
	cracktransformTRS( m, t, r, s, pivot, trs, xyz );
    }
    void fromMatrix(const matrix m;
		   const vector pivot; const vector pivotrot; const int trs; const int xyz)
    {
	cracktransformTRS( m, t, r, s, pivot, pivotrot, trs, xyz );
    }

    void fromMatrix(const matrix m)
    {
	vector pivot = 0;
	int trs = 0;
	int xyz = 0;
	cracktransformTRS( m, t, r, s, pivot, trs, xyz );
    }

    matrix toMatrix(const int trs; const int xyz; const vector p;)
    {
        return maketransform(trs, xyz, t, r, s, p);
    }

    matrix toMatrix()
    {
	int trs=0;
	int xyz=0;
	vector p=0;
        return maketransform(trs, xyz, t, r, s, p);
    }
}

struct chopConstraintContext
{
    vector t; // tx ty tz channels
    vector r; // rx ry rz channels
    vector s; // sx sy sz channels

    float inputs[];
    int connected_inputs[];
    int use_inputs;
    string prefix;

    void init()
    {
	use_inputs = 0;
	t = 0;
	r = 0;
	s = 1;
	prefix = "";
	resize(inputs,0);
	resize(connected_inputs,0);
    }

    void init( vector t0; vector r0; vector s0; const string prefix0 )
    {
	use_inputs = 1;
	t = t0;
	r = r0;
	s = s0;
	prefix = prefix0;
	resize(inputs,9*1);
	resize(connected_inputs,1);
	inputs[0] = t0.x;
	inputs[1] = t0.y;
	inputs[2] = t0.z;
	inputs[3] = r0.x;
	inputs[4] = r0.y;
	inputs[5] = r0.z;
	inputs[6] = s0.x;
	inputs[7] = s0.y;
	inputs[8] = s0.z;
    }

    void fromMatrix(const matrix m;
		   const vector pivot; const int trs; const int xyz)
    {
	cracktransformTRS( m, t, r, s, pivot, trs, xyz );
    }

    void fromMatrix(const matrix m;
		   const vector pivot; const vector pivotrot;
                   const int trs; const int xyz)
    {
	cracktransformTRS( m, t, r, s, pivot, pivotrot, trs, xyz );
    }

    void fromMatrix(const matrix m)
    {
	vector pivot = 0;
	int trs = 0;
	int xyz = 0;
	cracktransformTRS( m, t, r, s, pivot, trs, xyz );
    }

    chopTRS fetchInput( const int i )
    {
	if( use_inputs )
	{
	    int i9 = 9;
	    i9*=i;

	    chopTRS c;
	    c.t.x = inputs[i9+0];
	    c.t.y = inputs[i9+1];
	    c.t.z = inputs[i9+2];
	    c.r.x = inputs[i9+3];
	    c.r.y = inputs[i9+4];
	    c.r.z = inputs[i9+5];
	    c.s.x = inputs[i9+6];
	    c.s.y = inputs[i9+7];
	    c.s.z = inputs[i9+8];
	    return c;
	}
	else
	{
	    chopTRS c;
	    chinputTRS( i, c.t, c.r, c.s );
	    return c;
	}
    }

    float fetchInput( const int i; const string name; int result )
    {
	if( use_inputs )
	{
	    printf("Error chopConstraintContext::fetchInput not implemented\n");
	    result = 0;
	    return 0;
	}
	else
	{
	    result = chindex( i, name ) >=0;
	    return result ? chinput(i, name, I) : 0.0;
	}
    }

    float fetchInput( const int i; const int index; int result )
    {
	if( use_inputs )
	{
	    printf("Error chopConstraintContext::fetchInput not implemented\n");
	    result = 0;
	    return 0;
	}
	else
	{
	    result = index>=0 && index < chnumchan(i);
	    return result ? chinput(i, index, I) : 0.0;
	}
    }

    matrix fetchInputMatrix( const int i )
    {
        chopTRS c;
        c = this->fetchInput(i);

	vector p = 0;
	int trs = 0;
	int xyz = 0;
        return maketransform(trs, xyz, c.t, c.r, c.s, p);
    }

    int isConnected( const int i )
    {
        int ret = 0;
	if( use_inputs )
            ret = i<len( connected_inputs) && connected_inputs[i]>=0;
	else
            ret = isconnected(i);

        return ret;
    }

    int numInputs()
    {
	if( use_inputs )
	    return len(connected_inputs);
	else
	    return ninputs();
    }

    // Evaluate a channel input by name
    float chinput( int i; const string name; int ret )
    {
	if( use_inputs )
	{
	    ret = 0;
	    return 0.0;
	}
	else
	{
	    int j =-1;
	    j = chindex( i, name );

	    if( j != -1 )
	    {
		ret = j!=-1;
		return chinput(i, j, I);
	    }
	    else
	    {
		ret = 0;
		return 0.0;
	    }
	}
    }

    // Evaluate a float parameter on the current CHOP node
    float chf( const string parm )
    {
	if( use_inputs )
	    return chf( prefix + parm );
	else
	    return chf(parm);
    }

    // Evaluate an integer parameter on the current CHOP node
    int chi( const string parm )
    {
	if( use_inputs )
	    return chi( prefix + parm );
	else
	    return chi(parm);
    }

    // Evaluate a vector parameter on the current CHOP node
    vector chv( const string parm )
    {
	if( use_inputs )
	    return chv( prefix + parm );
	else
	    return chv(parm);
    }
}

#define FIX_ANGLES(A) \
if (A.r.x < -180.0f)	   A.r.x += 360.0f; \
else if (A.r.x >  180.0f)  A.r.x -= 360.0f; \
if (A.r.y < -180.0f)	   A.r.y += 360.0f; \
else if (A.r.y >  180.0f)  A.r.y -= 360.0f; \
if (A.r.z < -180.0f)	   A.r.z += 360.0f; \
else if (A.r.z >  180.0f)  A.r.z -= 360.0f

#define BEGIN_CONSTRAINT(c) \
	chopConstraintContext c;   \
	c->init()

#define END_CONSTRAINT(c) \
	chwritebufTRS( c.t, c.r, c.s )

void
constraintbegin( chopConstraintContext c; const string obj_path; const int trs, xyz, mode)
{
    matrix m = 1.0;

    if( mode == 0 || mode ==3 )
	m = oppreconstrainttransform(obj_path);
    else if( mode == 1 )
	m = opparmtransform(obj_path);
    else if( mode == 2 )
	m = opparenttransform(obj_path);

    vector p=0;
    c->fromMatrix( m, p, trs, xyz );
}

void
constraintobject( chopConstraintContext c; const string obj_path, ref_path; const int trs, xyz; )
{
    matrix m=1;
    m = optransform( obj_path );

    if( ref_path!="" )
    {
	matrix rm=1;
	rm = optransform(ref_path);
	m *= invert(rm);
    }

    c->fromMatrix(m);
}

void
constraintobjectpretransform( chopConstraintContext c; const string obj_path; const int trs, xyz; )
{
    matrix m=1;
    m = oppretransform( obj_path );
    c->fromMatrix(m);
}

void
constrainttransform( chopConstraintContext c;
                         const int trs;
                         const int xyz;
                         const vector t;
                         const vector r;
                         const vector s;
                         const vector p;
                         const vector pr;
                         const int invert;
                         const int mode;
                         const int pmode;
                         )
{
    matrix m0=1;
    matrix m1=1;
    matrix m=1;

    // Fetch first input or use identity
    if( mode != 3 )
    {
        if( c->isConnected(0) )
            m0 = c->fetchInputMatrix(0);
    }

    // Fetch second input or use the parameters
    if( mode != 2 )
    {
        if( c->isConnected(1) )
        {
            m1 = c->fetchInputMatrix(1);
        }
        else
        {
            // Transform the pivot by the pivot matrix if it is connected
            vector pivot = p;
	    vector pivot_rot = pr;
	    vector shears = 0;
            if( c->isConnected(2) )
            {
		matrix m2m0i = c->fetchInputMatrix(2) * invert(m0);

		if( pmode==0 || pmode==2 )
		{
		    pivot = 0;
		    pivot *= m2m0i;
		}

		if (pmode == 1 || pmode==2)
		{
		    pivot_rot = 0;
		    vector a0 = {0,0,0};
		    vector ay = {0,1,0};
		    vector az = {0,0,-1};
		    a0 *= m2m0i;
		    ay *= m2m0i;
		    az *= m2m0i;
		    ay = ay - a0;
		    pivot_rot = lookat( a0, az, ay, xyz );
		}
            }

            m1 = maketransform(trs, xyz, t, r, s, pivot, pivot_rot, shears);
        }

        if( invert )
            m1 = invert( m1 );
    }

    // Post multiply
    if( mode == 0 )
        m = m1 * m0;
    // Pre multiply
    else if( mode == 1 )
        m = m0 * m1;
    // First Input
    else if( mode == 2 )
        m = m0;
    // Second Input
    else if( mode == 3 )
        m = m1;

    c->fromMatrix(m);
}

void
constraintpose( chopConstraintContext c;
                         const int trs;
                         const int xyz;
                         const vector t;
                         const vector r;
                         const vector s;
                         )
{
    matrix m=1;
    vector pivot = 0;
    m = maketransform(trs, xyz, t, r, s, pivot);
    c->fromMatrix(m);
}


void
constraintinvert( chopConstraintContext c )
{
    matrix m=1;

    // Fetch first input or use identity
    if( c->isConnected(0) )
    {
	m = c->fetchInputMatrix(0);
        m = invert( m );
    }

    c->fromMatrix(m);
}

void
constraintidentity( chopConstraintContext c )
{
    matrix m=1;
    c->fromMatrix(m);
}

void
constraintmultiply( chopConstraintContext c)
{
    matrix m=1;
    for( int i=0; i<100; i++ )
    {
	if( !c->isConnected(i) )
	    break;
	m = c->fetchInputMatrix(i) * m;
    }
    c->fromMatrix(m);
}

void
constraintlookat( chopConstraintContext c;
			 const int lookataxis;
			 const int lookupaxisx;
			 const int lookupaxisy;
			 const int lookupaxisz;
                         const vector lookat;
                         const int mode;
                         const vector uppos;
                         const vector upvec;
			 const float twist;
                         )
{

    matrix m0 = 1;
    vector p0 = 0;

    // Fetch first input or use identity
    if( c->isConnected(0) )
    {
        m0 = c->fetchInputMatrix(0);
	p0 *= m0;

	chopTRS c0;
	c0 = c->fetchInput(0);
	c.t = c0.t;
	c.r = c0.r;
	c.s = c0.s;
    }

    // Transform the lookat by second input if it is connected
    vector l = 0;
    if( c->isConnected(1) )
    {
	matrix m1=1;
	m1 = c->fetchInputMatrix(1);
	l *= m1;
    }
    else
    {
	l = lookat;
    }

    // Transform the lookat by second input if it is connected
    vector u = 0;
    if( c->isConnected(2) )
    {
	matrix m1=1;
	m1 = c->fetchInputMatrix(2);
	u *= m1;
	u -= p0;
    }
    else if( mode == 0 )
    {
	u = uppos-p0;
    }
    else
    {
	u = upvec;
    }
    u = normalize(u);

    matrix mlookat=1;
    mlookat = lookat(p0, l, u);

    {
	vector axis_z = {0,0,1};
	vector axis_y = {0,1,0};
	if( lookataxis == 0 ) // X-
	{
	    if( lookupaxisx == 0 ) // Y-
	    {
		axis_y = {0,-1,0};
		axis_z = {1,0,0};
	    }
	    else if( lookupaxisx == 1 ) // Z-
	    {
		axis_y = {-1,0,0};
		axis_z = {0,-1,0};
	    }
	    else if( lookupaxisx == 2 ) // Y+
	    {
		axis_y = {0,1,0};
		axis_z = {-1,0,0};
	    }
	    else if( lookupaxisx == 3 ) // Z+
	    {
		axis_y = {1,0,0};
		axis_z = {0,1,0};
	    }
	}
	else if( lookataxis == 1 ) // Y-
	{
	    if( lookupaxisy == 0 ) // X-
	    {
		axis_y = {0,0,1};
		axis_z = {-1,0,0};
	    }
	    else if( lookupaxisy == 1 ) // Z-
	    {
		axis_y = {0,0,1};
		axis_z = {0,-1,0};
	    }
	    else if( lookupaxisy == 2 ) // X+
	    {
		axis_y = {0,0,1};
		axis_z = {1,0,0};
	    }
	    else if( lookupaxisy == 3 ) // Z+
	    {
		axis_y = {0,0,1};
		axis_z = {0,1,0};
	    }
	}
	else if( lookataxis == 2 ) // Z-
	{
	    if( lookupaxisz == 0 ) // X-
	    {
		axis_y = {1,0,0};
		axis_z = {0,0,1};
	    }
	    else if( lookupaxisz == 1 ) // Y-
	    {
		axis_y = {0,-1,0};
		axis_z = {0,0,1};
	    }
	    else if( lookupaxisz == 2 ) // X+
	    {
		axis_y = {-1,0,0};
		axis_z = {0,0,1};
	    }
	    else // Y+
	    {
		axis_y = {0,1,0};
		axis_z = {0,0,1};
	    }

	}
	else if( lookataxis == 3 ) // X+
	{
	    if( lookupaxisx == 0 ) // Y-
	    {
		axis_y = {0,-1,0};
		axis_z = {-1,0,0};
	    }
	    else if( lookupaxisx == 1 ) // Z-
	    {
		axis_y = {1,0,0};
		axis_z = {0,-1,0};
	    }
	    else if( lookupaxisx == 2 ) // Y+
	    {
		axis_y = {0,1,0};
		axis_z = {1,0,0};
	    }
	    else if( lookupaxisx == 3 ) // Z+
	    {
		axis_y = {-1,0,0};
		axis_z = {0,1,0};
	    }
	}
	else if( lookataxis == 4 ) // Y+
	{
	    if( lookupaxisy == 0 ) // X-
	    {
		axis_y = {0,0,-1};
		axis_z = {1,0,0};
	    }
	    else if( lookupaxisy == 1 ) // Z-
	    {
		axis_y = {0,0,-1};
		axis_z = {0,-1,0};
	    }
	    else if( lookupaxisy == 2 ) // X+
	    {
		axis_y = {0,0,-1};
		axis_z = {-1,0,0};
	    }
	    else if( lookupaxisy == 3 ) // Z+
	    {
		axis_y = {0,0,-1};
		axis_z = {0,1,0};
	    }
	}
	else if( lookataxis == 5 ) // Z+
	{
	    if( lookupaxisz == 0 ) // X-
	    {
		axis_y = {-1,0,0};
		axis_z = {0,0,-1};
	    }
	    else if( lookupaxisz == 1 ) // Y-
	    {
		axis_y = {0,-1,0};
		axis_z = {0,0,-1};
	    }
	    else if( lookupaxisz == 2 ) // X+
	    {
		axis_y = {1,0,0};
		axis_z = {0,0,-1};
	    }
	    else if( lookupaxisz == 3 ) // Y+
	    {
		axis_y = {0,1,0};
		axis_z = {0,0,-1};
	    }
	}
	vector axis_xyz[];
	resize( axis_xyz, 3 );
	axis_xyz[0] = cross( axis_y, axis_z );
	axis_xyz[1] = axis_y;
	axis_xyz[2] = axis_z;

	matrix3 mat3 = 1;
	mat3 = set(axis_xyz);

	float twist1 = twist;
	if( c->isConnected(3) )
	    twist1 = chinput( 3, C/9, I );

	if( twist1 != 0 )
	{
	    matrix3 twistmat3 = 1;
	    rotate(twistmat3, radians(twist1), {0,0,1} );

	    mat3 = mat3 * twistmat3;
	}

	mlookat = matrix( mat3 ) * mlookat;
    }

    c.r = cracktransform(0,0,1,0,mlookat);
}

int
extractAllRot( const int mask; )
{
    return (mask&8) && (mask&16) && (mask&32);
}
void
extractMaskWeights( const int mask; float w; vector tw, rw, sw )
{
    tw = set( mask& 1?w:0.0, mask&  2?w:0.0, mask&  4?w:0.0 );
    rw = set( mask& 8?w:0.0, mask& 16?w:0.0, mask& 32?w:0.0 );
    sw = set( mask&64?w:0.0, mask&128?w:0.0, mask&256?w:0.0 );
}

int
safeDivide( vector v; const vector divider; )
{
    v.x = divider.x>0.0 ? v.x/divider.x : 0.0;
    v.y = divider.y>0.0 ? v.y/divider.y : 0.0;
    v.z = divider.z>0.0 ? v.z/divider.z : 0.0;

    return v.x!=0 || v.y!=0.0 || v.z!=0.0;
}

// Component-wise vector selection
vector vectorselect( vector a; vector b; vector c;)
{
    vector ret = 0;
    ret.x = c.x==0 ? a.x : b.x;
    ret.y = c.y==0 ? a.y : b.y;
    ret.z = c.z==0 ? a.z : b.z;
    return ret;
}

void
constraintblend( chopConstraintContext c;  const int method; const int rotblend; const int writemask; const int numblends; )
{
    DECLARE_C_SUFFIX

    int rorder = 0;

    int result = 0;
    int mode = 0;
    float input_blend = 0.0;

    mode = method;

    vector write_tw, write_rw, write_sw;
    int _rotblend = rotblend;
    int allrot = _rotblend==0 ? 0 : extractAllRot(writemask);
    extractMaskWeights( writemask, 1.0, write_tw, write_rw, write_sw );

    float weights[];
    int masks[];
    int allrots[];
    resize( weights, numblends );
    resize( masks, numblends );
    resize( allrots, numblends );


    // Difference mode
    if( mode==1 )
    {
	for( int i=1; i<numblends; i++ )
	{
	    weights[i] = c->fetchInput(i, "blend"+C_suffix, result);
	    if( result==0 )
		weights[i] = c->chf("blend"+itoa(i));

	    masks[i] = c->chi("mask"+itoa(i));
	}

	if( numblends>0 )
	{
	    chopTRS c0;
            if( c->isConnected(0) )
                c0 = c->fetchInput(0);
            else
                c0->fromIdentity();
	    c.t = c0.t;
	    c.r = c0.r;
	    c.s = c0.s;

	    // We need to fallback to Euler blending if we don't blend all the rotation channels
	    if( allrot )
	    {
		for( int i=1; i<numblends; i++ )
		{
		    if( extractAllRot(masks[i]) == 0 )
		    {
			allrot = 0;
			_rotblend = 0;
			break;
		    }
		}
	    }
	    if( _rotblend==1 && numblends<2 )
	    {
		_rotblend = 0;
	    }
	    else if( _rotblend==2 )
	    {
		if( numblends > 2)
		    _rotblend = 1;
		else if( numblends < 2)
		    _rotblend = 0;
	    }

	    if( _rotblend==0 )
	    {
                FIX_ANGLES(c0);
	    }

	    vector4 q = 0.0;
	    vector4 q0inv = 0;
	    vector4 qident = {0,0,0,1};

	    if( _rotblend==1 )
	    {
		q = eulertoquaternion(radians(c.r), rorder);
		q0inv = qinvert(q);
	    }

	    for( int i=1; i<numblends; i++ )
	    {
		vector tw, rw, sw;
                extractMaskWeights( masks[i], weights[i], tw, rw, sw );
		
		chopTRS ci;
		ci = c->fetchInput( i );
		
		c.t += (ci.t-c0.t)*tw*write_tw;

		if( _rotblend==1 )
		{
		    vector4 qi = eulertoquaternion(radians(ci.r), rorder);
		    qi = qmultiply( q0inv, qi );
		    qi = slerp(qident, qi, rw.x );
		    q = qmultiply( q, qi );
		}
		else if( _rotblend==2 )
		{
		    if( i==1 )
		    {
			// Convert euler to quaternion
			vector4 q0 = eulertoquaternion( radians(c.r), rorder );
			vector4 q1 = eulertoquaternion( radians(ci.r), rorder );

			// Interpolate quaternions
			q = slerp( q0, q1, rw.x );

			// Convert back to euler
			matrix qm = matrix( qconvert(q) );
			vector pivot = 0.0;
			c.r = cracktransform( 0, rorder, 1, pivot, qm );
		    }
		}
		else // Euler blending
		{
                    FIX_ANGLES(ci);
		    ci.r -= c0.r;
		    c.r += (ci.r)*rw*write_rw;
		}

		c.s += (ci.s-c0.s)*sw*write_sw;
	    }

	    if( _rotblend==1 )
	    {
		q = normalize(q);

		// Convert back to euler
		matrix3 qm = qconvert(q);
		c.r = cracktransform( 0, rorder, 1, 0, matrix(qm) );
	    }
	}
    }
    // Propertionnal mode
    else
    {

	// We need to keep per channel weight sums because of the per channel weight flags
	vector sum_tw = 0.0;
	vector sum_rw = 0.0;
	vector sum_sw = 0.0;
	int has_w = 0;

	for( int i=0; i<numblends; i++ )
	{
	    weights[i] = c->fetchInput(i, "blend"+C_suffix, result);
	    if( result==0 )
		weights[i] = c->chf("blend"+itoa(i));

	    masks[i] = c->chi("mask"+itoa(i));

	    vector tw, rw, sw;
	    extractMaskWeights( masks[i], weights[i], tw, rw, sw );

	    sum_tw += tw;
	    sum_rw += rw;
	    sum_sw += sw;

	    // Sum w is only a indicator that we have at least one non-zero weight
	    has_w = has_w ||
		    tw.x!=0.0 || tw.y!=0.0 || tw.z!=0.0 ||
	            rw.x!=0.0 || rw.y!=0.0 || rw.z!=0.0 ||
	            sw.x!=0.0 || sw.y!=0.0 || sw.z!=0.0;
	}


	if( has_w )
	{
	    vector zero = {0,0,0};

            chopTRS c0;
            if( c->isConnected(0) )
                c0 = c->fetchInput(0);
            else
                c0->fromIdentity();

            c.t = vectorselect( c0.t, zero, write_tw );
            c.r = vectorselect( c0.r, zero, write_rw );
            c.s = vectorselect( c0.s, zero, write_sw );

	    // We need to fallback to Euler blending if we don't blend all the rotation channels
	    if( allrot )
	    {
		for( int i=0; i<numblends; i++ )
		{
		    if( extractAllRot(masks[i]) == 0 )
		    {
			allrot = 0;
			_rotblend = 0;
			break;
		    }
		}
	    }
	    if( _rotblend==1 && numblends<2 )
	    {
		_rotblend = 0;
	    }
	    else if( _rotblend==2 )
	    {
		if( numblends > 2)
		    _rotblend = 1;
		else if( numblends < 2)
		    _rotblend = 0;
	    }

	    vector4 q = 0.0;
	    vector4 qbase;

	    for( int i=0; i<numblends; i++ )
	    {
		chopTRS ci;
		ci = c->fetchInput( i );

		vector tw, rw, sw;
		extractMaskWeights( masks[i],  weights[i], tw, rw, sw );

		tw *= write_tw;
		rw *= write_rw;
		sw *= write_sw;

		if( safeDivide(tw, sum_tw) )
		    c.t += ci.t*tw;

		if( _rotblend==1 )
		{
		    vector4 qi = eulertoquaternion(radians(ci.r), rorder);
		    if( i==0 )
			qbase = qi;
		    else if( dot(qbase, qi) < 0.0 )
			qi *= -1.0;

		    if( safeDivide(rw, sum_rw) )
		    {
			qi *= rw.x;
			q += qi;
		    }
		}
		else if( _rotblend==2 )
		{
		    if( i==0 )
		    {
			c.r = ci.r;
		    }
		    else if( i==1 )
		    {
			if( safeDivide(rw, sum_rw) )
			{
			    // Convert euler to quaternion
			    vector4 q0 = eulertoquaternion( radians(c.r), rorder );
			    vector4 q1 = eulertoquaternion( radians(ci.r), rorder );

			    // Interpolate quaternions
			    q = slerp( q0, q1, rw.x );

			    // Convert back to euler
			    matrix qm = matrix( qconvert(q) );
			    vector pivot = 0.0;
			    c.r = cracktransform( 0, rorder, 1, pivot, qm );
			}
		    }
		}
		else // Euler blending
		{
		    if( safeDivide(rw, sum_rw) )
		    {
                        FIX_ANGLES(ci);
			c.r += ci.r*rw;
		    }
		}

		if( safeDivide(sw, sum_sw) )
		    c.s += ci.s*sw;

	    }

	    if( _rotblend==1 )
	    {
		q = normalize(q);

		// Convert back to euler
		matrix3 qm = qconvert(q);
		c.r = cracktransform( 0, rorder, 1, 0, matrix(qm) );
	    }
	}
        else
        {
            chopTRS c0;
            if( c->isConnected(0) )
                c0 = c->fetchInput(0);
            else
                c0->fromIdentity();
            c.t = c0.t;
            c.r = c0.r;
            c.s = c0.s;
        }
    }

}


void
constraintsequence( chopConstraintContext c;  const float blend; const int rotblend; const int writemask; )
{
    int rorder = 0;
    chopTRS c0;
    chopTRS c1;

    int num = c->numInputs();

    vector write_tw, write_rw, write_sw;
    int _rotblend = rotblend;
    int allrot = _rotblend==0 ? 0 : extractAllRot(writemask);
    extractMaskWeights( writemask, 1.0, write_tw, write_rw, write_sw );

    int result = 0;
    float _blend = blend;

    if( num == 0 )
    {
        c0->fromIdentity();
        c.t = c0.t;
        c.r = c0.r;
        c.s = c0.s;
	return;
    }
    else if( num == 1 || _blend<=0.0 )
    {
        if( c->isConnected(0) )
            c0 = c->fetchInput(0);
        else
            c0->fromIdentity();

        c.t = c0.t;
        c.r = c0.r;
        c.s = c0.s;
	return;
    }
    else if( _blend>=(num-1) )
    {
        if( c->isConnected(0) )
            c0 = c->fetchInput(0);
        else
            c0->fromIdentity();

	if( c->isConnected(num-1) )
	    c1 = c->fetchInput(num-1);
        else
            c1->fromIdentity();

        c.t = vectorselect( c0.t, c1.t, write_tw );
        c.r = vectorselect( c0.r, c1.r, write_rw );
        c.s = vectorselect( c0.s, c1.s, write_sw );
	return;
    }

    int i0 = clamp( floor( _blend ), 0, num-1 );
    int i1 = clamp( i0+1, 0, num-1 );
    float f = frac( _blend );

    if( c->isConnected(i0) )
	c0 = c->fetchInput(i0);

    if( c->isConnected(i1) )
	c1 = c->fetchInput(i1);

    if( allrot==0 )
	_rotblend = 0;


    c.t = vectorselect( c.t, lerp( c0.t, c1.t, f ), write_tw );
    c.s = vectorselect( c.s, lerp( c0.s, c1.s, f ), write_sw );

    if( _rotblend==0 )
    {
        FIX_ANGLES(c0);
        FIX_ANGLES(c1);
	c.r = vectorselect( c.r, lerp( c0.r, c1.r, f ), write_rw );
    }
    else
    {
	// Convert euler to quaternion
	vector4 q0 = eulertoquaternion( radians(c0.r), rorder );
	vector4 q1 = eulertoquaternion( radians(c1.r), rorder );

	// Interpolate quaternions
	vector4 q = slerp( q0, q1, f );

	// Convert back to euler
	matrix qm = matrix( qconvert(q) );
	vector pivot = 0.0;
	c.r = cracktransform( 0, rorder, 1, pivot, qm );
    }
}

void
constraintsimpleblend( chopConstraintContext c;  const float blend; const int rotblend; const int writemask; )
{
    DECLARE_C_SUFFIX
    int result =0;

    int rorder = 0;
    chopTRS c0;
    chopTRS c1;

    vector write_tw, write_rw, write_sw;
    int _rotblend = rotblend;
    int allrot = _rotblend==0 ? 0 : extractAllRot(writemask);
    extractMaskWeights( writemask, 1.0, write_tw, write_rw, write_sw );

    float _blend = blend;
    _blend = c->fetchInput(1, "blend"+C_suffix, result);
    if( result==0 )
	_blend = blend;

    float f = clamp( _blend, 0, 1 );

    if( c->isConnected(0) )
	c0 = c->fetchInput(0);
    else
	c0->fromIdentity();

    if( c->isConnected(1) )
	c1 = c->fetchInput(1);
    else
	c1->fromIdentity();

    if( allrot==0 )
	_rotblend = 0;

    c.t = vectorselect( c0.t, lerp( c0.t, c1.t, f ), write_tw );
    c.s = vectorselect( c0.s, lerp( c0.s, c1.s, f ), write_sw );

    if( f<=0.0 )
    {
	c.r = c0.r;
    }
    else if( f>=1.0 )
    {
        c.r = vectorselect( c0.r, c1.r, write_rw );
    }
    else if( _rotblend==0 )
    {
        FIX_ANGLES(c0);
        FIX_ANGLES(c1);
	c.r = vectorselect( c0.r, lerp( c0.r, c1.r, f ), write_rw );
    }
    else
    {
	// Convert euler to quaternion
	vector4 q0 = eulertoquaternion( radians(c0.r), rorder );
	vector4 q1 = eulertoquaternion( radians(c1.r), rorder );

	// Interpolate quaternions
	vector4 q = slerp( q0, q1, f );

	// Convert back to euler
	matrix qm = matrix( qconvert(q) );
	vector pivot = 0.0;
	c.r = cracktransform( 0, rorder, 1, pivot, qm );
    }

}

void
constraintparent( chopConstraintContext c; )
{
    int rorder = 0;
    int trsorder = 0;
    vector pivot = 0;

    chopTRS c0; // First input, world space transform of the object to reparent.
    chopTRS c1; // Second input, world space tranform of the current parent.
    chopTRS c2; // Third Input, world space transform of the new parent.


    if( c->isConnected(0) )
        c0 = c->fetchInput(0);

    if( c->isConnected(1) )
        c1 = c->fetchInput(1);

    if( c->isConnected(2) )
        c2 = c->fetchInput(2);

    // Compute the Local Matrix by removing the original parent contribution.
    matrix mlocal = 1.0;
    {
        matrix m0 = 1;
        m0 = maketransform(trsorder, rorder, c0.t, c0.r, c0.s, pivot);

        matrix m1 = 1;
        m1 = maketransform(trsorder, rorder, c1.t, c1.r, c1.s, pivot);

        matrix m1inv = 0;
        m1inv = invert( m1 );

        mlocal = m0 * m1inv;
    }

    // Extract matrix for the new parent
    matrix m2 = 1;
    m2 = maketransform(trsorder, rorder, c2.t, c2.r, c2.s, pivot);

    // Apply the parent to the local matrix.
    matrix m;
    m = mlocal * m2;

    c->fromMatrix(m);
}

void
splitFloat(float val; int ival)
{
    val += 1; // Deal with numbers between (-1, 1)
    float fval = floor(val);
    ival = ((int)fval) - 1;
    val  = val - fval;
}

void
constraintpath( chopConstraintContext c;
		const int uparmtype;
		const float pos;
		const string soppath;
		const int lookatmode;
		const int lookupmode;
		const int lookataxis;
		const int lookupaxisx;
		const int lookupaxisy;
		const int lookupaxisz;
		const string dirattribute;
		const string upattribute;
		const vector upvector;
		const float roll;
		)
{

    string st = "op:" + soppath;
    string s = "time:" + itoa(I) +" op:" + soppath;
    int nprim = nprimitives(s);

    if( nprim == 0)
    {
	if( c->isConnected(0) )
	{
	    chopTRS c0;
	    c0 = c->fetchInput(0);
	    c.t = c0.t;
	    c.r = c0.r;
	    c.s = c0.s;
	}
	return;
    }

    matrix opm4 = optransform(st);
    matrix3 opm3 = matrix3(opm4);

    float arclengths[]; // Keep the arclength of all the curves in case we go past the last curve
    resize( arclengths, nprim );



    // Compute the primitive inde and adjust the coordinate
    //pos is a global position along the curves. posprim is a local position on the curve.

    float _roll = roll;
    float posprim = pos;
    if( c->isConnected(3) )
    {
	int result = 0;
	int i = C/9;
	float f;

	f = c->fetchInput(3, "pos"+itoa(i), result);
	if( result==1 )
	    posprim = f;

	f = c->fetchInput(3, "roll"+itoa(i), result);
	if( result==1 )
	    _roll = f;
    }

    int pnum = 0;

    // Normalized Len or Uniform mode
    if( uparmtype == 0 || uparmtype == 1 )
    {
	if (posprim != nprim)
	{
	    splitFloat(posprim, pnum);

	    // If we wrap the fractional part (which we do) we must
	    // also wrap the integer part.
	    pnum %= nprim;
	    if (pnum < 0)
		pnum += nprim;
	}
	else
	{
	    // We explicitly interpolate the last point of the last primitive.
	    // If it is wrapped, this is the same as wrapping.  If not, we
	    // actually can go all the way to the end of the path.
	    pnum = nprim - 1;
	    posprim = 1.0f;
	}
    }
    else if( uparmtype == 2 || uparmtype == 3 )
    {
	float accum_lengths[]; // Keep the arclength of all the curves in case we go past the last curve
	resize( accum_lengths, nprim );

	float lastlen = 0.0;
	float nextlen = 0.0;

	int found = 0;

	if( uparmtype == 3 )
	    posprim = -posprim;

	if( posprim > 0.0 )
	{
	    for( int i=0; i<nprim; i++ )
	    {
		arclengths[i] = primintrinsic(s, "arclength", i);
		nextlen = lastlen + arclengths[i];
		accum_lengths[i] = nextlen;

		if( posprim < nextlen )
		{
		    found = 1;
		    posprim -= lastlen;
		    pnum = i;
		    break;
		}
		lastlen = nextlen;
	    }
	    if( !found )
	    {
		float totallen = nextlen;

		posprim -= floor( posprim / totallen ) * totallen;

		lastlen = nextlen = 0.0;
		for( int i=0; i<nprim; i++ )
		{
		    nextlen = accum_lengths[i];

		    if( posprim < nextlen )
		    {
			posprim -= lastlen;
			pnum = i;
			break;
		    }

		    lastlen = nextlen;
		}
	    }
	}
	else
	{
	    float reversed_lengths[];
	    resize( reversed_lengths, nprim );

	    for( int i=0; i<nprim; i++ )
	    {
		arclengths[i] = primintrinsic(s, "arclength", i);
		nextlen = lastlen + arclengths[i];
		accum_lengths[i] = nextlen;
		lastlen = nextlen;
	    }
	    float totallen = nextlen;

	    lastlen = nextlen = 0.0;
	    for( int i=0; i<nprim; i++ )
	    {
		int j=nprim-i-1;
		nextlen = lastlen - arclengths[j];
		reversed_lengths[j] = nextlen;
		lastlen = nextlen;
	    }

	    if( posprim < -totallen )
		posprim -= ceil( posprim / totallen ) * totallen;

	    lastlen = nextlen = 0.0;
	    lastlen = 0.0;
	    for( int i=0; i<nprim; i++ )
	    {
		int j=nprim-i-1;
		nextlen = reversed_lengths[j];

		if( posprim > nextlen )
		{
		    posprim -= nextlen;
		    pnum = j;
		    break;
		}

		lastlen = nextlen;
	    }
	}

	// If we wrap the fractional part (which we do) we must
	// also wrap the integer part.
	pnum %= nprim;
	if (pnum < 0)
	    pnum += nprim;
    }


    vector2 u = 0;
    if( uparmtype == 0 ) // Convert from unit length to uniform
    {
	u = primuvconvert(s,posprim,pnum,PRIMUV_UNITLEN_TO_UNIT);

    }
    else if( uparmtype == 1 ) // Uniform
    {
	u = posprim;
    }
    else // Distance based
    {
	float total = arclengths[pnum];
	u = primuvconvert(s, posprim/total, pnum, PRIMUV_UNITLEN_TO_UNIT);
    }
    u.y = 0.0;

    if( c->isConnected(0) )
    {
	chopTRS c0;
	c0 = c->fetchInput(0);
	c.t = c0.t;
	c.r = c0.r;
	c.s = c0.s;
    }

    c.t = primuv(s,"P",pnum, u);
    c.t *= opm4;

    int connected1 = c->isConnected(1);
    if( connected1 || lookatmode == 1 || lookatmode == 2 )
    {
	vector lookatpos = 0;
	
	if( connected1 )
	    lookatpos *= c->fetchInputMatrix(1);
	else if( lookatmode == 1 )
	    lookatpos = c.t + normalize( primduv(s,pnum, u, 1,0) )*opm3;
	else if( lookatmode == 2 )
	    if( dirattribute!="" )
		lookatpos = c.t + normalize( primuv(s,dirattribute,pnum, u) )*opm3;

	vector lookupvec = 0;
	vector lookuppos = 0;

	int lmode = 1;

	if( c->isConnected(2) )
	{
	    lookuppos *= c->fetchInputMatrix(2);
	    lmode = 0;
	}
	else if( lookupmode == 0 )
	{
	    lookupvec = upvector;
	}
	else if( lookupmode == 1 )
	{
	    if( upattribute!="" )
	    {
		lookupvec = normalize(primuv(s,upattribute,pnum, u))*opm3;
	    }
	}
	else if( lookupmode == 2 )
	{
	    lookupvec = normalize( primduv(s,pnum, u, 1,0) )*opm3;
	}

	// Calling constraint lookat using a dummy context
	chopConstraintContext clookat;
	clookat->init( c.t, 0, 1, c.prefix );

	constraintlookat( clookat,
	    lookataxis, lookupaxisx, lookupaxisy, lookupaxisz,
	    lookatpos, lmode, lookuppos,lookupvec, _roll );

	c.r = clookat.r;
    }
}

int
searchsetup( chopConstraintContext c; const matrix opm4; vector searchpos; float searchdistance;)
{
    vector searchradius = {0,0,1};
    searchradius.z = searchdistance;
    searchpos = 0;
    searchdistance = 0.0;
    int doSearch = c->isConnected(3);
    if( doSearch )
    {
	matrix m3 = c->fetchInputMatrix(3) * invert(opm4);
	searchpos *= m3;
	searchradius *= m3;

	searchdistance = distance( searchpos, searchradius );
    }

    return doSearch;
}

void
searchpoints( chopConstraintContext c; const string s; const matrix opm4; const string group; const float dist; const int maxpnts; int ids[]; )
{
    vector searchpos;
    float searchdistance = dist;
    int doSearch = searchsetup(c, opm4, searchpos, searchdistance);

    // Use Points
    if( doSearch )
    {
	if( group == "*" || group == "" )
	    ids = nearpoints(s, searchpos, searchdistance, maxpnts);
	else
	    ids = nearpoints(s, group, searchpos, searchdistance, maxpnts);
    }
    else
    {
	ids = expandpointgroup( s, group );
    }
}

void
searchprims( chopConstraintContext c; const string s; const matrix opm4; const string group; const float dist; const int maxpnts; int ids[]; )
{
    vector searchpos;
    float searchdistance = dist;
    int doSearch = searchsetup(c, opm4, searchpos, searchdistance);

    // Use Primitives
    {
	int sorted_ids[];
	int primids[] = expandprimgroup( s, group );

	// Now we get the points from the primitives
	// and remove the duplicates
	foreach( int i; primids )
	    append( sorted_ids, primpoints(s,i) );

	// Sort, to make it easier to make the array unique
	sorted_ids = sort( sorted_ids );

	int len_sorted_ids = len( sorted_ids );
	int len_new = 0;

	// Grow the array to its maximum size
	resize( ids, len_sorted_ids );
	for( int i=0; i<len_sorted_ids; ++i )
	{
	    // Don't add if the previous is the same
	    if( i>0 && sorted_ids[i-1] == sorted_ids[i] )
		continue;
	    ids[ len_new++ ] = sorted_ids[i];
	}

	// Shrink to the real size.
	resize( ids, len_new );


	if( doSearch )
	{
	    if( group == "*" || group == "" )
		ids = nearpoints(s, searchpos, searchdistance, maxpnts);
	    else
	    {
		string pnt_group = "";
		for( int i=0; i<len_new; i++ )
		{
		    pnt_group += itoa( ids[i] );
		    pnt_group += " ";
		}
		ids = nearpoints(s, pnt_group, searchpos, searchdistance, maxpnts);
	    }
	}
    }
}
void
constraintsurface( chopConstraintContext c;
		const string soppath;
		const int subdi;
		const int mode;
		const string group;
		const string uvattribute;
		const vector uv;
		const string pattribute;
		const vector pv;
		const int lookatmode;
		const int lookupmode;
		const int lookataxis;
		const int lookupaxisx;
		const int lookupaxisy;
		const int lookupaxisz;
		const string dirattribute;
		const string upattribute;
		const vector upvector;
		const float roll;
		const float searchdist;
		const int   searchmaxpnt;
		)
{
    int MODE_UV = 0;
    int MODE_PRIMUV = 1;
    int MODE_POINTS = 2;
    int MODE_PRIMS = 3;
    int MODE_CLOSEST = 4;

    string st = "op:" + soppath;
    string s = "time:" + itoa(I) +" op:" + soppath;

    int nprim = nprimitives(s);

    if( nprim == 0)
    {
	if( c->isConnected(0) )
	{
	    chopTRS c0;
	    c0 = c->fetchInput(0);
	    c.t = c0.t;
	    c.r = c0.r;
	    c.s = c0.s;
	}
	return;
    }


    matrix opm4 = optransform(st);
    matrix3 opm3 = matrix3(opm4);

    vector primuv = 0;
    vector2 primuv2 = 0;
    int pnum = -1;
    int ids[];

    // UV Sample based, only works with unique 2D flat uvs.
    if( mode == MODE_UV )
    {
	vector _uv = 0;
	if( c->isConnected(3) )
	{
	    int ret = 0;
	    int ci = C/9 *3;

	    float f = 0;
	    f = c->fetchInput( 3, ci,ret );
	    if(ret)
		_uv.x = f;

	    f = c->fetchInput( 3, ci+1,ret );
	    if(ret)
		_uv.y = f;

	    f = c->fetchInput( 3, ci+2,ret );
	    if(ret)
		_uv.z = f;
	}
	else
	    _uv = uv;

	float d = uvdist(s, uvattribute, _uv, pnum, primuv);
    }
    else if( mode==MODE_POINTS || mode == MODE_POINTS )
    {
	// Use Points
	if( mode == MODE_POINTS )
	{
	    searchpoints(c, s, opm4, group, searchdist, searchmaxpnt, ids );
	}
	// Use Primitives
	else if( mode == MODE_PRIMS )
	{
	    searchprims(c, s, opm4, group, searchdist, searchmaxpnt, ids );
	}
    }
    // Use PrimUV directly
    else if( mode == MODE_PRIMUV )
    {
	if( c->isConnected(3) )
	{
	    primuv = 0;
	    pnum = -1;

	    int ret = 0;
	    int ci = C/9 *3;

	    float f = 0;
	    f = c->fetchInput( 3, ci,ret );
	    if(ret)
		pnum = int(f);

	    f = c->fetchInput( 3, ci+1,ret );
	    if(ret)
		primuv.x = f;

	    f = c->fetchInput( 3, ci+2,ret );
	    if(ret)
		primuv.y = f;
	}
	else
	{
	    primuv.x = uv.x;
	    primuv.y = uv.y;
	    primuv.z = 0;
	    pnum = (int)uv.z;
	}
    }
    else if( mode == MODE_CLOSEST )
    {
	vector _pv = 0;
	if( c->isConnected(3) )
	{
	    _pv *= c->fetchInputMatrix(3) * invert(opm4);
	}
	else
	{
	    _pv = pv;
	    _pv *= invert(opm4);
	}

	float d = uvdist(s, pattribute, _pv, pnum, primuv);
    }

    primuv2.x = primuv.x;
    primuv2.y = primuv.y;

    if( c->isConnected(0) )
    {
	chopTRS c0;
	c0 = c->fetchInput(0);
	c.t = c0.t;
	c.r = c0.r;
	c.s = c0.s;
    }

    // Points Mode finds the closest intersection on the geometry
    if( mode == MODE_POINTS || mode == MODE_PRIMS )
    {
	int numpoints = len(ids);
	c.t = 0;

	if( numpoints>0 )
	{
	    vector p=0;

	    foreach( int i; ids )
		p += point(s,"P", i);

	    p /= vector(numpoints);

	    xyzdist(s, p, pnum, primuv);
	}
    }

    c.t = primuv(s,"P",pnum, primuv);
    c.t *= opm4;

    int connected1 = c->isConnected(1);
    if( connected1 || lookatmode != 0 )
    {
	vector lookatpos = 0;
	
	if( connected1 )
	    lookatpos *= c->fetchInputMatrix(1);
	else
	{
	    if( lookatmode == 1 )
	    {
		if( dirattribute!="" )
		    lookatpos = c.t + normalize( primuv(s,dirattribute,pnum, primuv) )*opm3;
	    }
	    else
	    {
		if( lookatmode == 2 )
		    lookatpos = c.t + normalize( primduv(s,pnum, primuv2, 1,0) )*opm3;
		else if( lookatmode == 3 )
		    lookatpos = c.t + normalize( primduv(s,pnum, primuv2, 0,1) )*opm3;
	    }
	}

	vector lookupvec = 0;
	vector lookuppos = 0;

	int lmode = 1;

	if( c->isConnected(2) )
	{
	    lookuppos *= c->fetchInputMatrix(2);
	    lmode = 0;
	}
	else
	{
	    if( lookupmode == 0 )
	    {
		lookupvec = upvector;
	    }
	    else if( lookupmode == 1 )
	    {
		if( upattribute!="" )
		    lookupvec = normalize(primuv(s,upattribute,pnum, primuv))*opm3;
	    }
	    else
	    {
		if( lookupmode == 2 )
		    lookupvec = normalize( primduv(s,pnum, primuv2, 1,0) )*opm3;
		else if( lookupmode == 3 )
		    lookupvec = normalize( primduv(s,pnum, primuv2, 0,1) )*opm3;
	    }

	}

	// Calling constraint lookat using a dummy context
	chopConstraintContext clookat;
	clookat->init( c.t, 0, 1, c.prefix );

	constraintlookat( clookat,
	    lookataxis, lookupaxisx, lookupaxisy, lookupaxisz,
	    lookatpos, lmode, lookuppos,lookupvec, roll );

	c.r = clookat.r;
    }
}

void
constraintpoints( chopConstraintContext c;
		const string soppath;
		const int mode;
		const string group;
		const vector weights;
		const int lookatmode;
		const int lookupmode;
		const int lookataxis;
		const int lookupaxisx;
		const int lookupaxisy;
		const int lookupaxisz;
		const string dirattribute;
		const string upattribute;
		const vector upvector;
		const float roll;
		const float searchdist;
		const int   searchmaxpnt;
		)
{
    int MODE_POINTS = 0;
    int MODE_PRIMS = 1;

    string st = "op:" + soppath;
    string s = "time:" + itoa(I) +" op:" + soppath;

    int npoints = npoints(s);

    if( npoints == 0)
    {
	if( c->isConnected(0) )
	{
	    chopTRS c0;
	    c0 = c->fetchInput(0);
	    c.t = c0.t;
	    c.r = c0.r;
	    c.s = c0.s;
	}
	return;
    }

    matrix opm4 = optransform(st);
    matrix3 opm3 = matrix3(opm4);

    int ids[];

    // Use Points
    if( mode == MODE_POINTS )
    {
	searchpoints(c, s, opm4, group, searchdist, searchmaxpnt, ids );
    }
    // Use Primitives
    else if( mode == MODE_PRIMS )
    {
	searchprims(c, s, opm4, group, searchdist, searchmaxpnt, ids );
    }

    int numpoints = len(ids);

    vector points0 = 0;
    vector points1 = 0;
    vector points2 = 0;

    if( numpoints>=1 )
	points0 = point(s,"P", ids[0]);

    if( numpoints>=2 )
	points1 = point(s,"P", ids[1]);

    if( numpoints>=3 )
	points2 = point(s,"P", ids[2]);

    if( c->isConnected(0) )
    {
	chopTRS c0;
	c0 = c->fetchInput(0);
	c.t = c0.t;
	c.r = c0.r;
	c.s = c0.s;
    }

    c.t = 0;

    float total_w = 0.0;
    if( numpoints>0 )
    {
	int k = 0;
	foreach( int i; ids )
	{
	    float w =1.0;
	    if( k==0 )
		w = weights.x;
	    else if( k==1 )
		w = weights.y;
	    else if( k==2 )
		w = weights.z;

	    total_w += w;
	    c.t += point(s,"P", i) * vector(w);
	    k++;
	}

	c.t /= vector(total_w);
	c.t *= opm4;
    }

    int connected1 = c->isConnected(1);
    if( connected1 || lookatmode != 0 )
    {
	vector lookatpos = 0;
	
	if( connected1 )
	    lookatpos *= c->fetchInputMatrix(1);
	else
	{
	    if( lookatmode == 1 )
	    {
		if( dirattribute!="" )
		{
		    vector v =0;
		    if( numpoints>0 )
		    {
			int k = 0;
			foreach( int i; ids )
			{
			    float w =1.0;
			    if( k==0 )
				w = weights.x;
			    else if( k==1 )
				w = weights.y;
			    else if( k==2 )
				w = weights.z;

			    v += point(s,dirattribute, i) * vector(w);
			    k++;
			}

			v /= vector(total_w);
			v *= opm3;
		    }
		    lookatpos = c.t + v;
		}
	    }
	    else if( lookatmode == 2 )
	    {
		vector v = normalize(points1-points0);
		lookatpos = c.t + v*opm3;
	    }
	    else if( lookatmode == 3 )
	    {
		vector v10 = normalize(points1-points0);
		vector v20 = normalize(points2-points0);
		lookatpos = c.t + cross( v10, v20 )*opm3;
	    }
	}

	vector lookupvec = 0;
	vector lookuppos = 0;

	int lmode = 1;

	if( c->isConnected(2) )
	{
	    lookuppos *= c->fetchInputMatrix(2);
	    lmode = 0;
	}
	else
	{
	    if( lookupmode == 0 )
	    {
		lookupvec = upvector;
	    }
	    else if( lookupmode == 1 )
	    {
		if( upattribute!="" )
		{
		    vector v =0;
		    if( numpoints>0 )
		    {
			int k = 0;
			foreach( int i; ids )
			{
			    float w =1.0;
			    if( k==0 )
				w = weights.x;
			    else if( k==1 )
				w = weights.y;
			    else if( k==2 )
				w = weights.z;

			    v += point(s,upattribute, i);
			    k++;
			}

			v /= vector(total_w);
			v *= opm3;
		    }
		    lookupvec = v;
		}
	    }
	    else if( lookupmode == 2 )
	    {
		vector v = normalize(points1-points0);
		lookupvec = v*opm3;
	    }
	    else if( lookupmode == 3 )
	    {
		vector v10 = normalize(points1-points0);
		vector v20 = normalize(points2-points0);
		lookupvec = cross( v10, v20 )*opm3;
	    }
	}

	// Calling constraint lookat using a dummy context
	chopConstraintContext clookat;
	clookat->init( c.t, 0, 1, c.prefix );

	constraintlookat( clookat,
	    lookataxis, lookupaxisx, lookupaxisy, lookupaxisz,
	    lookatpos, lmode, lookuppos,lookupvec, roll );

	c.r = clookat.r;
    }
}

#endif //__chop_constraints__


#ifdef CONSTRAINTOBJECTVEX
chop
constraintobjectvex( int __iterate_over_channels__=9;
		     int __remove_time_dependent__=1;
		     string obj_path="";
		     string ref_path="";
		     int trs=0;
		     int xyz=0; )
{
    BEGIN_CONSTRAINT(c);
    constraintobject( c, obj_path, ref_path, trs, xyz );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTOBJECTPRETRANSFORMVEX
chop
constraintobjectpretransformvex( int __iterate_over_channels__=9;
		     int __remove_time_dependent__=1;
		     string obj_path="";
		     int trs=0;
		     int xyz=0; )
{
    BEGIN_CONSTRAINT(c);
    constraintobjectpretransform( c, obj_path, trs, xyz );
    END_CONSTRAINT(c);
}
#endif


#ifdef CONSTRAINTTRANSFORMVEX
chop
constrainttransformvex( int __iterate_over_channels__=9;
		     int trs=0;
		     int xyz=0;
                     vector t=0;
                     vector r=0;
                     vector s=1;
                     vector p=0;
                     vector pr=0;
                     const int invert=0;
                     const int mode=0;
                     const int pmode=0;
                     )
{
    BEGIN_CONSTRAINT(c);
    constrainttransform( c, trs, xyz, c->chv("t"), c->chv("r"), c->chv("s"), c->chv("p"), c->chv("pr"), invert, mode, pmode );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTPOSEVEX
chop
constraintposevex( int __iterate_over_channels__=9;
		     int update=0; // button
		     int clear=0;  // button
		     int mask=0;   // data, but used only by the button callbacks
		     int trs=0;
		     int xyz=0;
                     vector t=0;
                     vector r=0;
                     vector s=1;
                     )
{
    BEGIN_CONSTRAINT(c);
    constraintpose( c, trs, xyz, c->chv("t"), c->chv("r"), c->chv("s") );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTIDENTITYVEX
chop
constraintidentityvex( int __iterate_over_channels__=9;)
{
    BEGIN_CONSTRAINT(c);
    constraintidentity( c );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTINVERTVEX
chop
constraintinvertvex( int __iterate_over_channels__=9;)
{
    BEGIN_CONSTRAINT(c);
    constraintinvert( c );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTMULTIPLYVEX
chop
constraintnultiplyvex( int __iterate_over_channels__=9;)
{
    BEGIN_CONSTRAINT(c);
    constraintmultiply( c );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTBLENDVEX
chop
constraintblendvex( int __iterate_over_channels__=9;
		    int method=0;
		    int rotblend=0;
                    int writemask=511;
		    int numblends=0; )
{
    BEGIN_CONSTRAINT(c);
    constraintblend( c, method, rotblend, writemask, numblends );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTSEQUENCEVEX
chop
constraintsequencevex( int __iterate_over_channels__=9;
		    float blend=0.0;
		    int rotblend=0;
                    int writemask=511;)
{
    BEGIN_CONSTRAINT(c);
    constraintsequence( c, c->chf("blend"), rotblend, writemask );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTSIMPLEBLENDVEX
chop
constraintsimpleblendvex( int __iterate_over_channels__=9;
		    float blend=0.0;
		    int rotblend=0;
                    int writemask=511;)
{
    BEGIN_CONSTRAINT(c);
    constraintsimpleblend( c, c->chf("blend"), rotblend, writemask );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTPARENTVEX
chop
constraintparentvex( int __iterate_over_channels__=9; )
{
    BEGIN_CONSTRAINT(c);
    constraintparent( c );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTBEGINVEX
chop
constraintbeginvex( int __iterate_over_channels__=9;
		    string obj_path="";
		    int mode=0;
		    int trs=0;
		    int xyz=0;
		    )
{
    BEGIN_CONSTRAINT(c);
    constraintbegin( c, obj_path, trs, xyz, mode);
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTGETWORLDSPACEVEX
chop
constraintgetworldspacevex( int __iterate_over_channels__=9;
		    string obj_path="";
		    int trs=0;
		    int xyz=0;
		    )
{
    BEGIN_CONSTRAINT(c);
    vector p=0;
    c->fromMatrix( oppreconstrainttransform(obj_path), p, trs, xyz );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTGETLOCALSPACEVEX
chop
constraintgetlocalspacevex( int __iterate_over_channels__=9;
		    string obj_path="";
		    int trs=0;
		    int xyz=0;
		    )
{
    BEGIN_CONSTRAINT(c);
    vector p=0;
    c->fromMatrix( opparmtransform(obj_path), p, trs, xyz );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTGETPARENTSPACEVEX
chop
constraintgetparentspacevex( int __iterate_over_channels__=9;
		    string obj_path="";
		    int trs=0;
		    int xyz=0;
		    int parent_bone=0;
		    )
{
    BEGIN_CONSTRAINT(c);
    vector p=0;
    if( parent_bone==1 )
	c->fromMatrix( opparenttransform(obj_path), p, trs, xyz );
    else
	c->fromMatrix( opparentbonetransform(obj_path), p, trs, xyz );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTLOOKATVEX
chop
constraintlookatvex( int __iterate_over_channels__=9;
		     int lookataxis=0;
		     int lookupaxisx=0;
		     int lookupaxisy=0;
		     int lookupaxisz=0;
                     vector lookat=0;
		     int mode=0;
                     vector uppos=0;
                     vector upvec=1;
		     float twist=0;
                     )
{
    BEGIN_CONSTRAINT(c);
    constraintlookat( c, lookataxis, lookupaxisx, lookupaxisy, lookupaxisz, c->chv("lookat"), mode, c->chv("uppos"), c->chv("upvec"), c->chf("twist") );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTPATHVEX
chop
constraintpathvex( int __iterate_over_channels__=9;
		   int uparmtype=0;
                   float pos=0;
		   string soppath="";
		   const int lookatmode=0;
		   const int lookupmode=0;
		   const int lookataxis=0;
		   const int lookupaxisx=0;
		   const int lookupaxisy=0;
		   const int lookupaxisz=0;
		   const string dirattribute="";
		   const string upattribute="";
		   const vector upvector=0;
		   const float roll=0;
                   )
{
    BEGIN_CONSTRAINT(c);
    constraintpath( c, uparmtype, c->chf("pos"), soppath , lookatmode, lookupmode, lookataxis, lookupaxisx, lookupaxisy, lookupaxisz, dirattribute, upattribute, c->chv("upvector"), c->chf("roll") );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTSURFACEVEX
chop
constraintsurfacevex( int __iterate_over_channels__=9;
		   const string soppath="";
		   const int subdi=0;
		   const int mode=0;
		   const string group="";
		   const string uvattribute="";
                   const vector uv=0;
		   const string pattribute="";
                   const vector p=0;
		   const int lookatmode=0;
		   const int lookupmode=0;
		   const int lookataxis=0;
		   const int lookupaxisx=0;
		   const int lookupaxisy=0;
		   const int lookupaxisz=0;
		   const string dirattribute="";
		   const string upattribute="";
		   const vector upvector=0;
		   const float roll=0;
		   const float searchdist=1.0;
		   const int searchmaxpnt=1;
                   )
{
    BEGIN_CONSTRAINT(c);

    constraintsurface( c, soppath, subdi, mode, group, uvattribute, c->chv("uv"), pattribute, c->chv("p"), lookatmode, lookupmode, lookataxis, lookupaxisx, lookupaxisy, lookupaxisz, dirattribute, upattribute, c->chv("upvector"), c->chf("roll"), c->chf("searchdist"), c->chi("searchmaxpnt")
 );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTPOINTSVEX
chop
constraintpointsvex( int __iterate_over_channels__=9;
		   const string soppath="";
		   const int mode=0;
		   const string group="";
                   const vector weights=0;
		   const int lookatmode=0;
		   const int lookupmode=0;
		   const int lookataxis=0;
		   const int lookupaxisx=0;
		   const int lookupaxisy=0;
		   const int lookupaxisz=0;
		   const string dirattribute="";
		   const string upattribute="";
		   const vector upvector=0;
		   const float roll=0;
		   const float searchdist=1.0;
		   const int searchmaxpnt=1;
                   )
{
    BEGIN_CONSTRAINT(c);

    constraintpoints( c, soppath, mode, group, c->chv("weights"), lookatmode, lookupmode, lookataxis, lookupaxisx, lookupaxisy, lookupaxisz, dirattribute, upattribute, c->chv("upvector"), c->chf("roll"), c->chf("searchdist"), c->chi("searchmaxpnt")
 );
    END_CONSTRAINT(c);
}
#endif

#ifdef CONSTRAINTEXPORTVEX
chop
constraintexportvex( int __iterate_over_channels__=2;
		   const int constraints_on=0;
		   const string constraints_path="";
                   )
{
    chresizebuf(2);
    chwritebuf(0, constraints_on == 1 ? 1.0 : 0.0 );
    chwritebuf(1, constraints_path != "" ? 1.0 : 0.0 );
    chsetattr("","export", 1, -1, opfullpath(constraints_path) );
}
#endif

#ifdef CONSTRAINTCOMPOUNDEDVEX
// Example code of an expanded graph
chop
constraintexpandedvex( int __iterate_over_channels__=9 )
{
    BEGIN_CONSTRAINT(__start__);
    constraintbegin(__start__, /*trs=*/0, /*xyz=*/0, /*mode=*/0);

    BEGIN_CONSTRAINT(__object1__);
    constraintobject(__object1__, /*obj_path=*/"/obj/box_object1", /*ref_path=*/"", /*trs=*/0, /*xyz=*/0);

    BEGIN_CONSTRAINT(__object2__);
    constraintobject(__object2__, /*obj_path=*/"/obj/box_object2", /*ref_path=*/"", /*trs=*/0, /*xyz=*/0);

    BEGIN_CONSTRAINT(__blend1__);
    __blend1__.use_inputs = 1;
    resize( __blend1__.inputs, 3*3 )
    __blend1__.inputs[0] = __start__.t;
    __blend1__.inputs[1] = __start__.r;
    __blend1__.inputs[2] = __start__.s;
    __blend1__.inputs[3] = __object1__.t;
    __blend1__.inputs[4] = __object1__.r;
    __blend1__.inputs[5] = __object1__.s;
    __blend1__.inputs[6] = __object2__.t;
    __blend1__.inputs[7] = __object2__.r;
    __blend1__.inputs[8] = __object2__.s;

    constraintblend( __blend1__, 3 );
    END_CONSTRAINT(__blend1__);
}
#endif


