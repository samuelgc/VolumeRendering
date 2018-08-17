#ifndef __math_h__
#define __math_h__

/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *	Side Effects Software Inc
 *	477 Richmond Street West
 *	Toronto, Ontario
 *	Canada   M5V 3E7
 *	416-504-9876
 *
 * NAME:	math.h (VEX)
 *
 * COMMENTS:	This include file contains useful math constants
 */

#define M_E		2.7182818
#define M_LN10		2.3025850
#define M_LN2		0.6931471
#define M_LOG10E	0.4342944
#define M_LOG2E		1.4426950
#define M_PI		3.1415926
#define M_TWO_PI    6.2831852
#define M_PI_2		1.5707963
#define M_PI_4		0.7853981
#define M_SQRT1_2	0.7071067
#define M_SQRT2		1.4142135
#define M_TOLERANCE	0.0001

#define M_2SQRT6_3  1.6329931618554518  // 2 * sqrt(6) / 3
#define M_SQRT3     1.7320508075688772  // sqrt(3)
#define M_1_SQRT3   0.5773502691896257  // 1 / sqrt(3)
#define M_SQRT_2_3  0.816496580927726   // sqrt(2 / 3)

#if !environment(VEX_STRICT_COMPILE) && !defined(VEX_STRICT_COMPILE)
    // If not strict compiling, then allow for simpler naming conventions for
    // the math constants.
    // However, we do not allow for the simpler form of E since this conflicts
    // with global variables in some contexts.

    #define LN10	M_LN10
    #define LN2		M_LN2
    #define LOG10E	M_LOG10E
    #define LOG2E	M_LOG2E
    #define PI		M_PI
    #define PI_2	M_PI_2
    #define PI_4	M_PI_4
    #define SQRT1_2	M_SQRT1_2
    #define SQRT2	M_SQRT2
    #define TOLERANCE	M_TOLERANCE
#endif

// Defines for the maketransform() VEX function.  The function takes two
// integer arguments which determine the order that the transform gets made.
#define XFORM_SRT	0	// Scale, Rotate, Translate
#define XFORM_STR	1	// Scale, Translate, Rotate
#define XFORM_RST	2	// Rotate, Scale, Translate
#define XFORM_RTS	3	// Rotate, Translate, Scale
#define XFORM_TSR	4	// Translate, Scale, Rotate
#define XFORM_TRS	5	// Translate, Rotate, Scale

#define XFORM_XYZ	0	// Rotate order X, Y, Z
#define XFORM_XZY	1	// Rotate order X, Z, Y
#define XFORM_YXZ	2	// Rotate order Y, X, Z
#define XFORM_YZX	3	// Rotate order Y, Z, X
#define	XFORM_ZXY	4	// Rotate order Z, X, Y
#define XFORM_ZYX	5	// Rotate order Z, Y, X

// Defines for the primuvconvert() VEX function mode.
// REAL domain is based on the number of curve segments. (0 to nSegments)
//  A segment can hold multiple control points based on the curve degree.
// UNIT domain is the REAL domain normalized to fit in 0 to 1.
// UNITLEN domain maps the curve based on its length but normalized. (0..1)
// LEN domain maps the curve based on its length. (0..CurveLength)
#define PRIMUV_REAL_TO_UNIT	    0
#define PRIMUV_REAL_TO_UNITLEN	    1
#define PRIMUV_REAL_TO_LEN	    2
#define PRIMUV_UNIT_TO_REAL	    3
#define PRIMUV_UNIT_TO_UNITLEN	    4
#define PRIMUV_UNIT_TO_LEN	    5
#define PRIMUV_UNITLEN_TO_REAL	    6
#define PRIMUV_UNITLEN_TO_UNIT	    7
#define PRIMUV_UNITLEN_TO_LEN	    8
#define PRIMUV_LEN_TO_REAL	    9
#define PRIMUV_LEN_TO_UNIT	    10
#define PRIMUV_LEN_TO_UNITLEN	    11


#endif
