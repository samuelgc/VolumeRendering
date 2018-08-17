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
 * NAME:	voptype.h ( VOP Library, C++)
 *
 * COMMENTS:
 *	This header file defines some "generic" types which can be used in
 *	either VEX or RSL code.  There is a corresponding voptype.h file
 *	included in ri_shaders (for compiling RSL code).
 */

#ifndef __voptype__
#define __voptype__

//
// Some type defines which can be used in code intended to be compiled for both
// VEX and RSL.
//
#define VOPvector4	vector4
#define VOPvector	vector
#define VOPfloat	float
#define VOPint		int
#define VOPstring	string
#define VOPmatrix3	matrix3
#define VOPmatrix4	matrix
#define VOPpoint	vector
#define VOPnormal	vector
#define VOPcolor	vector
#define VOPbsdf		bsdf

#define VOPuniform
#define VOPoutput

// How to declare a constant
#define VOP_C_COLOR(r,g,b)	{r,g,b}
#define VOP_C_VECTOR(x,y,z)	{x,y,z}
#define VOP_C_POINT(x,y,z)	{x,y,z}
#define VOP_C_NORMAL(x,y,z)	{x,y,z}

// How to declare a varying value
#define VOP_V_COLOR(r,g,b)	set(r,g,b)
#define VOP_V_VECTOR(x,y,z)	set(x,y,z)
#define VOP_V_POINT(x,y,z)	set(x,y,z)
#define VOP_V_NORMAL(x,y,z)	set(x,y,z)

// How to cast a variable
#define VOP_CAST_FLOAT(x)	((float)(x))

#define VOP_TRUE(condition)	(condition)
#define VOP_FALSE(condition)	(!condition)

#endif
