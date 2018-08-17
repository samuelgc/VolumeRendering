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
 * NAME:	comp.h ( VOP Library Functions, VEX )
 *
 * COMMENTS:	Support functions for VOPS.  Most of these library functions
 *		are shading context only.
 */

#ifndef __dualrest__
#define __dualrest__

struct DualRest
{
    vector      rest;
    float       rest_ratio;
    vector      rest2;
    float       rest2_ratio;
}

struct DualRest4
{
    vector      rest;
    float       rest_ratio;
    vector      rest2;
    float       rest2_ratio;
    float       time;
}

#endif
