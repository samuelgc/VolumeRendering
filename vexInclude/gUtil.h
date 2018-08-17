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
 * NAME:	gUtil.h (Gallery Library, VEX)
 *
 * COMMENTS:	This include file contains many useful functions for writing
 *		shaders.
 */

#ifndef __gUtil_h__
#define __gUtil_h__

#ifndef G_MINFILTERSIZE
    #define G_MINFILTERSIZE	1.0e-6
#endif


#define CONST_COLOR(x,y,z)	{x,y,z}
#define CONST_VECTOR(x,y,z)	{x,y,z}
#define DYNAMIC_COLOR(x,y,z)	set(x,y,z)
#define DYNAMIC_VECTOR(x,y,z)	set(x,y,z)

float
gFilterWidthF(float x)
{
    //return max(abs(Du(x, "extrapolate", 1))+abs(Dv(x, "extrapolate", 1)), G_MINFILTERSIZE);
    return max(abs(Du(x))+abs(Dv(x)), G_MINFILTERSIZE);
}

float
gFilterWidthV(vector x)
{
    return max(sqrt(area(x)), G_MINFILTERSIZE);
}

int
gTileMap(float x, freq)
{
    return floor(x * freq);
}

#endif
