/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *	Side Effects Software Inc
 *	123 Front Street West, Suite 1401
 *	Toronto, Ontario
 *	Canada   M5J 2M2
 *	416-504-9876
 *
 * NAME:	pbrexports.h ( VEX )
 *
 * COMMENTS:	Macros to abstract the exports supported by pbr_pathtrace
 *		and the direct lighting VOP
 */

#ifndef __pbrexports_h
#define __pbrexports_h

vector[]
create_comp_array()
{
    vector	a[];
    resize(a, len(getcomponents()));
    for (int i = 0; i < len(a); i++)
	a[i] = {0,0,0};
    return a;
}

vector
sum_comp(const vector a[])
{
    vector sum = 0;
    for (int i = 0; i < len(getcomponents()); i++)
	sum += a[i];
    return sum;
}

vector[]
create_lcomp_array(int nlight)
{
    vector	a[];
    resize(a, nlight*len(getcomponents()));
    for (int i = 0; i < len(a); i++)
	a[i] = {0,0,0};
    return a;
}

#define LCOMP(array, lidx, cidx) array[lidx*len(getcomponents()) + cidx]

#define FOR_ALL_EXPORTS \
    foreach(int cidx; string cname; getcomponents())

#endif
