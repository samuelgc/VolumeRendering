/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *      Neil Dickson
 *      Side Effects Software Inc
 *      1401-123 Front Street West
 *      Toronto, Ontario
 *      Canada   M5J 2M2
 *      416-504-9876
 *
 * NAME:        subdcurve.h (VOP Library Functions, VEX)
 *
 * COMMENTS:    Functions for evaluating subdivision curves
 */

#ifndef __subdcurve_h__
#define __subdcurve_h__

//******************************************************************************
// Functions for evaluating a single segment p0->p1 of a subd curve.
// A t-value of 0 corresponds with the beginning of the segment,
// i.e. (2/3)*p0 + (1/3)*(average of neighbours of p0).
// A t-value of 1 corresponds with the end of the segment,
// i.e. (2/3)*p1 + (1/3)*(average of neighbours of p1).
//
// subd_curve_evaluate          - for middle segment p0->p1 of a curve with
//                                no shared points
// subd_curve_evaluate_first    - for first segment p0->p1 of a curve with
//                                no shared points
// subd_curve_evaluate_last     - for last segment p0->p1 of a curve with
//                                no shared points
// subd_curve_evaluate_general  - for segment p0->(p0+diff) of a curve
//                                with shared points.  C0 is the average of
//                                neighbours of p0, minus p0.  C1 is the
//                                average of neighbours of p1, minus p1.
//
// subd_curve_derivative        - derivative of middle segment of a curve with
//                                no shared points
// subd_curve_derivative_first  - derivative of first segment of a curve with
//                                no shared points
// subd_curve_derivative_last   - derivative of last segment of a curve with
//                                no shared points
// subd_curve_derivative_general- derivative of segment p0->(p0+diff) of a curve
//                                with shared points.  C0 is the average of
//                                neighbours of p0, minus p0.  C1 is the
//                                average of neighbours of p1, minus p1.
//******************************************************************************

#define SUBD_EVAL(TYPE) \
TYPE subd_curve_evaluate_general(const float t; const TYPE C0; const TYPE p0; const TYPE diff; const TYPE C1) \
{ \
    float t3 = t*t*t; \
    float ti = 1-t; \
    float ti3 = ti*ti*ti; \
    /* Order of addition should reduce roundoff in common cases. */ \
    return p0 + (diff*t + (C0*((1.0/3.0)*ti3) + C1*((1.0/3.0)*t3))); \
} \
TYPE subd_curve_evaluate(const float t; const TYPE pn1; const TYPE p0; const TYPE p1; const TYPE p2) \
{ \
    TYPE diff = (p1-p0);            /* Vector from p0 to p1 */ \
    TYPE C0 = 0.5*(diff + (pn1-p0));/* Average of neighbours of p0, minus p0 */ \
    TYPE C1 = 0.5*((p2-p1) - diff); /* Average of neighbours of p1, minus p1 */ \
    return subd_curve_evaluate_general(t, C0, p0, diff, C1); \
} \
TYPE subd_curve_evaluate_first(const float t; const TYPE p0; const TYPE p1; const TYPE p2) \
{ \
    TYPE diff = (p1-p0);            /* Vector from p0 to p1 */ \
    TYPE C1 = 0.5*((p2-p1) - diff); /* Average of neighbours of p1, minus p1 */ \
    return subd_curve_evaluate_general(t, 0, p0, diff, C1); \
} \
TYPE subd_curve_evaluate_last(const float t; const TYPE pn1; const TYPE p0; const TYPE p1) \
{ \
    return subd_curve_evaluate_first(1-t, p1, p0, pn1); \
} \
TYPE subd_curve_derivative_general(const float t; const TYPE C0; const TYPE p0; const TYPE diff; const TYPE C1) \
{ \
    float t2 = t*t; \
    float ti = 1-t; \
    float ti2 = ti*ti; \
    /* Order of addition should reduce roundoff in common cases. */ \
    return diff + (C1*t2 - C0*ti2); \
} \
TYPE subd_curve_derivative(const float t; const TYPE pn1; const TYPE p0; const TYPE p1; const TYPE p2) \
{ \
    TYPE diff = (p1-p0);            /* Vector from p0 to p1 */ \
    TYPE C0 = 0.5*(diff + (pn1-p0));/* Average of neighbours of p0, minus p0 */ \
    TYPE C1 = 0.5*((p2-p1) - diff); /* Average of neighbours of p1, minus p1 */ \
    return subd_curve_derivative_general(t, C0, p0, diff, C1); \
} \
TYPE subd_curve_derivative_first(const float t; const TYPE p0; const TYPE p1; const TYPE p2) \
{ \
    TYPE diff = (p1-p0);            /* Vector from p0 to p1 */ \
    TYPE C1 = 0.5*((p2-p1) - diff); /* Average of neighbours of p1, minus p1 */ \
    return subd_curve_derivative_general(t, 0, p0, diff, C1); \
} \
TYPE subd_curve_derivative_last(const float t; const TYPE pn1; const TYPE p0; const TYPE p1) \
{ \
    return -subd_curve_derivative_first(1-t, p1, p0, pn1); \
} \
/*empty line at end of macro*/
SUBD_EVAL(float)
SUBD_EVAL(vector2)
SUBD_EVAL(vector)
SUBD_EVAL(vector4)
SUBD_EVAL(matrix2)
SUBD_EVAL(matrix3)
SUBD_EVAL(matrix)
#undef SUBD_EVAL

//******************************************************************************
// Functions for finding the bounding box of a subrange of a single segment
// of a subd curve.
// NOTE: Before calling, you *must* set boxmin to a large positive number
//       like 1e30 and boxmax to a large negative number like -1e30,
//       otherwise the bounding box may not be tight.  Reusing existing
//       values in this manner allows for easy accumulation of bounding boxes.
// To find the bounding box of the full segment, set t0 to 0 and t1 to 1.
// subd_curve_bbox          - for middle segment p0->p1 of a curve with
//                            no shared points
// subd_curve_bbox_first    - for first segment p0->p1 of a curve with
//                            no shared points
// subd_curve_bbox_last     - for last segment p0->p1 of a curve with
//                            no shared points
// subd_curve_bbox_general  - for segment p0->(p0+diff) of a curve
//                            with shared points.  C0 is the average of
//                            neighbours of p0, minus p0.  C1 is the
//                            average of neighbours of p1, minus p1.
//******************************************************************************

void subd_curve_bbox_between(const float t0; const float t1; const float C0; const float p0; const float diff; const float C1; float boxmin; float boxmax)
{
    float t02 = t0*t0;
    float ti0 = 1-t0;
    float ti02 = ti0*ti0;
    float t12 = t1*t1;
    float ti1 = 1-t1;
    float ti12 = ti1*ti1;

    /* Solve for t between t0 and t1 where derivative is zero: */
    /* diff + C0*t^2 - C1*(1-t)^2 == 0 */
    /* diff-C1 + 2*C1*t + (C0-C1)*t^2 == 0 */

    /* First, we do relatively fast checks to see if there are no solutions */
    /* between t0 and t1, since that's fairly common. */
    /* If derivatives have same sign at t0 and t1, */
    /* can only have 0 or 2 zero-crossings. */
    int dsignt0 = (diff + C0*t02 - C1*ti02 < 0);
    int dsignt1 = (diff + C0*t12 - C1*ti12 < 0);
    if (dsignt0 == dsignt1) {
        /* If 2nd derivative, 2*C0*t + 2*C1*(1-t) has no zero */
        /* between t0 and t1, there is no zero-crossing. */
        if ((C0*t0 + C1*ti0 < 0) == (C0*t1 + C1*ti1 < 0)) {
            return;
        }
        /* Find t where 2nd derivative is zero. */
        /* 2*C0*t + 2*C1*(1-t) == 0 */
        /* t == (C1-C0)/C1 == 1 - (C0/C1) */
        float ti2dz = (C0/C1);
        float t2dz = 1 - ti2dz ;
        /* If derivative at t2dz is same sign as at t0, no zero-crossing. */
        if ((diff + C0*t2dz*t2dz - C1*ti2dz*ti2dz < 0) == dsignt0) {
            return;
        }
        /* 2 zero-crossings between t0 and t1, unless zero-touch */
    }
    else {
        /* Exactly 1 zero-crossing between t0 and t1. */
    }
    float ta;
    float tb;
    int nzeros = solvequadratic(diff-C1, 2*C1, C0-C1, ta, tb);
    if (nzeros == 0) {
        return;
    }
    if (ta > t0 && ta < t1) {
        float va = subd_curve_evaluate_general(ta, C0, p0, diff, C1);
        boxmin = min(boxmin,va);
        boxmax = max(boxmax,va);
    }
    if (nzeros == 1) {
        return;
    }
    if (tb > t0 && tb < t1) {
        float vb = subd_curve_evaluate_general(tb, C0, p0, diff, C1);
        boxmin = min(boxmin,vb);
        boxmax = max(boxmax,vb);
    }
}
void subd_curve_bbox_general(const float t0; const float t1; const float C0; const float p0; const float diff; const float C1; float boxmin; float boxmax)
{
    /* Start with endpoints */
    float v0 = subd_curve_evaluate_general(t0, C0, p0, diff, C1);
    boxmin = min(boxmin,v0);
    boxmax = max(boxmax,v0);
    float v1 = subd_curve_evaluate_general(t1, C0, p0, diff, C1);
    boxmin = min(boxmin,v1);
    boxmax = max(boxmax,v1);

    subd_curve_bbox_between(t0, t1, C0, p0, diff, C1, boxmin, boxmax);
}
void subd_curve_bbox(const float t0; const float t1; const float pn1; const float p0; const float p1; const float p2; float boxmin; float boxmax)
{
    float diff = (p1-p0);            /* Vector from p0 to p1 */
    float C0 = 0.5*(diff + (pn1-p0));/* Average of neighbours of p0, minus p0 */
    float C1 = 0.5*((p2-p1) - diff); /* Average of neighbours of p1, minus p1 */
    subd_curve_bbox_general(t0, t1, C0, p0, diff, C1, boxmin, boxmax);
}
void subd_curve_bbox_first(const float t0; const float t1; const float p0; const float p1; const float p2; float boxmin; float boxmax)
{
    float diff = (p1-p0);            /* Vector from p0 to p1 */
    float C1 = 0.5*((p2-p1) - diff); /* Average of neighbours of p1, minus p1 */
    subd_curve_bbox_general(t0, t1, 0, p0, diff, C1, boxmin, boxmax);
}
void subd_curve_bbox_last(const float t0; const float t1; const float pn1; const float p0; const float p1; float boxmin; float boxmax)
{
    subd_curve_bbox_first(1-t1, 1-t0, p1, p0, pn1, boxmin, boxmax);
}

#define SUBD_BBOX(TYPE,VEC_SIZE) \
void subd_curve_bbox_general(const float t0; const float t1; const TYPE C0; const TYPE p0; const TYPE diff; const TYPE C1; TYPE boxmin; TYPE boxmax) \
{ \
    /* Start with endpoints */ \
    TYPE v0 = subd_curve_evaluate_general(t0, C0, p0, diff, C1); \
    boxmin = min(boxmin,v0); \
    boxmax = max(boxmax,v0); \
    TYPE v1 = subd_curve_evaluate_general(t1, C0, p0, diff, C1); \
    boxmin = min(boxmin,v1); \
    boxmax = max(boxmax,v1); \
    for (int i = 0; i < VEC_SIZE; ++i) { \
        subd_curve_bbox_between(t0, t1, C0[i], p0[i], diff[i], C1[i], boxmin[i], boxmax[i]); \
    } \
} \
void subd_curve_bbox(const float t0; const float t1; const TYPE pn1; const TYPE p0; const TYPE p1; const TYPE p2; TYPE boxmin; TYPE boxmax) \
{ \
    TYPE diff = (p1-p0);            /* Vector from p0 to p1 */ \
    TYPE C0 = 0.5*(diff + (pn1-p0));/* Average of neighbours of p0, minus p0 */ \
    TYPE C1 = 0.5*((p2-p1) - diff); /* Average of neighbours of p1, minus p1 */ \
    subd_curve_bbox_general(t0, t1, C0, p0, diff, C1, boxmin, boxmax); \
} \
void subd_curve_bbox_first(const float t0; const float t1; const TYPE p0; const TYPE p1; const TYPE p2; TYPE boxmin; TYPE boxmax) \
{ \
    TYPE diff = (p1-p0);            /* Vector from p0 to p1 */ \
    TYPE C1 = 0.5*((p2-p1) - diff); /* Average of neighbours of p1, minus p1 */ \
    subd_curve_bbox_general(t0, t1, (TYPE)0, p0, diff, C1, boxmin, boxmax); \
} \
void subd_curve_bbox_last(const float t0; const float t1; const TYPE pn1; const TYPE p0; const TYPE p1; TYPE boxmin; TYPE boxmax) \
{ \
    subd_curve_bbox_first(1-t1, 1-t0, p1, p0, pn1, boxmin, boxmax); \
} \
/*empty line at end of macro*/
SUBD_BBOX(vector2,2)
SUBD_BBOX(vector,3)
SUBD_BBOX(vector4,4)
#undef SUBD_BBOX

//******************************************************************************
// Function for taking points that you want a subd curve to pass through
// and finding the CVs that will accomplish that.
//******************************************************************************
#define SUBD_INTERPOLATE(TYPE) \
void subd_curve_interpolate(const TYPE ptstointerpolate[]; const int isclosedloop; TYPE subdcurvepts[]) \
{ \
    /* Start with ptstointerpolate as initial guess */ \
    subdcurvepts = ptstointerpolate; \
    \
    /* Be careful not to overwrite values before they're done being used. */ \
    \
    int n = len(ptstointerpolate); \
    int starti = !isclosedloop; \
    int endi = n - !isclosedloop; \
    \
    /* It should have converged pretty well after 10 iterations. */ \
    /* 20 is just in case it needs to be extra accurate. */ \
    for (int iteration = 0; iteration < 20; ++iteration) \
    { \
        TYPE prev = subdcurvepts[isclosedloop ? (n-1) : 0]; \
        TYPE curr = subdcurvepts[starti]; \
        TYPE first = curr; \
        for (int i = starti; i < endi; ++i) \
        { \
            TYPE next = (isclosedloop && i == endi-1) ? first : subdcurvepts[i+1]; \
            \
            /* This is just Jacobi's method for solving a diagonally   */ \
            /* dominant system of linear equations, applied to:        */ \
            /* (1/6)*prev + (2/3)*curr + (1/6)*next = target, namely:  */ \
            /* curr' = (3/2)*target - (1/4)*prev - (1/4)*next          */ \
            /* It seems to converge pretty fast, and it's very simple. */ \
            TYPE val = 1.5f*ptstointerpolate[i] - 0.25f*(prev + next); \
            \
            subdcurvepts[i] = val; \
            \
            prev = curr; \
            curr = next; \
        } \
    } \
} \
/*empty line at end of macro*/
SUBD_INTERPOLATE(float)
SUBD_INTERPOLATE(vector2)
SUBD_INTERPOLATE(vector)
SUBD_INTERPOLATE(vector4)
SUBD_INTERPOLATE(matrix2)
SUBD_INTERPOLATE(matrix3)
SUBD_INTERPOLATE(matrix)
#undef SUBD_INTERPOLATE

#endif
