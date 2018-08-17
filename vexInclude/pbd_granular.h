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
 * NAME:	pbd_granular.h ( VOP Library Functions, VEX )
 *
 * COMMENTS:	Support functions for PBD.
 */
#ifndef __pbd_granular_h__
#define __pbd_granular_h__

#include <math.h>

// collide two spheres using DEM method
vector collideSpheres(vector posA, posB, velA, velB;
                      float radiusA, radiusB, spring, damping, shear, attraction)
{
    // calculate relative position
    vector relPos = posB - posA;

    float dist = length(relPos);
    float collideDist = radiusA + radiusB;

    vector force = 0;
    if (dist < collideDist)
    {
        vector norm = relPos / dist;

        // relative velocity
        vector relVel = velB - velA;

        // relative tangential velocity
        vector tanVel = relVel - (dot(relVel, norm) * norm);

        // spring force
        force = -spring*(collideDist - dist) * norm;
        // dashpot (damping) force
        force += damping*relVel;
        // tangential shear force
        force += shear*tanVel;
        // attraction
        force += attraction*relPos;
    }

    return force;
}

vector demForcesPC(string geo; vector pos, vel; int thispt; 
                   float radius, spring, damping, shear, attraction)
{
    int h = pcopen(geo, "P", pos, radius * 2, 1000, "preload", 1);
    int otherpt;
    vector otherP;
    vector otherv;
    vector force = 0;

    int npt = pcnumfound(h);
    for(int i=0; i < npt; ++i)
    {
        otherpt = pcimportbyidxi(h, "point.number", i);
        if (otherpt == thispt)
            continue;
        otherP = pcimportbyidxv(h, "P", i);
        otherv = pcimportbyidxv(h, "v", i);
        force += collideSpheres(pos, otherP, vel, otherv, radius, radius,
                                spring, damping, shear, attraction);
    }
    pcclose(h);
    return force;
}

vector demForcesNeighbors(string geo; vector pos, vel; int neighbors[];
                   float radius, spring, damping, shear, attraction)
{
    vector otherP;
    vector otherv;
    vector force = 0;

    foreach(int otherpt; neighbors)
    {
        otherP = point(geo, "P", otherpt);
        otherv = point(geo, "v", otherpt);

        force += collideSpheres(pos, otherP, vel, otherv, radius, radius,
                                spring, damping, shear, attraction);
    }
    return force;
}

vector []getPoints(int geo)
{
    vector curP[];
    int i, npts = npoints(geo);
    resize(curP, npts);
    for (i = 0; i < npts;  ++i)
        curP[i] = point(geo, "P", i);
    return curP;
}

void setPoints(int geoh; vector curP[])
{
    int npts = len(curP);
    for (int i = 0; i < npts;  ++i)
        setpointattrib(geoh, "P", i, curP[i]);
}

int[] getNeighbors(const int geo; const vector pos; const int thispt; 
                   const float radius; const int maxn)
{
    int h = pcopen(geo, "P", pos, radius, maxn);//, "preload", 1);
    int otherpt;
    int neighbors[] = {};
    int npt = pcnumfound(h);
    resize(neighbors, npt-1);
    int idx=0;
    for(int i=0; i < npt; ++i)
    {
        otherpt = pcimportbyidxi(h, "point.number", i);
        if (otherpt != thispt)
            neighbors[idx++] = otherpt;
    }
    pcclose(h);
    return neighbors;
}

// Fetches point number, position, and 1/mass for the given vertex.
// Looks up point position from curP if provided else gets the P attribute.
void getVertexConstraintAttrib(const int geo, prim, vtx; const vector curP[];
                               int pt; vector pos; float w)
{
    pt = vertexpoint(geo, vertexindex(geo, prim, vtx));
    w = 1.0 / point(geo, "mass", pt);
    // Lookup from curP array if provided.
    if (len(curP))
        pos = curP[pt];
    else
        pos = point(geo, "P", pt);
}

void polyDistanceConstraint(const int geo, geoh, prim;
                            const float weight, stiffness; const int niter;
                            vector curP[])
{
    int nvtx = primvertexcount(geo, prim);
    float kp = 1 - pow(1 - stiffness, 1.0 / niter);
    int pt1, pt2;
    float w1, w2;
    vector p1, p2;
    // Get vertex attributes.
    getVertexConstraintAttrib(geo, prim, 0, curP, pt1, p1, w1);
    getVertexConstraintAttrib(geo, prim, 1, curP, pt2, p2, w2);
    // Nothing to do for two infinite mass points.
    if (w1 == 0 && w2 == 0)
        return;

    float restdist = prim(geo, "restlength", prim);
    vector p1p2 = p1 - p2;
    float curdist = length(p1p2);

    // Constraint and gradient.
    float C = curdist - restdist;
    vector gradC = p1p2 / curdist;

    // Update vectors.
    vector4 dp1 = - weight * kp * w1 / (w1 + w2) * C * gradC;
    vector4 dp2 = + weight * kp * w2 / (w1 + w2) * C * gradC;
    dp1.w = dp2.w = weight;

    // Update curP if provided.
    if (len(curP))
    {
        curP[pt1] += vector(dp1);
        curP[pt2] += vector(dp2);
    }
    else
    {
        // Add deltas if finite mass particle.
        if (w1 > 0)
            setpointattrib(geoh, "dPw", pt1, dp1, "add");
        if (w2 > 0)
            setpointattrib(geoh, "dPw", pt2, dp2, "add");
    }
}

void tetVolumeConstraint(const int geo, geoh, prim;
                         const float weight, stiffness; const int niter;
                         vector curP[])
{
    if (primvertexcount(geo, prim) != 4)
        return;

    float kp = 1 - pow(1 - stiffness, 1.0 / niter);
    int pt1, pt2, pt3, pt4;
    float w1, w2, w3, w4;
    vector p1, p2, p3, p4;

    getVertexConstraintAttrib(geo, prim, 0, curP, pt1, p1, w1);
    getVertexConstraintAttrib(geo, prim, 1, curP, pt2, p2, w2);
    getVertexConstraintAttrib(geo, prim, 2, curP, pt3, p3, w3);
    getVertexConstraintAttrib(geo, prim, 3, curP, pt4, p4, w4);

    // Nothing to do for four infinite mass points.
    if (w1 == 0 && w2 == 0 && w3 == 0 && w4 == 0)
        return;

    float restvol = prim(geo, "restvolume", prim);
    float curvol = dot(cross(p1 - p2, p1 - p3), p1 - p4) / 6;

    // Constraint and gradients.
    float C = curvol - restvol;
    // Check for inverted tet.
    if (curvol < 0)
    {
        // TODO - turn this into constraint to un-invert.
        //printf("inverted tet!\n");
        return;
    }
    //printf("restvol = %g, curvol = %g, C = %g\n", restvol, curvol,C);
    vector gradC1 = cross(p3 - p2, p4 - p2) / 6;
    vector gradC2 = -cross(p4 - p3, p1 - p3) / 6;
    vector gradC3 = cross(p1 - p4, p2 - p4) / 6;
    vector gradC4 = -cross(p2 - p1, p3 - p1) / 6;

    // Scale factor.
    float s = weight * kp * C / (w1 * length2(gradC1) + 
                                 w2 * length2(gradC2) +
                                 w3 * length2(gradC3) +
                                 w4 * length2(gradC4));

    // Update vectors.
    vector4 dp1 = -w1 * s * gradC1;
    vector4 dp2 = -w2 * s * gradC2;
    vector4 dp3 = -w3 * s * gradC3;
    vector4 dp4 = -w4 * s * gradC4;
    dp1.w = dp2.w = dp3.w = dp4.w = weight;

    // Update curP if provided.
    if (len(curP))
    {
        curP[pt1] += vector(dp1);
        curP[pt2] += vector(dp2);
        curP[pt3] += vector(dp3);
        curP[pt4] += vector(dp4);        
    }
    else 
    {
        // Add deltas for all finite mass particles.
        if (w1 > 0)
            setpointattrib(geoh, "dPw", pt1, dp1, "add");
        if (w2 > 0)
            setpointattrib(geoh, "dPw", pt2, dp2, "add");
        if (w3 > 0)
            setpointattrib(geoh, "dPw", pt3, dp3, "add");
        if (w4 > 0)
            setpointattrib(geoh, "dPw", pt4, dp4, "add");
    }
}

float
__taylor_exp(const float is)
{
    // This becomes increasingly unstable as we move away from 0,
    // however we are supposed to be in the range -.3 to .3 here.
    // To account for people dialing shockscale to 11, we want to
    // clamp the incoming s.
    float s = clamp(is, -1.0, 1.0);
    float s2 = s * s;
    float s3 = s2 * s;
    // Taylor expansion of e^s
    return (1 + s + s2/2 + s3 / 6);
}

void
computePointDistanceDelta(vector4 dPr, dPa;
			  const vector pi, pj;
			  const float mass, massj;
			  const float curdist;
			  const float restdist;
			  const float kpr, kpa;
			  const float wr, wa;
			  const int shocktype;
			  const vector shockaxis;
			  const float shockscale)
{
    vector r = pi - pj;

    // Constraint and gradient.
    float C = curdist - restdist;
    vector gradC = r / curdist;

    float ks = 1;
    if (shocktype == 2)		// local
    {
	ks = __taylor_exp( shockscale * dot(r, shockaxis));
    }

    float weight = ks * massj / (mass + ks * massj); 

    // Handle opposing weight being 0, ie, infinite
    weight = select(massj == 0, 1, weight);

    // Use weighted distance constraint if within attract distance.
    vector4 dpj = -kpa * C * gradC * (weight * wa);
    dpj.w = wa;
    dPa += dpj;
    
    // Use weighted (inequality) repel constraint if within repel distance.
    dpj = -kpr * C * gradC * weight;
    dpj.w = wr;
    dpj = select(C < 0, dpj, 0);
    dPr += dpj;
}

vector4 computeFrictionDelta(const vector pi, pj;
			     const vector dpi, dpj;
			     const float mass, massj;
			     const float curdist;
			     const float restdist;
			     const float mus, muk;
			     const int shocktype;
			     const vector shockaxis;
			     const float shockscale)
{
    vector4		dP = 0;

    vector r = pj - pi;
    vector relvel = dpi - dpj;

    float ks = 1;
    if (shocktype == 2)		// local
    {
	// r has the opposite sign here than in computePoint..
	ks = __taylor_exp( - shockscale * dot(r, shockaxis) );
    }

    // Find the direction between the particles.
    // This is the normal of the impact.
    r /= curdist;

    // The idea is to approximate our normal force by
    // how much we interpenetrated.
    float nmlforce = restdist - curdist;

    // Find tangent component of our relative velocity
    vector dPtan = relvel - dot(relvel, r) * r;
    float ldPtan = length(dPtan);

    // Compute kinetic friction scalar
    // Clamped to brining to a full stop.
    float fkin = muk * nmlforce;
    fkin = select(fkin > ldPtan, 1, fkin / ldPtan);

    // Both particles receive friction, so we bias
    // the application according to mass.
    float weight = ks * massj / (mass + ks * massj);
    weight = select(massj == 0, 1, weight);

    // Friction coefficient.
    // If we are in the static cone, apply the full friction,
    // ie, attempt to bring the object to a stop.
    // Otherwise use the kinetic friction.
    float fcoeff = select(ldPtan < mus * nmlforce, 1, fkin);

    dP -= weight * dPtan * fcoeff;
    dP.w = 1;

    return dP;
}

// #define DISABLE_CONSTRAINT_AVERAGING

vector pointDistanceUpdate(const int geo; const vector pi; 
			    const int neighbors[];
                            const float pscale, mass;
			    const string massname;
			    const float sor;
			    const float kr, wr, ka, wa;
			    const float timestep; 
			    const int niter;
			    const vector pprevious;
			    const float frscale, mu_s, mu_k;
			    const int shocktype;
			    const vector shockaxis;
			    const float shockscale
			    )
{
    vector4 dPa = 0, dPr = 0;
    vector4 dF = 0;

    // If we are an infinite mass, we don't listen to constraints.
    if (mass == 0)
	return vector(0);

    float kpr = 1 - exp(-kr * 24 * 3 * timestep / niter);
    float kpa = 1 - exp(-ka * 24 * 3 * timestep / niter);
    float other_wa = 1;
    int   haswa = haspointattrib(geo, 'attractionweight');

    foreach(int ptj; neighbors)
    {
        vector pj = point(0, "P", ptj);
        float massj = point(0, massname, ptj);
        float restdist = pscale + point(0, "pscale", ptj);
	if (haswa)
	    other_wa = point(geo, 'attractionweight', ptj);
	float curdist = length(pj - pi);

        computePointDistanceDelta(dPr, dPa, 
				    pi, pj, mass, massj, curdist, restdist,
				    kpr, kpa,
				    1, other_wa,
				    shocktype, shockaxis, shockscale);
	if (frscale > 0 && curdist < restdist)
	{
	    vector prevj = point(0, "pprevious", ptj);

	    dF += computeFrictionDelta(pi, pj, pi - pprevious, pj - prevj,
					mass, massj,
					curdist, restdist,
					mu_s, mu_k,
					shocktype, shockaxis, shockscale);
	}
    }

    vector nDP = 0;
#ifndef DISABLE_CONSTRAINT_AVERAGING
    nDP += wa * vector(dPa) / dPa.w;
    nDP += sor * wr * vector(dPr) / dPr.w;
    nDP += frscale * vector(dF) / dF.w;
#else
    nDP += wa * vector(dPa);
    nDP += sor * wr * vector(dPr);
    nDP += frscale * vector(dF);
#endif

    return nDP;
}

vector pointDistanceUpdateUniformPscale(const int geo; const vector pi; 
			    const int neighbors[];
                            const float pscale, mass;
			    const string massname;
			    const float sor;
			    const float kr, wr, ka, wa;
			    const float timestep; 
			    const int niter;
			    const vector pprevious;
			    const float frscale, mu_s, mu_k;
			    const int shocktype;
			    const vector shockaxis;
			    const float shockscale
			    )
{
    vector4 dPa = 0, dPr = 0;
    vector4 dF = 0;

    // If we are an infinite mass, we don't listen to constraints.
    if (mass == 0)
	return vector(0);

    float kpr = 1 - exp(-kr * 24 * 3 * timestep / niter);
    float kpa = 1 - exp(-ka * 24 * 3 * timestep / niter);
    float other_wa = 1;
    int   haswa = haspointattrib(geo, 'attractionweight');

    foreach(int ptj; neighbors)
    {
        vector pj = point(geo, "P", ptj);
        float massj = point(geo, massname, ptj);
	if (haswa)
	    other_wa = point(geo, 'attractionweight', ptj);
        float restdist = 2 * pscale;
	float curdist = length(pj - pi);

        computePointDistanceDelta(dPr, dPa, 
				  pi, pj, mass, massj, curdist, restdist,
				  kpr, kpa,
				  1, other_wa,
				  shocktype, shockaxis, shockscale);
	if (frscale > 0 && curdist < restdist)
	{
	    vector prevj = point(geo, "pprevious", ptj);

	    dF += computeFrictionDelta(pi, pj, pi - pprevious, pj - prevj,
					mass, massj,
					curdist, restdist,
					mu_s, mu_k,
					shocktype, shockaxis, shockscale);
	}
    }

    vector nDP = 0;
#ifndef DISABLE_CONSTRAINT_AVERAGING
    nDP += wa * vector(dPa) / dPa.w;
    nDP += sor * wr * vector(dPr) / dPr.w;
    nDP += frscale * vector(dF) / dF.w;
#else
    nDP += wa * vector(dPa);
    nDP += sor * wr * vector(dPr);
    nDP += frscale * vector(dF);
#endif

    return nDP;
}

vector pointDistanceUpdateNoMass(const int geo; const vector pi; 
			    const int neighbors[];
                            const float pscale;
			    const float sor;
			    const float kr, wr, ka, wa;
			    const float timestep; 
			    const int niter;
			    const vector pprevious;
			    const float frscale, mu_s, mu_k;
			    const int shocktype;
			    const vector shockaxis;
			    const float shockscale
			    )
{
    vector4 dPa = 0, dPr = 0;
    vector4 dF = 0;

    float kpr = 1 - exp(-kr * 24 * 3 * timestep / niter);
    float kpa = 1 - exp(-ka * 24 * 3 * timestep / niter);
    float other_wa = 1;
    int   haswa = haspointattrib(geo, 'attractionweight');

    foreach(int ptj; neighbors)
    {
        vector pj = point(geo, "P", ptj);
        float restdist = pscale + point(geo, "pscale", ptj);
	if (haswa)
	    other_wa = point(geo, 'attractionweight', ptj);
	float curdist = length(pj - pi);

        computePointDistanceDelta(dPr, dPa, 
				  pi, pj, 1.0, 1.0, curdist, restdist,
				  kpr, kpa,
				  1, other_wa,
				  shocktype, shockaxis, shockscale);
	if (frscale > 0 && curdist < restdist)
	{
	    vector prevj = point(geo, "pprevious", ptj);

	    dF += computeFrictionDelta(pi, pj, pi - pprevious, pj - prevj,
					1, 1,
					curdist, restdist,
					mu_s, mu_k,
					shocktype, shockaxis, shockscale);
	}
    }

    vector nDP = 0;
#ifndef DISABLE_CONSTRAINT_AVERAGING
    nDP += wa * vector(dPa) / dPa.w;
    nDP += sor * wr * vector(dPr) / dPr.w;
    nDP += frscale * vector(dF) / dF.w;
#else
    nDP += wa * vector(dPa);
    nDP += sor * wr * vector(dPr);
    nDP += frscale * vector(dF);
#endif

    return nDP;
}
 
vector pointDistanceUpdateNoMassUniformPscale(const int geo; const vector pi; 
			    const int neighbors[];
                            const float pscale;
			    const float sor;
			    const float kr, wr, ka, wa;
			    const float timestep; 
			    const int niter;
			    const vector pprevious;
			    const float frscale, mu_s, mu_k;
			    const int shocktype;
			    const vector shockaxis;
			    const float shockscale
			    )
{
    vector4 dPa = 0, dPr = 0;
    vector4 dF = 0;

    float kpr = 1 - exp(-kr * 24 * 3 * timestep / niter);
    float kpa = 1 - exp(-ka * 24 * 3 * timestep / niter);
    float other_wa = 1;
    int   haswa = haspointattrib(geo, 'attractionweight');

    foreach(int ptj; neighbors)
    {
        vector pj = point(geo, "P", ptj);
	if (haswa)
	{
	    other_wa = point(geo, 'attractionweight', ptj);
	}
        float restdist = 2 * pscale;
	float curdist = length(pj - pi);

        computePointDistanceDelta(dPr, dPa,
				  pi, pj, 1.0, 1.0, curdist, restdist,
				  kpr, kpa,
				  1, other_wa,
				  shocktype, shockaxis, shockscale);
	if (frscale > 0 && curdist < restdist)
	{
	    vector prevj = point(geo, "pprevious", ptj);

	    dF += computeFrictionDelta(pi, pj, pi - pprevious, pj - prevj,
					1, 1,
					curdist, restdist,
					mu_s, mu_k,
					shocktype, shockaxis, shockscale);
	}
    }

    vector nDP = 0;
#ifndef DISABLE_CONSTRAINT_AVERAGING
    nDP += wa * vector(dPa) / dPa.w;
    nDP += sor * wr * vector(dPr) / dPr.w;
    nDP += frscale * vector(dF) / dF.w;
#else
    nDP += wa * vector(dPa);
    nDP += sor * wr * vector(dPr);
    nDP += frscale * vector(dF);
#endif

    return nDP;
}
 
// Do we want to take mass into consideration?
vector4 explicitConstraintUpdate(const int geo; const vector pi; 
			const int neighbors[]; const float distances[]; 
			float strains[];
			const float mass, ke, we, 
			timestep; const int niter, inststrain)
{
    // 24 is for FPS.  3 is to give 95% contraction for a stiffness
    // of 1 in a standard FPS timestep, thereby allowing for an easy
    // rule of thumb to compute the maximum effective stiffnesss
    // normalize a stiffness of 1.
    float kpe = 1 - exp(-ke * 24 * 3 * timestep / niter);
    vector4 dP = 0;
    int has_mass = haspointattrib(geo, "mass");

    // Zero mass particles do not react to constraints;
    if (mass == 0)
	return dP;

    foreach (int i; int ptj; neighbors)
    {
        if (ptj < 0)
            continue;
        vector pj = point(geo, "P", ptj);
	float massj = 1;
	if (has_mass)
	    massj = point(geo, "mass", ptj);

        float edist = distances[i]; // the correct distance constraint
        
        vector r = pi - pj;
        float curdist = length(r);
        vector grad = r / curdist;

        float C = curdist - edist;

        float weight = massj / (mass + massj);
	weight = select(massj == 0, 1, weight);
        
        vector4 dpj = -kpe * C * grad * weight;
        dP += dpj;

	// Normalized strain
	float strain = max(curdist/distances[i] - 1, 0);

	// Instantaneous strain does a direct update.
        if(inststrain)
        {
            strains[i] = strain;
        }
        else
        {
            strains[i] += strain * kpe * timestep;
        }
    }
    return dP;
}

vector jacobiFriction(const int geo, neighbors[]; 
		      const vector pi, pprevi; const float rad,
		      mass, mus, muk)
{
    float nnbr = 0;
    vector dP = 0;
    vector dpi  = pi - pprevi;

    int		hasmass = haspointattrib(0, 'mass');
    if (mass == 0)
	return dP;

    // We intentionally conflate dP with vel here to make
    // it easier to read in terms of traditional physics.
    foreach (int ptj; neighbors)
    {
        vector pj = point(0, "P", ptj);
        float massj;
        float radj = point(0, "pscale", ptj);

        vector r = pj - pi;
        float rlen = length(r);
        float mindist = (rad+radj);

        if (rlen > mindist) continue;

	if (hasmass)
	    massj = point(0, "mass", ptj);
	else
	    massj = 1;
        vector pprevj = point(0, "pprevious", ptj);

        vector dpj = pj - pprevj;
        vector relvel = dpi - dpj;

	// Find the direction between the particles.
	// This is the normal of the impact.
        r /= rlen;

	// The idea is to approximate our normal force by
	// how much we interpenetrated.
        float nmlforce = mindist - rlen;

	// Find tangent component of our relative velocity
        vector dPtan = relvel - dot(relvel, r) * r;
        float ldPtan = length(dPtan);

	// Compute kinetic friction scalar
	// this appears to violate coulomb's law.
        float fkin = select(ldPtan < 0.001, 0, muk * nmlforce / ldPtan);
        fkin = min(fkin, 1);

	// Both particles receive friction, so we bias
	// the application according to mass.
	float weight = massj / (mass + massj);
	weight = select(massj == 0, 1, weight);

	// Friction coefficient.
	// If we are in the static cone, apply the full friction,
	// ie, attempt to bring the object to a stop.
	// Otherwise use the kinetic friction.
	float fcoeff = select(ldPtan < mus * nmlforce, 1, fkin);

        dP -= weight * dPtan * fcoeff;
        nnbr++;
    }

    return dP / nnbr;
}

vector pbdSurfaceCollision(const int geo; const vector pi, pprev, hitdp, hitnml; const float mus, muk)
{
    vector dpi = pi - pprev;
    vector dp = dpi - hitdp;
    vector dPnml = dot(dp, hitnml) * hitnml;
    vector dPtan = dp - dPnml;
    float ldPtan = length(dPtan);
    float ldPnml = length(dPnml);

    // ldPnml is our approximate normal force, so we scale our tangent
    // by the ratio, clamping at 1.
    float fkin = muk * ldPnml;
    fkin = select(fkin > ldPtan, 1, fkin / ldPtan);

    // Check if we lie in the in the static cone, if so bring to
    // a full stop.
    float fcoeff = select(ldPtan < mus * ldPnml, 1, fkin);
    vector dPout = -dPtan * fcoeff;
    return dPout;
}

// Version of pointDistance that takes an array of restlengths.
vector4 pointDistanceUpdate(const int geo; const vector pi; const int neighbors[];
                           const float restlengths[], kr, wr, ka, wa; const int niter)
{
    float kpr = 1 - pow(1 - kr, 1.0 / niter);
    float kpa = 1 - pow(1 - ka, 1.0 / niter);
    vector4 dP = 0;
    foreach(int idx; int ptj; neighbors)
    {
        vector r = pi - point(0, "P", ptj);
        float curdist = length(r);
        // Constraint and gradient.
        float restdist = restlengths[idx];
        float C = curdist - restdist;
        vector gradC = r / curdist;
        // Use weighted distance constraint if within attract distance.
        if (wa > 0)
        {
            vector4 dpj = -kpa * 0.5 * C * gradC;
            dpj.w = 1;
            dP += wa * dpj;
        }
        // Use weighted (inequality) repel constraint if within repel distance.
        if (C < 0 && wr > 0)
        {
            vector4 dpj = -kpr * 0.5 * C * gradC;
            dpj.w = 1;
            dP += wr * dpj;
        }
    }
    return dP;
}

// Gauss-Seidel version
void pointDistanceUpdate(const int geo; const int pti, neighbors[];
                         const float pscale, kpr, kpa; vector curP[])
{
    float restdist = pscale * 2;
    vector4 dP = 0;
    vector pi = curP[pti];
    foreach(int ptj; neighbors)
    {
        if (ptj > pti)
            continue;
        vector r = pi - curP[ptj];
        float curdist = length(r);
        // Constraint and gradient.
        float C = curdist - restdist;
        vector gradC = r / curdist;
        // Use repel or attract stiffness.
        float kp = (curdist <= restdist ? kpr : kpa);
        vector dpj = -kp * 0.5 * C * gradC;
        pi += dpj;
        curP[ptj] -= dpj;
    }
    curP[pti] = pi;
}

// Spiky SPH kernel.
float Wspiky(const vector r; const float h)
{
    float rlen = length(r);
    if (rlen >= h || rlen < 1e-5)
        return 0;
    float h2 = h * h;
    float s = 15.0 / (M_PI * h2 * h2 * h2);
    float hr = h - rlen;
    return s * hr * hr * hr;
}

// Gradient of Spiky SPH kernel.
vector gradWspiky(const vector r; const float h)
{
    float rlen = length(r);
    if (rlen >= h || rlen < 1e-5)
        return 0;
    float h2 = h * h;
    float s = -45.0 / (M_PI * h2 * h2 * h2);
    float hr = h - rlen;
    return s * hr * hr * r / rlen;
}

// Poly6 SPH kernel.
float Wpoly6(const vector r; const float h)
{
    float r2 = length2(r);
    float h2 = h * h;
    if (r2 >= h2)
        return 0;
    float h3 = h2 * h;
    float s = 315.0 / (64.0 * M_PI * h3 * h3 * h3);
    float hr2 = h2 - r2;
    return s * hr2 * hr2 * hr2;
}

// Density kernel. Note the PBF paper uses different kernels
// for density constraint and the gradient for that constraint,
// which is pretty weird, but it looks better than using Wspiky
// for both.
float Wdensity(const vector r; const float h)
{
#if 1
    return Wpoly6(r, h);
#else
    return Wspiky(r, h);
#endif
}


// Calc SPH density for particles, assuming uniform mass.
float calcSPHDensity(const int geo; const vector pi; const int neigbors[];
                     const float h, mass)
{
    // Start with this point's density.
    float density = Wdensity(0, h);
    foreach(int ptj; neigbors)
        density += Wdensity(pi - point(geo, "P", ptj), h);
    return mass * density;
 }

// Calc scaling factor for each SPH density constraint.
float calcSPHLambda(const int geo; const vector pi; const int neighbors[];
                    const float h, mass, restdensity, eps)
{
    // Start with this point's density.
    float density = Wdensity(0, h);
    vector gradi = 0;
    float denom = 0;
    float h2 = h * h;
    foreach(int ptj; neighbors)
    {
        vector r = pi - point(geo, "P", ptj);
        // Neighbors array might be out of date.
        if (length2(r) > h2)
            continue;
        density += Wdensity(r, h);
        // Mass and h can be very small, in which case grad is very large.
        // For the moment, group these operations to attempt some numerical
        // sanity.
        vector gradj = (mass * gradWspiky(r, h)) / restdensity;
        denom += length2(gradj);
        gradi += gradj;
    }
    denom += length2(gradi);
    // denom *= (mass * mass) / (restdensity * restdensity)
    denom += eps;
    density *= mass;
    float C = density / restdensity - 1;
    float lambda = -C / denom;
    // printf("C = %g, lambda = %g, denom = %g, density = %g, restdensity = %g, gradi = %g\n", C, lambda, denom, density, restdensity, gradi);
    return lambda;
}

// Update vector from SPH density constraint.
vector calcSPHUpdate(const int geo; const vector pi; const int neighbors[];
                     const float si, h, mass, restdensity, k)
{
    vector dP = 0;
    // Aritificial pressure terms.
    // Vector with length 0.2 * h
#if 0
    vector dq = h * 0.2  / 1.7320508075688772;
#else
    vector dq = 0;  //Purely repulsive.
#endif
    int n = 4;
    float h2 = h * h;
    float Wq = Wdensity(dq, h);
    float scorr = 0;
    foreach(int ptj; neighbors)
    {
        vector r = pi - point(geo, "P", ptj);
        // Neighbors array might be out of date.
        if (length2(r) > h2)
            continue;
        float sj = point(geo, "scale", ptj);
        if (k > 0)
        {
            //scorr = -k * pow(Wdensity(r, h) / Wq, n);
            // Assume n = 4 always and avoid the pow call.
            scorr = Wdensity(r, h) / Wq;
            scorr *= scorr;
            scorr *= scorr;
            scorr *= -k;
        }
        vector gradj = (mass * gradWspiky(r, h)) / restdensity;
        dP += (si + sj + scorr) * gradj;
    }
    return dP;
}

vector calcXSPHViscosity(const int geo; const vector pi, vi; const int neighbors[];
                         const float h, mass, viscosity, confinement)
{
    vector dv = 0;
    foreach(int ptj; neighbors)
    {
        vector pj = point(geo, "P", ptj);
        vector vj = point(geo, "v", ptj);
        vector dj = point(geo, "density", ptj);
        vector vij = vi - vj;
        if (viscosity > 0)
        {
            // XSPH viscosity
            dv -= viscosity * (mass * vij) * (Wdensity(pi - pj, h) / dj);
        }
        if (confinement > 0)
        {
            //vector theta = cross(vij, gradWspiky(pi - pj, h));
            //dV += 
        }
    }
    return dv;
}

// Stays within 0.1% of raised inverted cosine sigmoid function
float blinnWyvillSigmoid(const float x)
{
    float x2 = x*x;
    float x4 = x2*x2;
    float x6 = x4*x2;

    float fa = ( 4.0/9.0);
    float fb = (17.0/9.0);
    float fc = (22.0/9.0);

    float y = fa*x6 - fb*x4 + fc*x2;
    return y;
}

// Stays within 2.8% of raised inverted cosine sigmoid function
float symmetricDoubleQuadraticSigmoid(const float x)
{
    float _2x = 2 * x;
    float y1 = _2x * x;
    float y2 = -y1 + _2x - 1;
    float y = select(x < 0.5, y1, y2);
    return y;
}

// use this in cases when a sigmoid might have inputs not in [0,1]
float clampSigmoid(const float x, y)
{    
    float out = select(x < 0, 0, y);
    out = select(x > 1, 1, out);
    return out;
}

// Simple stiff, one-way distance inequality constraint that can be used
// to limit the amount a particle can travel in one timestep, useful
// for speed limiting among other things.
vector limitDistance(const vector pos, prev; const float maxdist)
{
    vector r = pos - prev;
    float C = length(r) - maxdist;
    if (C > 0)
        return -C * normalize(r);
    return 0;
}

vector collisionUpdate(const vector pos, hitpos, hitnml)
{
    // Assumes hitnml is normalized.
    float C = dot(pos - hitpos, hitnml);
    if (C < 0)
        return -C * hitnml;
    return 0;
}

// Calc velocity based on current and previous
vector frictionUpdate(const vector v, hitnml; float friction)
{
    // Assumes hitnml is normalized.
    vector tanv = v - dot(v, hitnml) * hitnml;
    return -friction * tanv;
}

// Collision response. Handles friction but zero restitution for now.
vector collisionResponse(const vector v, hitv, hitnml; float friction)
{
    // Assumes hitnml is normalized.
    vector relv = v - hitv;
    vector normv = dot(relv, hitnml) * hitnml;
    vector tanv = relv - normv;
    vector collnormv = dot(hitv, hitnml) * hitnml;
    // We only include friction in this update, the clamping of velocity
    // is implicitly handled by updating velocity to match the actual
    // change in position, and we have already clamped position.
    // Trying to clamp velocity here may actually speed things up sticking
    // particles to trailing edges.
    return -friction * tanv; // + collnormv - normv;
}

float airFriction(const int geo, neighbors[]; const vector pi, vi; const float rad, mass, airresist)
{
    // calculate the degree of isolation of the sand particle
    int nnbr = 0;
    float subtend_area = 0; // not really, but kinda
    foreach(int ptj; neighbors)
    {
        vector r = pi - point(geo, "P", ptj);
        float rlen = length(r);
        subtend_area += cos(asin(rad / rlen));
        ++nnbr; // assume each neighbor provides equal cover
    }
    // float exp_spd = 1 + spd;// + spd * spd * (1 + spd / 3) / 2;
    float isolation = 1 / (1 + subtend_area);
    return airresist * isolation;// * exp_spd;
}

// Simple inter-particle drag friction velocity adjustment.
vector particleDragFriction(const int geo; const vector pi, vi; const int neighbors[];
                            const float rad, friction)//, staticthreshhold)
{
    vector dV = 0;
    int nnbr = 0;
    foreach(int ptj; neighbors)
    {
        vector r = pi - point(geo, "P", ptj);
        float rlen = length(r);
        if (rlen > rad)
            continue;
        r /= rlen;
        vector vj = point(geo, "v", ptj);
        vector vij = vi - vj;
        vector tanv = vij - dot(vij, r) * r;
        //float fric = select(length(tanv) > staticthreshhold, friction, 1);
        //dV -= fric * tanv;
        dV -= friction * tanv;
        ++nnbr;
    }
    return dV / nnbr;
}

// computes the center of mass of the passed in points
vector computeCenterOfMass(int geo; int points[])
{
    vector com = 0;
    float total_mass = 0;
    int has_mass = haspointattrib(geo, "mass");

    foreach(int idx; int pti; points)
    {
        vector P = point(geo, "P", pti);
	float mass = 1;
	if (has_mass)
	    mass = point(geo, "mass", pti);
        
        com += mass * P;
        total_mass += mass;
    }

    com /= total_mass;

    return com;
}

// for each passed in point, returns the offset from
// the center of mass of the collection of points
vector[] computeCenterOfMassOffset(int geo; int points[])
{
    vector out[];
    vector com = computeCenterOfMass(geo, points);

    foreach(int idx; int pti; points)
    {
        vector P = point(geo, "P", pti);
        
        out[idx] = P - com;
    }

    return out;
}

vector[] matchShape(vector pos[], restpos[]; float mass[])
{
    vector dP[];
    matrix3 A = {0,0,0,0,0,0,0,0,0};
    
    vector c = 0, cr = 0;
    int npt = len(pos);
    float mass_total = 0.0; 

    for(int i=0; i<npt; i++)
    {
        vector pi = pos[i];
        vector ri = restpos[i];
        mass_total += mass[i];
        
        c += mass[i] * pi;
        cr += mass[i] * ri;
    }
    
    c /= mass_total;
    cr /= mass_total;
    
    for(int i=0; i<npt; i++)
    {
        vector pi = pos[i] - c;
        vector ri = restpos[i] - cr;
        
        A += outerproduct(pi, ri);
    }
    
    matrix3 Q = polardecomp(A);
    
    for(int i=0; i<npt; i++)
    {
        vector pi = pos[i] - c;
        vector ri = restpos[i] - cr;
        
        dP[i] = (Q * ri + c) - pos[i];
    }

    return dP;
}

// return an array of dPs such that pos[i]+dP[i] will 
// put the points[] back in their rest shape,
// while mainting overall translation and rotation
vector[] matchShape(int geo; int points[])
{
    float mass[];
    vector pos[];
    vector restpos[];

    int has_mass = haspointattrib(0, "mass");
    
    foreach(int idx; int pti; points)
    {
        mass[idx] = has_mass ? point(0, "mass", pti) : 1.0;
        pos[idx] = point(0, "P", pti);
        restpos[idx] = point(0, "rest", pti);
    }
    
    vector dP[] = matchShape(pos, restpos, mass);
    
    return dP;
}

#endif
