#ifndef __pbd_constraints_h
#define __pbd_constraints_h


void
createDistanceConstraint(const int geo; const int ptnum; const string srcgrp;
                         const int outgeo; const string outgrp)
{
    int nbrs[] = neighbours(geo, ptnum);
    vector ptpos = point(geo, "P", ptnum);
    foreach(int n; nbrs)
    {
        if (n <= ptnum || !inpointgroup(geo, srcgrp, n))
            continue;
        int prim = addprim(outgeo, "polyline", array(ptnum, n));
        setprimgroup(outgeo, outgrp, prim, 1);
        setprimattrib(outgeo, "restlength", prim, distance(ptpos, point(geo, "P", n)));
    }
}

int
oppositepoint(const int geo; const int hedge)
{
    return hedge_dstpoint(geo, hedge_next(geo, hedge));
}

void
createBendDistanceConstraint(const int geo; const int primnum)
{
    int starthedge = primhedge(geo, primnum);

    if (starthedge < 0)
        return;

    int h = starthedge;
    int p = primnum;
    int type = 11;
    while (1)
    {
        // Giving half edge h, add a bend polygon.
        int oh = hedge_nextequiv(geo, h);
        if (h != oh && oh >= 0)
        {
            int op = hedge_prim(geo, oh);
            // Always build in ascending direction.
            if (op >= 0 && p < op)
            {
                int pt0 = oppositepoint(geo, h);
                int pt1 = oppositepoint(geo, oh);
                int prim = addprim(geoself(), 'polyline', array(pt0, pt1));
                setprimattrib(geoself(), "type", prim, type);
            }
        }

        int nh = hedge_next(geo, h);
        
        if (nh == starthedge)
            break;
        h = nh;
    }
}

void
addTriangleFEMConstraint(const int geo; const int primnum)
{
	int type = 6;
	int pts[] = primpoints(geo, primnum);
	vector p0 = point(geo, "P", pts[0]);
	vector p1 = point(geo, "P", pts[1]);
	vector p2 = point(geo, "P", pts[2]);

	vector normal0 = cross(p1 - p0, p2 - p0);
	float area = length(normal0) * 0.5;

	vector axis0_1 = normalize(p1 - p0);
	vector axis0_2 = normalize(cross(normal0, axis0_1));

	vector2 p[];
	resize(p, 3);
	p[0] = set(dot(p0, axis0_2), dot(p0, axis0_1));
	p[1] = set(dot(p1, axis0_2), dot(p1, axis0_1));
	p[2] = set(dot(p2, axis0_2), dot(p2, axis0_1));

	matrix2 P;
	setcomp(P, p[0][0] - p[2][0], 0, 0);
	setcomp(P, p[0][1] - p[2][1], 1, 0);
	setcomp(P, p[1][0] - p[2][0], 0, 1);
	setcomp(P, p[1][1] - p[2][1], 1, 1);

	matrix2 invRestMat = invert(P);
    setprimattrib(geoself(), "restlength", primnum, area);            
	setprimattrib(geoself(), "restmatrix", primnum, invRestMat);
    setprimattrib(geoself(), "type", primnum, type);

}

void
createDihedralConstraint(const int geo; const int primnum;
                         const int outgeo; const string outgrp)
{
    int starthedge = primhedge(geo, primnum);

    if (starthedge < 0)
        return;

    // Ignore open curves as the hedge function won't give us
    // the setup we expect with them.
    if (primintrinsic(geo, 'closed', primnum) == 0)
        return;

    int h = starthedge;
    int p = primnum;
    int type = 1;
    int hasgrp = strlen(outgrp) > 0;
    while (1)
    {
        // Giving half edge h, add a bend polygon.
        int oh = hedge_nextequiv(geo, h);
	// Skip boundary edges, invalid hedges, and non-manifold
	// edges
        if (h != oh && oh >= 0 && h == hedge_nextequiv(geo, oh))
        {
            int op = hedge_prim(geo, oh);
            // Always build in ascending direction.
            if (op >= 0 && p < op)
            {
                int pt0 = oppositepoint(geo, h);
                int pt1 = oppositepoint(geo, oh);
                int pt2 = hedge_srcpoint(geo, h);
                int pt3 = hedge_dstpoint(geo, h);

                vector p0 = point(geo, "P", pt0);
                vector p1 = point(geo, "P", pt1);
                vector p2 = point(geo, "P", pt2);
                vector p3 = point(geo, "P", pt3);
                
                vector e = p3 - p2;
                float elen = length(e);
                if (elen < 1e-6)
                    continue;
                float invElen = 1 / elen;
                
                // Find initial rest angle.
                vector n1 = cross(p3 - p0, p2 - p0);
                vector n2 = cross(p2 - p1, p3 - p1);
                float d = dot(normalize(n1), normalize(n2));
                d = clamp(d, -1, 1);
                float phi = acos(d);
                // We want to xpress phi a -PI..PI
                if (dot(cross(n1, n2), e) < 0)
                	phi = -phi;
                int prim = addprim(outgeo, 'polyline', array(pt0, pt1, pt2, pt3));
                setprimattrib(outgeo, "restlength", prim, phi);
                setprimattrib(outgeo, "type", prim, type);
                if (hasgrp)
                    setprimgroup(outgeo, outgrp, prim, 1);
            }
        }

        int nh = hedge_next(geo, h);
        
	// Stop the loop when we complete the polygon
	// Closed polygons won't complete.
        if (nh == starthedge || nh < 0)
            break;
        h = nh;
    }
}

void
createDihedralConstraint(const int geo; const int primnum)
{
    createDihedralConstraint(geo, primnum, geoself(), "");
}

void
computeOrientInertia(const int geo; const int primnum; const float defmass; const string srcgrp, pingrp;
                          const int outgeo)
{
    // Ignore anything but open polylines.
    if (primintrinsic(geo, "typename", primnum) != "Poly" || 
        primintrinsic(geo, "closed", primnum) == 1)
        return;
    // Only test group if it's not all points.
    int hasgrp = npointsgroup(geo, srcgrp) < npoints(geo);
    vector from = {0, 0, 1};
    int pts[] = primpoints(geo, primnum);
    int npts = len(pts);
    vector4 orients[];
    float rodlens[];
    resize(orients, npts - 1);
    resize(rodlens, npts - 1);
    for(int i=0; i < npts - 1; i++)
    {
        vector d = point(0, "P", pts[i + 1]) - point(0, "P", pts[i]);
        vector to = normalize(d);
        vector4 dq = dihedral(from, to);
        if (i == 0)
            orients[i] = dq;
        else
            orients[i] = qmultiply(dq, orients[i-1]);
        rodlens[i] = length(d);
        from = to;
    }

    for(int i=0; i < npts - 1; i++)
    {
        if (hasgrp && !inpointgroup(geo, srcgrp, pts[i]))
            continue;
        setpointattrib(geoself(), "orient", pts[i], orients[i]);
        float mass = defmass;
        if (inpointgroup(geo, pingrp, pts[i]))
            mass = 0;
        setpointattrib(geoself(), "inertia", pts[i], mass * 2 * rodlens[i] / 5);
    }
    // Set inertia for last point as well.
    int lastpt = pts[npts-1];
    if (npts > 1 && inpointgroup(geo, srcgrp, lastpt))
    {
        float mass = defmass;
        if (inpointgroup(geo, pingrp, lastpt))
            mass = 0;
        setpointattrib(geoself(), "inertia", lastpt, mass * 2 * rodlens[npts - 2] / 5);    
    }
}

void
createBendTwistConstraint(const int geo; const int ptnum; const string srcgrp;
                          const int outgeo; const string outgrp)
{
    int prims[] = pointprims(geo, ptnum);
    // Ignore anything but open polylines.
    if (len(prims) > 1)
        return;
    int primnum = prims[0];        
    if (primintrinsic(geo, "typename", primnum) != "Poly" || 
        primintrinsic(geo, "closed", primnum) == 1)
        return;
    int nbrs[] = neighbours(geo, ptnum);
    int n = max(nbrs);
    // Stop at end of line.
    if (n < ptnum || !inpointgroup(geo, srcgrp, n))
        return;
    vector4 q0 = point(geo, "orient", ptnum);
    vector4 q1 = point(geo, "orient", n);    
    vector4 restDarbeaux = qmultiply(set(-q0.x, -q0.y, -q0.z, q0.w), q1);
    vector4 omegaplus = restDarbeaux + {0, 0, 0, 1};
    vector4 omegaminus = restDarbeaux - {0, 0, 0, 1};
    if (dot(omegaminus, omegaminus) > dot(omegaplus, omegaplus))
        restDarbeaux *= -1;

    int prim = addprim(outgeo, "polyline", array(ptnum, n));
    setprimattrib(outgeo, "restvector", prim, restDarbeaux);
    setprimgroup(outgeo, outgrp, prim, 1);
}

float cotTheta(const vector v, w)
{
    float cosTheta = dot(v, w);
    float sinTheta = length(cross(v, w));
    return (cosTheta / sinTheta);
}

void
createIsometricConstraint(const int geo; const int primnum)
{
    int starthedge = primhedge(geo, primnum);

    if (starthedge < 0)
        return;

    int h = starthedge;
    int p = primnum;
    int type = 2;
    while (1)
    {
        // Giving half edge h, add a bend polygon.
        int oh = hedge_nextequiv(geo, h);
        if (h != oh && oh >= 0)
        {
            int op = hedge_prim(geo, oh);
            // Always build in ascending direction.
            if (op >= 0 && p < op)
            {
                int pt0 = oppositepoint(geo, h);
                int pt1 = oppositepoint(geo, oh);
                int pt2 = hedge_srcpoint(geo, h);
                int pt3 = hedge_dstpoint(geo, h);

                vector p0 = point(geo, "P", pt0);
                vector p1 = point(geo, "P", pt1);
                vector p2 = point(geo, "P", pt2);
                vector p3 = point(geo, "P", pt3);
                
                vector x[] = array(p2, p3, p0, p1);

                vector e0 = x[1] - x[0];
                vector e1 = x[2] - x[0];
                vector e2 = x[3] - x[0];
                vector e3 = x[2] - x[1];
                vector e4 = x[3] - x[1];

                float c01 = cotTheta(e0, e1);
                float c02 = cotTheta(e0, e2);
                float c03 = cotTheta(-e0, e3);
                float c04 = cotTheta(-e0, e4);
                
                float A0 = 0.5 * length(cross(e0, e1));
                float A1 = 0.5 * length(cross(e0, e2));

                float coef = -3.f / (2.f*(A0 + A1));
#if 0
                float K[] = array( c03 + c04, c01 + c02, -c01 - c03, -c02 - c04 );
                float K2[] = array(  coef*K[0], coef*K[1], coef*K[2], coef*K[3] );

                matrix Q = 0;
                for (int j = 0; j < 4; j++)
                {
                    for (int k = 0; k < j; k++)
                    {
                        float KK2 = K[j] * K2[k];
                        setcomp(Q, KK2, k, j);
                        setcomp(Q, KK2, j, k);
                    }
                    setcomp(Q, K[j] * K2[j], j, j);
                }
#else
                vector4 K = set(c03 + c04, c01 + c02, -c01 - c03, -c02 - c04 );
				matrix Q = coef * outerproduct(K, K);
#endif
                int prim = addprim(geoself(), 'polyline');
                addvertex(geoself(), prim, pt0);
                addvertex(geoself(), prim, pt1);
                addvertex(geoself(), prim, pt2);
                addvertex(geoself(), prim, pt3);
                setprimattrib(geoself(), "restmatrix", prim, Q);            
                setprimattrib(geoself(), "type", prim, type);
            }
        }

        int nh = hedge_next(geo, h);
        
        if (nh == starthedge)
            break;
        h = nh;
    }
}


int[]
findAttachments(const int geo; const string pingrp; const string classattr;
				const int ptnum)
{
	int attachments[];
	int pins[] = expandpointgroup(geo, pingrp);
	int type = 7;
	if (find(pins, ptnum) >= 0)
		return attachments;
	foreach(int pin; pins)
	{
		if (pin == ptnum)
			continue;
		// Make sure connectivity class is equal.
		int pinclass = point(geo, classattr, pin);
		int pclass = point(geo, classattr, ptnum);
		if (pinclass != pclass)
			continue;
		append(attachments, pin);
	}
	return attachments;
}

void
createAttachmentConstraint(const int geo; const string pingrp; const string classattr;
						   const int ptnum; float kstiff)
{
	int attachments[] = findAttachments(geo, pingrp, classattr, ptnum);
	int type = 7;
	foreach(int pin; attachments)
	{
	    int prim = addprim(geoself(), 'polyline', array(pin, ptnum));
	    setprimattrib(geoself(), "stiffness", prim, kstiff);
	    setprimattrib(geoself(), "type", prim, type);
	    vector pinp = point(geo, "P", pin);
	    vector p = point(geo, "P", ptnum);
	    float restlen = distance(pinp, p);
	    setprimattrib(geoself(), "restlength", prim, restlen);
	}
}

int []findConnectedComponents(const int geo; const string group)
{
	int pins[] = expandpointgroup(geo, group);
	int comp[];
	resize(comp, npoints(geo));
	int npins = len(pins);
	int curcomp = 1;
	for(int i=0; i < npins; i++)
	{
	    int pin = pins[i];
	    // Already in connected component.
	    if (comp[pin])
	        continue;
	    int searches[];
	    push(searches, pin);

	    while(len(searches))
	    {
	        int cur = pop(searches);
	        comp[cur] = curcomp;
	        int ns[] = neighbours(geo, cur);
	        foreach(int n; ns)
	        {
	            if(comp[n])
	                continue;
	            if(!inpointgroup(geo, group, n))
	                continue;
	            push(searches, n);            
	        }
	    }
	    curcomp++;
	}
	return comp;
}

#endif
