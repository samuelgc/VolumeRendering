#ifndef __groom_h__
#define __groom_h__

/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *    Side Effects Software Inc
 *    477 Richmond Street West
 *    Toronto, Ontario
 *    Canada   M5V 3E7
 *    416-504-9876
 *
 * NAME:    groom.h ( VOP Library Functions, VEX )
 *
 * COMMENTS:    Support functions for fur grooming tools.
 */
#include <math.h>
#include <file.h>

void adjustPrimLength(const int geo, prim; const float currentlength, targetlength)
{
    float diff = targetlength - currentlength;

    int points[] = primpoints(geo, prim);

    if(diff > 0)
    {
        vector posA = point(geo, "P", points[-2]);
        vector posB = point(geo, "P", points[-1]);

        vector posC = diff * normalize(posB-posA) + posB;

        int ptnum = addpoint(geo, posC);

        addvertex(geo, prim, ptnum);
    }
    else if(diff < 0)
    {
        vector lastpos = 0;
        float length = 0;
        int cut = 0;

        foreach(int idx; int pt; points)
        {
            vector pos = point(geo, "P", pt);
            if(idx > 0 && !cut)
            {
                float seglength = distance(pos, lastpos);

                if(length+seglength > targetlength)
                {
                    float frac = (targetlength-length)/seglength;

                    vector newpos = lerp(lastpos, pos, frac);

                    int ptnum = addpoint(geo, newpos);

                    for(int i=idx; i<len(points); i++)
                        removepoint(geo, points[i]);

                    addvertex(geo, prim, ptnum);

                    break;
                }

                length += seglength;
            }

            lastpos = pos;
        }
    }
}

// creates a guide at 'sourcepos' using 'surfacenormal' as the skin normal.
// the guide is created by interpolating between the guides in refinput.
// the guides and weights parameters define which guides from the refinput to use
// and how to weigh the interpolation
// the refinput is expected to have a 'skinN' primitive attribute defining
// the skin normal at the primitive's root.
void createGuideFromReference(int input; vector sourcepos, sourcenormal; int refinput; int guides[]; float weights[])
{
    int inherit_Cd = haspointattrib(refinput, "Cd");
    if(inherit_Cd)
        addpointattrib(input, "Cd", {0,0,0});

    // get the vertex count of the first guide
    // we'll assume that all guides have the same vertex count.
    int vtxcount = primvertexcount(refinput, guides[0]);

    vector primcolor = 0;
    float wsum = 0;

    int prim = addprim(input, "polyline");

    for(int i=0; i<vtxcount; i++)
    {
        wsum = 0;
        vector hairpos = {0,0,0};
        vector haircolor = {0,0,0};

        foreach(int g; int gprim; guides)
        {
            int guidept = primpoint(refinput, gprim, i);
            vector pos = point(refinput, "P", guidept);
            vector color = point(refinput, "Cd", guidept);

            int rootpt = primpoint(refinput, gprim, 0);
            vector rootpos = point(refinput, "P", rootpt);
            vector skinN = prim(refinput, "skinN", gprim);

            matrix3 rot = dihedral(skinN, sourcenormal);

            vector localpos = (pos - rootpos) * rot;

            hairpos += weights[g] * localpos;
            haircolor += weights[g] * color;
            wsum += weights[g];
        }

        float invweight = 1.0/wsum;
        int pt = addpoint(input, sourcepos + invweight * hairpos);
	if(inherit_Cd)
	    setpointattrib(input, "Cd", pt, invweight * haircolor);

        addvertex(input, prim, pt);
    }

}

// find the maxguides closest guides within maxradius
// and create a new guide using createGuideFromReference
void createGuideFromNeighbors(int input; vector sourcepos, sourcenormal; int refinput; float maxradius; int maxguides)
{
    int nearpts[] = nearpoints(refinput, "roots", sourcepos, maxradius, maxguides);
    int nearpt = nearpts[0];

    if(len(nearpts) > 0)
    {
        float weights[];
        int guides[];

        float wsum = 0;
        foreach(int idx; int pt; nearpts)
        {
            append(guides, pointprims(refinput, pt)[0]);

            vector pj = point(refinput, "P", pt);
            float dist = distance(sourcepos, pj);

            float w = 1.0/(dist*dist);
            append(weights, w);
            wsum += w;
        }

        foreach(int idx; float w; weights)
        {
            weights[idx] /= wsum;
        }

        createGuideFromReference(input, sourcepos, sourcenormal, refinput, guides, weights);
    }
}

// find the 4 closest guides
// and create a new guide using createGuideFromReference
void createGuideFromNeighbors(int input; vector sourcepos, sourcenormal; int refinput)
{
    createGuideFromNeighbors(input, sourcepos, sourcenormal, refinput, 1e12, 4);
}

vector constrain_surface(const vector pos, rest; const int volume_input; float lift; float strength; float error)
{
    float locallift = lift;
    float sdist = volumesample(volume_input, "surface", pos);
    float restlift = volumesample(volume_input, "surface", rest);

    if(locallift > restlift)
        locallift = lerp(restlift, locallift, min(strength, 1.0));

    sdist -= locallift;

    if(sdist < 0)
    {

        vector dir = normalize(volumegradient(volume_input, "surface", pos));

        error = -sdist;//abs(sdist/locallift);

        /*if(error < 1)
            error = 1.0/error-1;
        else
            error = error-1;*/
        return -sdist * dir;
    }
    else
    {
        error = 0;
        return 0;
    }

}

vector constrain_len(vector pos, lastpos, rest, lastrest; float error)
{
    float restlen = distance(rest, lastrest );

    vector dir = normalize(pos - lastpos);

    return lastpos + restlen * normalize(dir) - pos;
}


void constrain(vector positions[]; vector rest[]; vector moved[];
               const int volume_input; int do_length; float lift)
{
    float len = 0;
    vector dp;

    float error;
    foreach(int i; vector pos; positions)
    {
        if(i==0)
            continue;

        len += distance(rest[i-1], rest[i]);

        int niter = 20;
        for(int j=0; j<niter; j++)
        {
            error = 0;
            // length first, otherwise first iteration of surface constraint is
            // usually useless
            if(do_length)
                positions[i] += constrain_len(positions[i], positions[i-1], rest[i], rest[i-1], error);

            float pointlift = len * sin(0.5 * M_PI * lift);
            float strength = 2 * length(moved[i]);
            positions[i] += constrain_surface(positions[i], rest[i], volume_input, pointlift, strength, error);

            //setpointattrib(0, "Cd", points[i], (float)j/10 * {1,1,1});

            if(error < 0.01)
                break;
        }
        // constrain length at end, this enforces 0 stretch error
        if(do_length)
            positions[i] += constrain_len(positions[i], positions[i-1], rest[i], rest[i-1], error);

    }
}

float compute_blend_value(int stroke_curve_input; float in_u, strength)
{
        float u = clamp(in_u, 0, 1);
        float falloff[] = prim(stroke_curve_input, "stroke_falloff", 0);
        float blend = u >= 1.0 ? 0.0 : spline("linear", u, falloff);
        blend *= strength;
        blend = clamp(blend, 0, 1);

        return blend;
}

#define ADVECT_POINTS(inputtype) \
void advect_points(inputtype skin_file; vector positions[]; vector vel[]; const int points[]; \
        int lock_root; vector guideorigin) \
{ \
    foreach(int i; int pt; points) \
    { \
        if(i == 0 && lock_root) \
            continue; \
 \
        int visible = point(0, "visible", pt); \
 \
        if(!visible) \
            continue; \
 \
        positions[i] += vel[i]; \
 \
        if(i==0) \
        { \
            int prim; \
            vector primuv; \
            xyzdist(skin_file, positions[i], prim, primuv); \
            positions[i] = primuv(skin_file, "P", prim, primuv); \
            guideorigin = positions[i]; \
        } \
    } \
} \

ADVECT_POINTS(string)
ADVECT_POINTS(int)

void accumulate_segv(int stroke_curve_input; int stroke_geo_input; vector positions[];
        const int iter; const float strength, maxradius; vector segv[])
{
    string pattern = sprintf("\@seg==%d", iter);

    foreach(int i; vector pos; positions)
    {
        int prim;
        vector primuv;
        float dist = xyzdist(stroke_geo_input, pattern, pos, prim, primuv, maxradius);

        //grab values from start of segment
        primuv.y = 0;

        float radius = primuv(stroke_geo_input, "radius", prim, primuv);
        float blend = prim == -1 ? 0 : compute_blend_value(stroke_curve_input, dist/radius, strength);


        vector tempsegv = primuv(stroke_geo_input, "segv", prim, primuv);
        segv[i] += blend * tempsegv;

    }
}

// returns the parametric coordinate on the line segment (a,b)
// closest to point p
float get_pos_on_line(vector a, b, p)
{
  vector a_to_p = p - a;
  vector a_to_b = b - a;

  float atb2 = length2(a_to_b);

  float atp_dot_atb = dot(a_to_p, a_to_b);

  float t = atp_dot_atb / atb2;

  return clamp(t, 0, 1);
}

// On the given stroke segment, retuns the position closest to 'pos'.
// Sets p0-p3, which represent the corner points of the quad on which the
// closest position was found.
// Sets u and v to the parametric coordinates within the closest primitive.
vector get_pos_bilinear(vector positions[]; vector pos; int offset; int ncurvepts; int seg; float u; float v; int p0, p1, p2, p3)
{
    p0 = seg + offset;
    p1 = seg + ncurvepts + offset;
    p2 = p0 + 1;
    p3 = p1 + 1;

    vector pa1 = positions[p0];
    vector pa2 = positions[p1];
    vector pb1 = positions[p2];
    vector pb2 = positions[p3];

    vector pa = lerp(pa1, pa2, u);
    vector pb = lerp(pb1, pb2, u);

    v = get_pos_on_line(pa, pb, pos);

    return lerp(pa, pb, v);
}

// Returns the bilinear interpolation at (u,v) of values[] indexed by p0-p3
//
// Expected point order:
//
// p0---p1
// |     |
// p2---p3
//
// u=0, v=0 is the top left corner
vector interpolate_value(vector values[]; int p0, p1, p2, p3; float u; float v)
{
    vector va1 = values[p0];
    vector va2 = values[p1];
    vector vb1 = values[p2];
    vector vb2 = values[p3];

    vector va = lerp(va1, va2, u);
    vector vb = lerp(vb1, vb2, u);

    return lerp(va, vb, v);
}

// See above
float interpolate_value(float values[]; int p0, p1, p2, p3; float u; float v)
{
    float va1 = values[p0];
    float va2 = values[p1];
    float vb1 = values[p2];
    float vb2 = values[p3];

    float va = lerp(va1, va2, u);
    float vb = lerp(vb1, vb2, u);

    return lerp(va, vb, v);
}

// returns the normalized depth of 'pos' within the stroke geometry.
float get_stroke_u(int stroke_curve_input; int stroke_geo_input; int stroke_num; vector pos)
{
    vector stroke_projdir = prim(stroke_curve_input, "stroke_projdir", stroke_num);
    int first_pt = primpoint(stroke_curve_input, stroke_num, 0);
    vector stroke_orig = point(stroke_curve_input, "stroke_orig", first_pt);

    vector offset = pos - stroke_orig;
    float dp = dot(offset, normalize(stroke_projdir));

    float u = fit(dp, detail(stroke_geo_input, "mindist"), detail(stroke_geo_input, "maxdist"), 0, 1);

    return u;
}

// fetches positions and radii from the stroke geometry's points.
void get_stroke_data(int stroke_curve_input, stroke_geo_input; int stroke_num; vector pos; vector positions[]; float radii[]; int ncurvepts)
{
    for(int i=ncurvepts*2*stroke_num; i<ncurvepts*2*(stroke_num+1); i++)
    {
        positions[i] = point(stroke_geo_input, "P", i);
        radii[i] = point(stroke_geo_input, "radius", i);
    }
}

vector closest_stroke_pos(int stroke_curve_input, stroke_geo_input; vector pos; float in_maxradius; float extratol; int hitstroke, hitprim; float closestu, closestv)
{
    vector positions[];
    float radii[];
    float maxradius = in_maxradius + extratol;

    int ncurves = nprimitives(stroke_curve_input);
    int ncurvepts = npoints(stroke_geo_input) / ncurves / 2;
    int nsegs = ncurvepts - 1;

    vector closestpos = 0;
    int p0, p1, p2, p3;
    int cp0 = 0, cp1 = 0, cp2 = 0, cp3 = 0;
    int found = 0;
    float mindist = 1e9;

    for(int stroke_num=0; stroke_num<ncurves; stroke_num++)
    {
        int pointoffset = stroke_num * 2 * ncurvepts;

        float u = get_stroke_u(stroke_curve_input, stroke_geo_input, stroke_num, pos);
        get_stroke_data(stroke_curve_input, stroke_geo_input, stroke_num, pos, positions, radii, ncurvepts);

        for(int i=0; i<nsegs; i++)
        {
            float v;
            vector foundpos = get_pos_bilinear(positions, pos, pointoffset, ncurvepts, i, u, v, p0, p1, p2, p3);

            float dist = distance(pos, foundpos);

            int cond = dist < maxradius && dist < mindist;

            mindist = select(cond, dist, mindist);
            closestpos = select(cond, foundpos, closestpos);
            cp0 = select(cond, p0, cp0);
            cp1 = select(cond, p1, cp1);
            cp2 = select(cond, p2, cp2);
            cp3 = select(cond, p3, cp3);
            hitprim = select(cond, stroke_num * nsegs + i, hitprim);
            hitstroke = select(cond, stroke_num, hitstroke);
            closestu = select(cond, u, closestu);
            closestv = select(cond, v, closestv);
            found = select(cond, 1, found);
        }
    }

    float radius = interpolate_value(radii, cp0, cp1, cp2, cp3, closestu, closestv);

    radius += extratol;

    int final_cond = found && mindist < radius;
    hitstroke = select(final_cond, hitstroke, -1);
    hitprim = select(final_cond, hitprim, -1);
    return select(final_cond, closestpos, {0,0,0});
}

// old signature for backwards compatibility
int closest_stroke_pos(int stroke_curve_input, stroke_geo_input; vector pos; float in_maxradius; float extratol; vector closestpos; float closestu, closestv)
{
    int hitstroke = -1;
    int hitprim = -1;

    vector closest_pos = closest_stroke_pos(stroke_curve_input, stroke_geo_input, pos, in_maxradius, extratol, hitstroke, hitprim, closestu, closestv);

    return hitprim;
}

// Returns the length of the path defined by points[].
float pathlength(int input; int points[])
{
    float length = 0;
    vector lastpos = 0;
    foreach(int i; int pt; points)
    {
        vector pos = point(input, "P", pt);

        length += select(i > 0, distance(lastpos, pos), 0.0);

        lastpos = pos;
    }

    return length;
}

// Returns 1 if pos is in input's bounding box.
// Grows the bounding box by 'tol'.
int inbbox(int input; vector pos; float tol)
{
    vector bbmin;
    vector bbmax;
    getbbox(input, bbmin, bbmax);

    int inbb = pos.x > bbmin.x - tol && pos.x < bbmax.x + tol &&
       pos.y > bbmin.y - tol && pos.y < bbmax.y + tol &&
       pos.z > bbmin.z - tol && pos.z < bbmax.z + tol;

    return inbb;
}


vector blendsegmentdir(vector segvec; float tangentblend, normalblend; vector targetdir, planenormal)
{
    matrix3 rot = ident();

    // Compute direction and direction projected to skin plane
    vector segdir = normalize(segvec);

    // Use targetdir as upvector so it is the 0 direction
    // unless it's too close to the skin normal. in this case direction
    // blending will not do anything, because the 0 direction lines up
    // with the target direction.
    int usesegdir = dot(targetdir, planenormal) > 0.99999;
    matrix3 skinspace = lookat({0,0,0}, -planenormal, select(usesegdir, segdir, targetdir));
    matrix3 skinspace_inv = invert(skinspace);
    vector segdir_local = segdir * skinspace_inv;

    float currdirangle = atan2(-segdir_local.x, segdir_local.y);

    // Get shortest rotation
    float diff = currdirangle;

    if (segdir_local.z > 0.99999)
        diff = 0;
    else if (currdirangle > M_PI)
        diff = currdirangle - M_TWO_PI;
    else if (currdirangle < -M_PI)
        diff = currdirangle + M_TWO_PI;
    
    // rotate by direction factor
    rotate(rot, tangentblend * diff, {0,0,-1});

    vector targetdir_local = targetdir * skinspace_inv;

    float currliftangle = acos(segdir_local.z);
    float targetliftangle = acos(targetdir_local.z);

    // clamp lift angle to outer surface
    targetliftangle = min(targetliftangle, 0.49 * M_PI);

    // lift axis perpendicular to rotation
    vector liftaxis = set(cos(diff), sin(diff), 0.0);
    matrix3 liftrot = ident();
    rotate(liftrot, normalblend * (currliftangle-targetliftangle), liftaxis);
    rot = liftrot * rot;

    // build rotation around previous point
    matrix3 currtm = invert(skinspace) * rot * skinspace;

    // compute new position
    return segvec * currtm;
}


float getparmoverridevalue(string oppath; string parmname; string multiindex;
    int skininput; int skinprim; vector skinprimuv;
    vector uv;
    int curveinput; int curveprim; float curveu;
    int clumpinput; int clumpprim; float clumpu;
    int perskinpoint; int skinpoint)
{
    string override = chs(oppath + "/" + parmname + "override" + multiindex);
    if(override == "clumpattrib")
    {
	if(clumpinput == -1)
	    return 0.0;

        string attribname = chs(oppath + "/" + parmname + "clumpattrib" + multiindex);
        return select(
		hasvertexattrib(clumpinput, attribname) ||
		haspointattrib(clumpinput, attribname) ||
		hasprimattrib(clumpinput, attribname) ||
		hasdetailattrib(clumpinput, attribname),
		primuv(clumpinput, attribname, clumpprim, set(clumpu, 0, 0)),
		1.0);
    }
    if(override == "curveattrib")
    {
        string attribname = chs(oppath + "/" + parmname + "curveattrib" + multiindex);
        return select(
		hasvertexattrib(curveinput, attribname) ||
		haspointattrib(curveinput, attribname) ||
		hasprimattrib(curveinput, attribname) ||
		hasdetailattrib(curveinput, attribname),
		primuv(curveinput, attribname, curveprim, set(curveu, 0, 0)),
		1.0);
    }
    else if(override == "skinattrib")
    {
        string attribname = chs(oppath + "/" + parmname + "attrib" + multiindex);
        return select(
		perskinpoint,
		select(
		    haspointattrib(skininput, attribname),
		    point(skininput, attribname, skinpoint),
		    1.0),
		select(
		    hasvertexattrib(skininput, attribname) ||
		    haspointattrib(skininput, attribname) ||
		    hasprimattrib(skininput, attribname) ||
		    hasdetailattrib(skininput, attribname),
		    primuv(skininput, attribname, skinprim, skinprimuv),
		    1.0)
		);
    }
    else if(override == "texture")
    {
        string texture = chs(oppath + "/" + parmname + "texture" + multiindex);
	texture = expand_udim(uv.x, uv.y, texture);
        vector texval = texture(texture, uv.x, uv.y);
        return select(file_stat(texture)->isValid(), luminance(texval), 1.0);
    }
    else
    {
        return 1.0;
    }
}

float getparmoverridevalue(string oppath; string parmname; string multiindex;
	int skininput; int skinprim; vector skinprimuv;
	vector uv;
	int curveinput; int curveprim; float curveu;
	int clumpinput; int clumpprim; float clumpu)
{
    return getparmoverridevalue(oppath, parmname, multiindex,
	    skininput, skinprim, skinprimuv,
	    uv,
	    curveinput, curveprim, curveu,
	    clumpinput, clumpprim, clumpu,
	    0 /*perskinpoint*/, 0 /*skinpoint*/);
}

float getparmoverridevalue(string oppath; string parmname; string multiindex;
	int skininput; int skinprim; vector skinprimuv;
	vector uv;
	int curveinput; int curveprim; float curveu)
{
    return getparmoverridevalue(oppath, parmname, multiindex,
	    skininput, skinprim, skinprimuv,
	    uv,
	    curveinput, curveprim, curveu,
	    -1 /*clumpinput*/, -1 /*clumpprim*/, 0.0 /*clumpu*/,
	    0 /*perskinpoint*/, 0 /*skinpoint*/);
}

float evalparmoverridef(string oppath; string parmname; string multiindex;
    int skininput; int skinprim; vector skinprimuv;
    vector uv;
    int curveinput; int curveprim; float curveu;
    int clumpinput; int clumpprim; float clumpu;
    int perskinpoint; int skinpoint)
{
    float parmval = chf(oppath + "/" + parmname + multiindex);
    float overrideval = getparmoverridevalue(oppath, parmname, multiindex,
	    skininput, skinprim, skinprimuv,
	    uv,
	    curveinput, curveprim, curveu,
	    clumpinput, clumpprim, clumpu,
	    perskinpoint, skinpoint);

    return parmval * overrideval;
}

float evalparmoverridef(string oppath; string parmname; string multiindex;
    int skininput; int skinprim; vector skinprimuv;
    vector uv;
    int curveinput; int curveprim; float curveu;
    int clumpinput; int clumpprim; float clumpu)
{
    return evalparmoverridef(oppath, parmname, multiindex,
	    skininput, skinprim, skinprimuv,
	    uv,
	    curveinput, curveprim, curveu,
	    clumpinput, clumpprim, clumpu,
	    0 /*perskinpoint*/, 0 /*skinpoint*/);
}
float evalparmoverridef(string oppath; string parmname; string multiindex;
    int skininput; int skinprim; vector skinprimuv;
    vector uv;
    int curveinput; int curveprim; float curveu)
{
    return evalparmoverridef(oppath, parmname, multiindex,
	    skininput, skinprim, skinprimuv,
	    uv,
	    curveinput, curveprim, curveu,
	    -1 /*clumpinput*/, -1 /*clumpprim*/, 0.0 /*clumpu*/,
	    0 /*perskinpoint*/, 0 /*skinpoint*/);
}

#endif
