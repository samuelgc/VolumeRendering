/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *  Side Effects Software Inc
 *  123 Front Street West
 *  Toronto, Ontario
 *  Canada   M5V 3E7
 *  416-504-9876
 *
 * NAME:    agent_util.h
 */

#ifndef __agent_util__
#define __agent_util__

#include <math.h>

#define M_TAU 6.283185307179586

///
/// Draw a set of axes at a specific point. Useful for debugging in the viewport.
///
void drawStar(vector pos) 
{
    int origin = addpoint(geoself(), pos);
    
    int prim = addprim(geoself(), "polyline");
    addvertex(geoself(), prim, origin);
    int pt = addpoint(geoself(), pos + {-0.25, 0, 0});
    addvertex(geoself(), prim, pt);
    pt = addpoint(geoself(), pos + {0.25,0,0});
    addvertex(geoself(), prim, pt);
    
    prim = addprim(geoself(), "polyline");
    addvertex(geoself(), prim, origin);
    pt = addpoint(geoself(), pos + {0, 0.25, 0});
    addvertex(geoself(), prim, pt);
    pt = addpoint(geoself(), pos + {0,-0.25,0});
    addvertex(geoself(), prim, pt);
    
    prim = addprim(geoself(), "polyline");
    addvertex(geoself(), prim, origin);
    pt = addpoint(geoself(), pos + {0, 0, 0.25});
    addvertex(geoself(), prim, pt);
    pt = addpoint(geoself(), pos + {0,0,-0.25});
    addvertex(geoself(), prim, pt);
}

///
/// Draw a line from a position going in a direction. Useful for debugging in the
/// viewport.
void drawLine(vector pos; vector dir) 
{
    int prim = addprim(geoself(), "polyline");
    int pt = addpoint(geoself(), pos);
    addvertex(geoself(), prim, pt);
    pt = addpoint(geoself(), pos + dir);
    addvertex(geoself(), prim, pt);
}

/// 
/// Description: Get the position of a matrix transform. 
///
/// Parameters:
///     * m - the matrix
///
vector getMatrixPosition(matrix m) 
{
    vector v = {0,0,0};
    v.x = m.wx;
    v.y = m.wy;
    v.z = m.wz;
    return v;
}

/// 
/// Description: Set the position of a matrix transform.
///
/// Parameters:
///     * m - the matrix
///     * v - the new position
///
void setMatrixPosition(matrix m; vector v)
{
    m.wx = v.x;
    m.wy = v.y;
    m.wz = v.z;
}

///
/// Description: Project a vector onto another.
///
/// Parameters:
///     * v1 - the first vector
///     * v2 - the second vector to project onto
/// 
vector project(const vector v1; const vector v2)
{
    return (dot(v1, v2) / length2(v2)) * v2;
}

///
/// Description: Scalar projection
///
/// Parameters:
///     * v1 - the first vector
///     * v2 - the second vector to project onto
///
float projectedDistance(const vector v1; const vector v2)
{
    return dot(v1, v2) / length(v2);
}

///
/// Description: Perpendicular vector to projection of a vector onto another
///
/// Parameters:
///     * v1 - the first vector
///     * v2 - the second vector
///
float perpproject(const vector v1; const vector v2)
{
    float projdist = projectedDistance(v1, v2);
    return sqrt(length2(v1) - projdist * projdist);
}

///
/// Description: Find the acute angle between two vectors.
///
/// Parameters:
///     * v1 - the first vector
///     * v2 - the second vector
///
float angleBetween(vector v1; vector v2) 
{
    return abs(acos(dot(v1, v2) / (length(v1) * length(v2))));
}

///
/// Rotate a matrix by a quaternion.
///
/// Parameters:
///     * m - the matrix to rotate
///     * q - the quaternion
///     * p - the pivot
///
void rotate(matrix m; vector4 q; vector p) 
{
    translate(m, -p);
    matrix3 m3 = qconvert(q);
    m = m * matrix(m3);
    translate(m, p);
}

///
/// Rotate a matrix by a quaternion.
///
/// Parameters:
///     * m - the matrix to rotate
///     * q - the quaternion
///
void rotate(matrix m; vector4 q) 
{
    rotate(m, q, {0,0,0});
}

///
/// Description: Rotate a matrix about itself by a quaternion.
///
/// Parameters:
///     * m - the matrix to rotate
///     * q - the quaternion
///
void rotateAboutSelf(matrix m; vector4 q) 
{
    rotate(m, q, getMatrixPosition(m));
}

///
/// Test whether a direction lies within a cone.
///
/// Params:
///     * conedir - the cone's centre axis
///     * up - the up vector
///     * hozlimitangle - horizontal angle of the cone 
///     * vertlimitangle - vertical angle of the cone
///     * dir - the direction to test
///
int isDirectionWithinCone(vector conedir; vector up; float hozlimitangle;
                          float vertlimitangle; vector dir)
{
    vector conedirproj = conedir - project(conedir, up);
    vector dirproj = dir - project(dir, up);
    float hozangle = angleBetween(conedirproj, dirproj);
    vector4 hozq = dihedral(conedirproj, dirproj);  

    vector conedirrotated = qrotate(hozq, conedir); 
    float vertangle = angleBetween(conedirrotated, dir); 

    return hozangle <= hozlimitangle && vertangle <= vertlimitangle;
}

///
/// Description: Find an intersection along a line. Returns the 
///              primitive number of the intersection, or -1 if there is an error.
///
/// Parameters:
///	* terrain_input_path - the path to the terrain input object
///     * origin - the origin of the ray (in world space).
///     * direction - the direction of the ray
///     * hitpos - the computed world space intersection point.
///
int
lineIntersect(const string terrain_input_path; const vector origin;
              const vector direction; const int bidirectional; vector hitpos)
{
    // Compute where the origin point is in the space of the terrain object.
    vector origin_tspace = ptransform("space:world", terrain_input_path, origin);
    vector dir_tspace = vtransform("space:world", terrain_input_path, normalize(direction));
    float max_distance = 500;
    float u;
    float v;

    int hitprimid = intersect(terrain_input_path, origin_tspace,
                              dir_tspace * max_distance, hitpos, u, v);
    if (bidirectional)
    {
        // Search in both directions, and pick the closer point.
        vector hitpos_2;
        int hitprimid_2 = intersect(terrain_input_path, origin_tspace,
                                    -dir_tspace * max_distance, hitpos_2, u, v);

        if (hitprimid_2 >= 0 &&
            (hitprimid < 0 || (distance2(origin_tspace, hitpos_2) <
                               distance2(origin_tspace, hitpos))))
        {
            hitprimid = hitprimid_2;
            hitpos = hitpos_2;
        }
    }

    // Transform the intersection point back to world space.
    hitpos = ptransform(terrain_input_path, "space:world", hitpos);
    return hitprimid;
}

///
/// Description: Solve the 2-bone IK problem. Returns the midpoint between the
///              origin and target positions. 
///
/// Parameters:
///     * origin - the origin position
///     * target - the target position
///     * l0 - the length of bone 0 
///     * l1 - the length of bone 1
///     * dir - the direction vector that, along with (target - origin), defines
///             the plane of the solution.
///
vector ik2Solve(vector origin; vector target; float l1; float l2; vector dir) 
{
    // Calculate the x and y-axes of the plane that is formed by the vector 
    // from origin to target and dir that will form a new co-ordinate system.
    vector yaxis = normalize(target - origin);
    vector plane_norm = normalize(cross(yaxis, dir));
    vector xaxis = normalize(cross(yaxis, plane_norm));
    float dist = length(target - origin);

    // Flip the sign of the xaxis if doing so makes it more similar to dir.
    if (dot(-xaxis, dir) > dot(xaxis, dir))
        xaxis = -xaxis; 
    
    vector p1 = {0,0,0};
    float theta = abs(acos((l1*l1 + dist*dist - l2*l2) / (2 * l1 * dist)));
    p1.x = l1 * sin(theta);
    p1.y = l1 * cos(theta);
    
    // Convert the point back to the original co-ordinate system.
    p1 = p1.x * xaxis + p1.y * yaxis + origin;
    return p1;
}
/// 
/// Description: Find the cross product of the image of the vectors on the XZ plane
/// 
/// Parameters:
///     * a - first vector
///     * b - second vector
float det2d(const vector a; const vector b)
{
    return a.x*b.z - a.z*b.x;
}
///
/// Description: Find the normal of the image on the XZ planei
///
/// Parameters:
///     * a - vector to find normal of  
vector norm2d(const vector a)
{
    vector b = 0;
    b.x = -a.z;
    b.z = a.x;
    return normalize(b);
}
///
/// Description: Find the closest point on the line segment to the pos
///
/// Parameters:
///     * p1 - first point on line segment 
///     * p2 - second point on line segment
///     * pos - position to find closest point to
vector closestPt2d(const vector p1; const vector p2; const vector pos)
{
    // Return p1 if line segment negligible
    float len = length(p2 - p1);
    if(len < M_TOLERANCE)
	return p1;
    
    float t = dot(pos - p1, p2 - p1) / (len*len);
    if(t < 0)
	return p1;
    else if(t > 1)
	return p2;
    
    return p1 + t * (p2 - p1);
}
///
/// Description: Find the time to collision to a line segment
///
/// Parameters:
///     * relp - displacement vector between points on line segment 
///     * vel - velocity of agent
///     * po - relative position of agent to endpoint of line segment 
///     * invdet - inverse determinant of relp and vel
///     * tau - time to collision 
float AGENTsegCollisionTime(const vector relp; const vector vel; const vector po; const float invdet; const float tau)
{
    // Use Cramer's rule to solve for line of intersection
    float t = det2d(relp, po) * invdet;
    float s = det2d(vel, po) * invdet;
    if(t > 0 && s >= 0 && s <= 1)
    {
        if(tau < 0)
            return t;
        else
            return min(tau, t);
    }
    
    return tau;
}
///
/// Description: Find the collision time between two agents.
///
/// Parameters:
///     * relp - relative position of agent relative to neighbor
///     * relv - relative velocity of agent relative to neighbor
///     * rad - agent radius
///     * rvsq - relative sqeed squared
///     * rdsq - separation distance squared
///     * rbsq - dot product of relative velocity and relative position
float AGENTcollisionTime(const vector relp; const vector relv; const float rad; export float rvsq; export float rdsq; export float rbsq) 
{
    // Find the separation distance
    float dsq = dot(relp, relp);
    rbsq = dot(relp, relv);

    // Find parameters for quadratic
    // The squared distance at time t is length2(relp + relv * t), which gives
    // dot(relp, relp) + 2 * dot(relp, relv) * t + dot(relv, relv) * t^2.
    // We want to find when the distance is (2*r)^2 == 4 * r^2
    float qa, qb, qc;
    qa = dot(relv, relv);
    qb = 2 * rbsq;
    qc = dsq - 4 * rad * rad;

    // Export relative speed and separation distance to save computation later
    rvsq = qa;
    rdsq = qc;

    // If the agents are already colliding, the time to collision is zero.
    if (qc >= 0)
    {
	// Find time to collision if one exists.
	float t1, t2;
	if (!solvequadratic(qa, qb, qc, t1, t2))
	    return -1;
	if (t1 > 0)
	{
	    if (t2 > 0)
		return min(t1, t2);
	    return t1;
	}
	if (t2 > 0)
	    return t2;
	return -1;
    }
    return 0;
}
///
/// Description: Calculate interaction force between agent and neighbor (see: http://motion.cs.umn.edu/PowerLaw/KSG-PowerLaw.pdf)
///
/// Parameters: 
///     * relp - relative position of agent relative to neighbor
///     * relv - relative velocity of agent relative to neighbor
///     * rvsq - relative sqeed squared
///     * rbsq - dot product of relative velocity and relative position
///     * rdsq - separation distance squared
///     * rad - effective radius
///     * tau - time to collision
///     * k - interaction energy
///     * tau0 - population time constant 
vector AGENTinteractionForce(const vector relp; const vector relv; const float rvsq; const float rbsq; const float rdsq; const float rad; const float k; const float tau; const float tau0)
{
    if(tau < 0)
        return 0;
    else if (tau < 0.00001f)
    {
        // If the agents are already colliding, just apply a separating force.
        return k * normalize(relp);
    }

    // The agent interaction energy is proportional to 1/tau^2 up to tau0
    float qa, qb, qc, qd;
    qa = rvsq;
    qb = -rbsq;
    qc = rdsq;
    qd = qb*qb - qa*qc;

    // Calculate the interaction force if time to collision exists and relative velocity is non-negligible
    if(qd > .0f && (qa < -0.00001f || qa > 0.00001f))
    {
        float en = k * exp(-tau / tau0) / (tau * tau);
        float cof = -(en / qa) * (2 / tau + 1 / tau0);
        return cof * (relv - ((qa * relp + qb * relv)/sqrt(qd)));
    }

    return 0;
}
///
/// Description: Find the interaction force applied to an agent due to a static obstacle (line segment)
/// (see: http://motion.cs.umn.edu/PowerLaw/KSG-PowerLaw.pdf)
///
/// Parameters:
///     * p1 - position of first point on line segment
///     * p2 - position of second point on line segment
///     * pos - position of agent
///     * vel - velocity of agent
///     * rad - effective radius of the agent
///     * k - energy constant of power law
///     * tau0 - lower bound on unscreened interactions (time to collisions >> tau0 -> E(tau) = 0)
///     * tau - time to collision
vector AGENTobstacleForce(const vector p1; const vector p2; const vector pos; const vector vel; const float rad; const float k; const float tau0; export float tau)
{
    // Vector from p1 of line segment to p2
    vector relp = p2 - p1;
    
    // Ensure that the agent is not penetrating the obstacle (perturb endpoints along normal by radius of agent)
    // Pertubation in both directions must be checked for min collision time.
    vector norm = norm2d(relp);
    vector o1 = p1 + rad * norm;
    vector o2 = p1 - rad * norm;
    
    // Vector from first point on the line segment to the agent
    vector po1 = pos - o1;
    vector po2 = pos - o2;
   
    // Solving po = t(vx + rx) + s(vy + ry) using Cramer's rule
    float detv = det2d(vel, relp);
    tau = -1;
    if(detv > 0.00001f || detv < -0.00001f) 
    {
	// Find the minimum time to collision to the line segment (same model as agentInteractionForce), directed
	// on normal to the line segment
        float invdet = 1.0f/detv;
	tau = AGENTsegCollisionTime(relp, vel, po1, invdet, tau);
	tau = AGENTsegCollisionTime(relp, vel, po2, invdet, tau);
	if(tau > 0)
	{
	    vector dir = 0;
	    dir.x = -relp.z;
	    dir.z = relp.x;
	    float en = k * exp(-tau / tau0) / (tau * tau);
	    float cof = (en / detv) * (2 / tau + 1 / tau0);
	    return cof * dir;
	} 
    }
    
    return 0;
}

/// Return the world transform of the specified bone, using the cached
/// transforms from 'world_matrices' if possible. Otherwise,
/// 'agentworldtransform' is called and the result is cached in
/// 'world_matrices'.
matrix getworldtransform(int primnum; int bone; matrix world_matrices[];
                         int valid_matrices[])
{
    if (valid_matrices[bone] == 0)
    {
        world_matrices[bone] = agentworldtransform(0, primnum, bone);
        valid_matrices[bone] = 1;
    }

    return world_matrices[bone];
}

void
AGENTcomputeChildTransforms(const int primnum; const int xform_idx;
                            const int modified_xforms[];
                            matrix world_matrices[]; int valid_matrices[])
{
    int queue[] = array(xform_idx);
    while (len(queue) > 0)
    {
        int i = pop(queue, 0);

        // Update the world transforms for any other child bones.
        if (find(modified_xforms, i) < 0)
        {
            int parent = agentrigparent(0, primnum, i);
            matrix local_xform = agentlocaltransform(0, primnum, i);
            world_matrices[i] = local_xform * world_matrices[parent];
        }

        // Update our locally cached transforms as well as the agent's
        // transforms (the calls to setagentworldtransform append the edit to
        // the geometry command queue, so we need to update our cached
        // transforms in order to do further work with the new transforms).
        valid_matrices[i] = 1;
        setagentworldtransform(geoself(), primnum, world_matrices[i], i);

        append(queue, agentrigchildren(0, primnum, i));
    }
}

///
/// Description: Run 2-bone IK and adjust a 3 joint chain (along with any
/// intermediate bones).
///
/// Parameters:
///     * primnum - the primitive number of the agent
///     * dir_2 - second direction vector that defines the plane of the 2-bone
///               IK solution.
///     * target - the target position of the ankle
///     * damping_threshold - percentage of the leg length at which knee
///                           damping is activated.
///     * world_matrices - array containing the world transforms of the bones
///     * valid_matrices - specifies which entries in 'world_matrices' are valid
///     * top_bone - index of the top bone in the chain.
///     * middle_bone - index of the middle bone in the chain.
///     * bottom_bone - index of the last bone in the chain.
///
void AGENTsolve2BoneIK(
    const int primnum;
    const vector dir_2;
    const vector target;
    const float damping_threshold;
    matrix world_matrices[];
    int valid_matrices[];
    const int top_bone;
    const int middle_bone;
    const int bottom_bone)
{
    // Get the world matrices of the bones.
    matrix W0 = getworldtransform(primnum, top_bone, world_matrices,
                                  valid_matrices);
    matrix W1 = getworldtransform(primnum, middle_bone, world_matrices,
                                  valid_matrices);
    matrix W2 = getworldtransform(primnum, bottom_bone, world_matrices,
                                  valid_matrices);

    // Get the positions of the bones.
    vector p0 = getMatrixPosition(W0);
    vector p1 = getMatrixPosition(W1);
    vector p2 = getMatrixPosition(W2);

    // Get the lengths of the bones.
    float l1 = length(p1 - p0);
    float l2 = length(p2 - p1);

    // Adjust the end affector's position if damping is enabled.
    vector damped_target = target;
    if (damping_threshold < 1.0)
    {
        float max_dist = l1 + l2;
        float damping_start = damping_threshold * max_dist;
        float d = length(target - p0);

        if (d > damping_start)
        {
            float damping_dist = max_dist - damping_start;
            float falloff = exp(-(d - damping_start) / damping_dist);
            d = damping_start + damping_dist * (1 - falloff);

            damped_target = p0 +  d * normalize(target - p0);
        }
    }

    // Solve for the new bone positions.
    vector p1p = ik2Solve(p0, damped_target, l1, l2, dir_2);

    // Update the world transforms with the new positions (subtracting out the
    // P offset).
    setMatrixPosition(W1, p1p);
    setMatrixPosition(W2, target);

    // Find the quaternion for rotating upper bone.
    // v1 = original direction of upper bone.
    // v2 = new direction of upper bone.
    vector v1 = p1 - p0;
    vector v2 = p1p - p0;
    vector4 q = dihedral(v1, v2);
    rotateAboutSelf(W0, q);

    v1 = p2 - p1;
    v2 = target - p1p;
    q = dihedral(v1, v2);
    rotateAboutSelf(W1, q);

    // Update the world transforms of the bones.
    world_matrices[top_bone] = W0;
    world_matrices[middle_bone] = W1;
    world_matrices[bottom_bone] = W2;
    int modified_xforms[] = array(top_bone, middle_bone, bottom_bone);

    AGENTcomputeChildTransforms(primnum, top_bone, modified_xforms,
                                world_matrices, valid_matrices);
}

///
/// Description: Tilt the back of the agent depending on the slope of the terrain.
///
/// Parameters:
///     * ref_up - the up direction in the local space of the agent
///     * ref_dir - the direction the agent faces in local space
///     * tnorm - the terrain normal in local space
///     * tilt - the current tilt direction
///     * tilt_angle - the angle to tilt per frame
///     * tilt_min - the minimum tilt angle
///     * tilt_max - the maximum tilt angle
///     * lowerback - index of the lower back bone
///
vector AGENTtiltBack(
    const int primnum;
    const vector ref_up;
    const vector ref_dir;
    const vector tnorm;
    const vector tilt;
    const float tilt_angle;
    const float tilt_min;
    const float tilt_max;
    const int lowerback)
{
    float theta = min(angleBetween(ref_up, tnorm), PI_2);
    // The goal tilt angle follows a quadratic function where omega = 0 for
    // theta = 0 and theta = pi/2, and the maximum is omega = pi/8 at theta =
    // pi/4. This ensures that the agent doesn't tilt by unrealistic angles.
    float omega = theta - (theta * theta) / PI_2;

    if (dot(tnorm, ref_dir) >= 0)
        omega *= -1;

    omega = clamp(omega, tilt_min, tilt_max);

    vector norm = normalize(cross(tnorm, ref_dir));

    vector4 q = quaternion(omega, norm);
    vector targettilt = qrotate(q, ref_up);

    q = quaternion(min(tilt_angle, angleBetween(tilt, targettilt)),
                   normalize(cross(tilt, targettilt)));
    vector newtilt = qrotate(q, tilt);

    q = dihedral(ref_up, newtilt);

    matrix W = agentworldtransform(0, primnum, lowerback);
    rotateAboutSelf(W, q);

    int parent = agentrigparent(0, primnum, lowerback);
    matrix W_parent = agentworldtransform(0, primnum, parent);
    setagentlocaltransform(geoself(), primnum, W * invert(W_parent), lowerback);

    return newtilt;
}

/// Choose one of a set of possible look at targets for an agent.
void agentchooselookat(
    const int primnum;
    const vector up;
    const string head_bone;
    const float horiz_limit_angle;
    const float vert_limit_angle;
    const vector targets[];
    const int enabled_targets[];
    int current_target_idx;
    vector current_target)
{
    if (len(targets) == 0)
    {
        current_target_idx = -1;
        return;
    }

    matrix agent_xform = primintrinsic(0, "packedfulltransform", primnum);
    matrix agent_xform_inv = invert(agent_xform);

    matrix Wcurr = agentworldtransform(0, primnum, agentrigfind(0, primnum, head_bone));
    vector headpos = cracktransform(XFORM_SRT, XFORM_XYZ, 0, {0, 0, 0}, Wcurr);

    vector zaxis = {0, 0, 1};

    // Convert to the local space of the agent.
    vector targets_local[];
    foreach (vector target; targets)
        append(targets_local, target * agent_xform_inv);

    // Check if existing target is still in range.
    if (current_target_idx >= 0 && current_target_idx < len(targets) &&
        enabled_targets[current_target_idx] &&
        isDirectionWithinCone(zaxis, up, horiz_limit_angle, vert_limit_angle,
                              targets_local[current_target_idx] - headpos))
    {
        current_target = targets[current_target_idx];
    }
    else
    {
        current_target_idx = -1;
        current_target = {0, 0, 0};
    }

    // Find a new target if necessary.
    if (current_target_idx < 0)
    {
        float min_dist = -1;

        foreach (int i; vector target; targets)
        {
            if (enabled_targets[i] &&
                isDirectionWithinCone(zaxis, up, horiz_limit_angle,
                                      vert_limit_angle,
                                      targets_local[i] - headpos))
            {
                float dist = length(targets_local[i] - headpos);
                if (min_dist < 0 || dist < min_dist)
                {
                    min_dist = dist;
                    current_target_idx = i;
                    current_target = targets[i];
                }
            }
        }
    }
}

/// Modify an agent's head bone to look at a target position.
/// If the agent doesn't have a current target, it will gradually move back to
/// the head position from its current clip.
void agentlookattarget(
    const int primnum;
    const string head_bone;
    const string rest_clip;

    const int has_target1;
    vector target1;
    const int has_target2;
    vector target2;
    const float weight; // Optionally, blend between two target positions.

    const vector eye_offset;
    const int limit_head_turn;
    const float max_turn_angle;
    vector current_head_dir;
    )
{
    matrix agent_xform = primintrinsic(0, "packedfulltransform", primnum);
    matrix agent_xform_inv = invert(agent_xform);
    vector zaxis = {0, 0, 1};

    int headidx = agentrigfind(0, primnum, head_bone);
    matrix Wcurr = agentworldtransform(0, primnum, headidx);
    matrix Wparent = agentworldtransform(0, primnum, agentrigparent(0, primnum, headidx));
    vector headpos = getMatrixPosition(Wcurr) + eye_offset;

    if (find(agentclipcatalog(0, primnum), rest_clip) < 0)
        warning("Rest clip '%g' does not exist.", rest_clip);
    matrix Wrest = agentclipsampleworld(0, primnum, rest_clip, 0.0, headidx);

    vector headscale = cracktransform(XFORM_SRT, XFORM_XYZ, 2, {0,0,0}, Wcurr);
    vector headrotrest = radians(cracktransform(XFORM_SRT, XFORM_XYZ, 1, {0,0,0}, Wrest));
    vector headtrans = cracktransform(XFORM_SRT, XFORM_XYZ, 0, {0,0,0}, Wcurr); 

    Wrest = ident();
    scale(Wrest, headscale); 
    Wrest *= matrix(qconvert(eulertoquaternion(headrotrest, XFORM_XYZ)));
    translate(Wrest, headtrans);

    // Get the direction the agent is looking in its current clip.
    matrix Wrest_inv = invert(Wrest);
    vector localheadoffset = headpos * Wrest_inv;
    vector localzaxis = (headpos + zaxis) * Wrest_inv;
    vector clip_dir = localzaxis * Wcurr - localheadoffset * Wcurr;
    vector clip_target = (headpos + clip_dir) * agent_xform;

    // The first time around, our current head direction is whatever the
    // current clip uses.
    if (current_head_dir == {0, 0, 0})
        current_head_dir = clip_dir;

    // Find the target position.
    vector target;
    if (!has_target1 && !has_target2)
    {
        // If we don't have a target, move the head back to follow the clip.
        target = clip_target;
    }
    else
    {
        // Blend between the given targets.
        if (!has_target1)
            target1 = clip_target;
        if (!has_target2)
            target2 = clip_target;

        target = target1 + weight * (target2 - target1);
    }

    // Transform the target to be relative to the agent's default orientation.  
    target *= agent_xform_inv;

    vector target_dir = target - getMatrixPosition(Wcurr);
    vector target_lookat_dir;
    if (eye_offset == {0, 0, 0})
    {
        target_lookat_dir = target_dir;
    }
    else
    {
        float l = length(eye_offset);
        float d = length(target_dir);
        float theta = PI - angleBetween(eye_offset, zaxis);

        float t1, t2;
        solvequadratic(1, -2*l*cos(theta), l*l - d*d, t1, t2);
        float x = abs(max(t1, t2));
        float phi = acos((l*l + d*d - x*x) / (2*l*d));

        vector headoffsetaxis = normalize(cross(zaxis, eye_offset));
        target_lookat_dir = qrotate(quaternion(phi, headoffsetaxis), target_dir);
        vector lookataxis = normalize(cross(eye_offset, zaxis));
        target_lookat_dir = qrotate(quaternion(PI - theta, lookataxis), target_lookat_dir);
    }

    vector new_lookat_dir = {0, 0, 0};
    if (limit_head_turn == 0)
    {
        new_lookat_dir = target_lookat_dir;
    }
    else
    {
        vector rotaxis = normalize(cross(current_head_dir, target_lookat_dir));
        float rotangle = min(max_turn_angle, angleBetween(current_head_dir, target_lookat_dir));
        new_lookat_dir = qrotate(quaternion(rotangle, rotaxis), current_head_dir);
    }

    matrix W = Wrest;
    rotateAboutSelf(W, dihedral(zaxis, new_lookat_dir));

    matrix local_xform = W * invert(Wparent);
    setagentlocaltransform(geoself(), primnum, local_xform, headidx);

    current_head_dir = new_lookat_dir;  
}
/// 
/// Description: Obstacle avoidance force between a static point and an agent
///
/// Parameters: 
///     * pos - object space position of the agent
///     * vel - object space velocity of the agent
///     * rad - radius of the agent
///     * k - energy constant
///     * tau0 - upper bound for collision time 
///     * tau - collision time
vector AGENTobstaclePtForce(const vector pos; const vector vel; const vector hitpos; const float rad; const float k; const float tau0; export float tau)
{
    vector wPos = hitpos;
    float rvsq, rdsq, rbsq;
    float ctime = AGENTcollisionTime(pos - wPos, vel, rad, rvsq, rdsq, rbsq);
    vector force = 0;
    tau = -1;
    if(ctime > 0) 
    {
        // Consider the image of the wPos
        tau = ctime;
        force = AGENTinteractionForce(pos - wPos, vel, rvsq, rbsq, rdsq, rad, k, ctime, tau0);  
    }

    return force;
}
///
/// Description: Obstacle avoidance force for general geometry
/// 
/// Parameters:
///     * pos - position of the agent (world space)
///     * vel - velocity of the agent (world space)
///     * up - up vector of the agent (world space)
///     * samples - number of rays to cast
///     * geo - input geometry
///     * maxdist - maximum search distance for an intersection
///     * fovhorizontal - fovhorizontal
///     * fovvertical - fovvertical
///     * seed - seed for random sampling
///     * rad - radius of the agent
///     * k - energy constant 
///     * tau0 - upper bound for collision time
///     * sampleweightbias - sample weight bias for long range obstacle avoidance
///     * linforcescale - avoidance force scale for long range obstacle avoidance
vector 
AGENTobstacleForceEx(
	const vector pos; 
	const vector vel; 
	const vector up; 
	const int samples; 
	const string geo; 
	const float maxdist; 
	const float fovhorizontal; 
	const float fovvertical; 
	const float seed; 
	const float rad;
	const float k; 
	const float tau0; 
	float samplesweightbias; 
	float linforcescale;
	vector ringbuf_hitpos[];
	vector ringbuf_dir[];
	int    ringbuf_id[];
	int    offset) 
{   
    // Define frustum about the velocity/heading 
    vector force = 0;
    vector objectP = ptransform("space:world", geo, pos);  
    vector objectV = vtransform("space:world", geo, vel);
    vector objectU = vtransform("space:world", geo, up);
    vector optimaldirection = 0;
    int hitcount = 0;
    float avghitdist = 0;
    float minhitdist = maxdist;
    
    // Compute transformation matrix for each sample (first need to orient to z-axis)
    matrix tform = lookat(objectV, {0, 0, 0}, objectU);
    matrix zform = 1;
    rotate(zform, M_PI_2, {0, 1, 0});
    rotate(zform, M_PI_2, {0, 0, 1});
    tform = zform * tform;

    // Compute the parameterizations over the sphericial rectangle
    float theta = (fovhorizontal/360) * M_PI;
    float phi = (fovvertical/360) * M_PI;
    float sphi = sin(phi);
    theta /= M_TAU;

    // NB: Amortize the cost of ray casting over time (keep a ring buffer that 
    // gets pushed with 1/10 samples)
    int sub_samples = samples;

    // Initialize the ring buffer on first frame
    if (samples > 10 && len(ringbuf_hitpos) > 0)
	sub_samples /= 10;

    for(int i = 0; i < sub_samples; ++i)
    {
	float a = rand(2*i + seed + offset);
        float b = rand(2*i + 1 + seed + offset);
	vector2 uv = set(0.5 +     (2*a-1)*theta,
			 0.5 + 0.5*(2*b-1)*sphi);
        vector sampledir = sample_direction_uniform(uv) * tform;

        /// Intersect sample    
        vector hitpos = 0;
        float u, v;
        int hitid = intersect(geo, objectP, sampledir * maxdist, hitpos, u, v);
	
	push(ringbuf_hitpos, hitpos);
	push(ringbuf_dir, sampledir);
	push(ringbuf_id, hitid);
	if (len(ringbuf_hitpos) > samples)
	{
	    pop(ringbuf_hitpos, 0);
	    pop(ringbuf_dir, 0);
	    pop(ringbuf_id, 0);
	}
    }

    offset = (offset + sub_samples) % samples;
          
    // Uniform random sampling about the velocity based on the FOV parameters
    for(int i = 0; i < len(ringbuf_hitpos); ++i) 
    {
	int hitid = ringbuf_id[i];
	vector sampledir = ringbuf_dir[i];
	vector hitpos = ringbuf_hitpos[i];
    
        if(hitid != -1) 
        {
            // Determine short range power law behavior (~= 0 for tau > tau0)
            float tau;
            force += AGENTobstaclePtForce(objectP, objectV, hitpos, rad, k, tau0, tau);
            
            // Determine long range best path behavior
            hitcount++;
            float hitdistance = length(hitpos - objectP);
            float hitdistanceratio = hitdistance / maxdist;
            avghitdist += hitdistance;
            if (hitdistance < minhitdist) 
                minhitdist = hitdistance;

            // Longer path without collisions
            optimaldirection += sampledir * pow(hitdistanceratio, samplesweightbias); 
        }
        else  
            optimaldirection += sampledir * samplesweightbias;
               
    }

    force = vtransform(geo, "space:world", force);

    if (hitcount > 0)
    {
        optimaldirection =  vtransform(geo, "space:world", normalize(optimaldirection));
        vector linforce = optimaldirection;
        float hitfalloff = 1;

        avghitdist /= hitcount;
        hitfalloff = (1 - avghitdist/maxdist);

        // Combine total force from linear avoidance force and power law avoidance force
        force += optimaldirection * linforcescale * hitfalloff;
    }

    return force;
}

float[] 
AGENTgetTimeToCollisionFOV(const float tau0; const float fov; const int samples; const vector pos; const vector vel; const vector up; const float rad; const string geo; const float maxdist)
{
    float tau_bins[];
    resize(tau_bins, samples);
    
    // Transform to object space
    vector sampledir = normalize(vel);
    vector objectP = ptransform("space:world", geo, pos);  
    vector objectV = vtransform("space:world", geo, vel);
    vector objectU = vtransform("space:world", geo, up);
    
    for (int i = 0; i < samples; ++i)
    {
	matrix3 tform = 1;
	float rvsq, rdsq, rbsq, ctime, u, v;
	vector hitpos;
	int hitid;
	float theta = fit(i, 0, samples - 1, -fov/360 * PI, fov/360 * PI);

	// Fill in each bin
	rotate(tform, theta, objectU);
	sampledir = normalize(objectV) * tform;
	hitid = intersect(geo, objectP, sampledir * maxdist, hitpos, u, v);
	if (hitid != -1)
	    tau_bins[i] = AGENTcollisionTime(objectP - hitpos, objectV, rad, rvsq, rdsq, rbsq);
	else 
	    tau_bins[i] = -1;

	if (tau_bins[i] > tau0 || tau_bins[i] < 0)
	    tau_bins[i] = tau0;
    }

    return tau_bins;
}

float
agent_integratespring(const float stiffness; const float damping;
                  const float max_accel; const float timestep; const float x_0;
                  const float v_0)
{
    float h = timestep;
    float h_2 = h * h;
    float k = stiffness;
    float c = damping;

    // Acceleration follows a(x, v) = -c * v - k * x
    // Obtain the new velocity and position using a backwards Euler step:
    // x_1 = x_0 + h * v_1
    // v_1 = v_0 + h * a(x_1, v_2)
    float v_1 = (v_0 - h * k * x_0) / (1 + c * h + h_2 * k);

    // Optionally limit maximum acceleration.
    if (max_accel >= 0)
    {
        float d_v = v_1 - v_0;
        d_v = min(max_accel * timestep, abs(d_v)) * sign(d_v);
        v_1 = v_0 + d_v;
    }

    float d_x = v_1 * h;
    return d_x;
}

#endif
