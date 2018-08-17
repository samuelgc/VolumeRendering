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
 * NAME:	advect.h ( VOP Library Functions, VEX )
 *
 * COMMENTS:	Multi-step advection functions.
 */
#ifndef __advect_h_
#define __advect_h_

#define MAX_TRACE 32

// Single Euler step.
void stepEuler(vector pos; float dt; vector vel)
{
    pos += dt * vel;
}

// Single step of Midpoint.  This midpoint formula matches the GasAdvect DOP.
void stepMidpoint(vector pos; const string input; const string vname; const float dt; const vector vel)
{
    vector midpos = pos + dt * vel;
    // Take another step of the same time value.
    vector V1 = volumesamplev(input, vname, midpos);
    vector endpos = midpos + dt * V1;
    // Average the last pos with end pos to get the midpoint
    // answer for a single step.
    pos += endpos;
    pos *= 0.5;
}

// Single step of RK3.  This RK3 formula matches vdb::tools::VelocityIntegrator
// and the VDB Advect Points SOP as of VDB 3.0.
void stepRK3(vector pos; const string input; const string vname; const float dt; const vector V0)
{
    vector V1 = volumesamplev(input, vname, pos + (0.5 * dt) * V0);
    vector V2 = volumesamplev(input, vname, pos + dt * (2 * V1 - V0));
    pos += (dt / 6) * (V0 + 4 * V1 + V2);
}

// Single step of RK4.  This RK4 formula matches vdb::tools::VelocityIntegrator
// and the VDB Advect Points SOP as of VDB 3.0.
void stepRK4(vector pos; const string input; const string vname; const float dt; const vector V0)
{
    vector V1 = volumesamplev(input, vname, pos + (0.5 * dt) * V0);
    vector V2 = volumesamplev(input, vname, pos + (0.5 * dt) * V1);
    vector V3 = volumesamplev(input, vname, pos + dt * V2);
    pos += (dt / 6) * (V0 + 2 * (V1 + V2) + V3);
}

// Advect the position in the specified volume according to the supplied
// cfl condition using the specified integration method.
void advectTrace(vector pos; const int method; const string input; const string vname; const float totaltime; const float cfl)
{
    vector		vel, invcfl;
    float		dt, timestep, multiplier, ratio, time = totaltime;
    int			iter = 0;

    // Look up volume prim for voxel size calc.
    int prim = findattribval(input, "primitive", "name", sprintf("%s.x", vname));
    if (prim < 0)
	prim = findattribval(input, "primitive", "name", vname);
    // No integration functions will work if we can't find any volumes.
    if (prim < 0)
	return;
    // Approximate voxel size.
    vector voxsize = volumeindextopos(input, prim, 1) - volumeindextopos(input, prim, 0);

    multiplier = (time < 0) ? -1 : 1;
    time *= multiplier;
    invcfl = 1.0F / (voxsize * cfl);

    while (time > 0 && iter < MAX_TRACE)
    {
    	// Assume the full timestep to start.
	timestep = time;
	// Get velocity at current position.
	vel = volumesamplev(input, vname, pos);

	// Determine our timestep to obey CFL.
	ratio = max(abs(timestep * vel * invcfl));
	if (ratio > 1.0)
	    timestep /= ratio;

	// Integrate.
	dt = timestep * multiplier;
	if (method == 1)
	    stepEuler(pos, dt, vel);
	else if (method == 2)
	    stepMidpoint(pos, input, vname, dt, vel);
	else if (method == 3)
	    stepRK3(pos, input, vname, dt, vel);
	else if (method == 4)
	    stepRK4(pos, input, vname, dt, vel);

	time -= timestep;
	iter++;
    }
}

// Advect the position in the specified volume according to the supplied
// method and CFL:
// 0 = Single Step (CFL is ignored)
// 1 = Trace (Euler)
// 2 = Trace Midpoint
// 3 = Trace RK3
// 4 = Trace RK4
vector advectbyvolumes(const string input; const string vname; const vector pos; const int method; const float dt; const float cfl)
{
    vector advectpos = pos;
    if (method == 0)
	stepEuler(advectpos, dt, volumesamplev(input, vname, pos));
    else
	advectTrace(advectpos, method, input, vname, dt, cfl);
    return advectpos;
}

// Move the position to the specified iso on the input SDF volume.
// Returns the fraction of iso-distance actually moved towards goaliso.
float movepostoiso(const int input; string vname; vector pos; const float goaliso; const float maxdist; const float tol)
{
    vector              grad;
    vector              step, voxelsize;
    float               cfl = 0.9;
    float               movedist, dist, gradlen, movetime;
    
    float               maxtime = maxdist;

    voxelsize = volumeindextopos(input, vname, 1) - volumeindextopos(input, vname, 0);
    voxelsize *= cfl;

    while (maxtime >= 0)
    {
        // Current distance value.
        dist = volumesample(input, vname, pos);
        // How far we want to move in the direction of the gradient.
        movedist = goaliso - dist;

        // See if we got to our goal.
        if (abs(movedist) < tol)
            return 1;

        grad = volumegradient(input, vname, pos);
        gradlen = length(grad);
        if (gradlen < 1e-5)
        {
            // Zero gradient leaves us no where to go, and we
            // are currently too far out to pass.
            return 0;
        }

        movetime = movedist / gradlen;
        step = grad * movetime;

        // Clamp the step by the voxel size (which includes CFL
        // condition).
        if (abs(step.x) > voxelsize.x)
        {
            movetime *= voxelsize.x / abs(step.x);
            step = grad * movetime;
        }
        if (abs(step.y) > voxelsize.y)
        {
            movetime *= voxelsize.y / abs(step.y);
            step = grad * movetime;
        }
        if (abs(step.z) > voxelsize.z)
        {
            movetime *= voxelsize.z / abs(step.z);
            step = grad * movetime;
        }

        pos += step;
        // Subtract the time used.
        maxtime -= abs(movetime);
    }
    return 0;
}

#endif
