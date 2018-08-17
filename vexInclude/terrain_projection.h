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
 * NAME:    terrain_projection.h
 */

#ifndef __terrain_projection__
#define __terrain_projection__

struct TerrainProjectParms
{
    string terrain_path; /// Path to the terrain geometry.
    int mode; /// 0 = specified direction, 1 = up attribute.
    vector projection_dir;
    float offset;
    int project_forces;
};

vector
terrain_getprojectiondir(const TerrainProjectParms parms; const vector up)
{
    vector projection_dir;
    if (parms.mode == 0)
        projection_dir = parms.projection_dir;
    else
        projection_dir = up;

    projection_dir = normalize(projection_dir);
    return projection_dir;
}

void
terrain_projectparticle(const TerrainProjectParms parms; const vector up;
                        vector P; vector terrain_normal; vector terrain_uvw;
                        int terrain_prim; vector terrain_offset; vector force)
{
    // Choose projection direction, and transform into the terrain object's
    // space.
    vector projection_dir = terrain_getprojectiondir(parms, up);
    projection_dir =
        vtransform("space:world", parms.terrain_path, projection_dir);

    vector hitpos;
    vector hituvw;
    int hitprim = -1;
    vector object_P = ptransform("space:world", parms.terrain_path, P);

    hitprim = intersect(parms.terrain_path, object_P, -projection_dir * 1e8,
                        hitpos, hituvw);
    {
        // Search in both directions, and take the closer point.
        vector hitpos_2;
        vector hituvw_2;
        int hitprim_2 = -1;

        hitprim_2 = intersect(parms.terrain_path, object_P,
                              projection_dir * 1e8, hitpos_2, hituvw_2);

        if (hitprim_2 >= 0 && (hitprim < 0 || (distance2(object_P, hitpos_2) <
                                               distance2(object_P, hitpos))))
        {
            hitprim = hitprim_2;
            hitpos = hitpos_2;
            hituvw = hituvw_2;
        }
    }

    // Record primitive and UV coordinates at the intersection point.
    if (hitprim >= 0)
    {
        terrain_offset = projection_dir * -1 * parms.offset;

        hitpos += terrain_offset;
        P = ptransform(parms.terrain_path, "space:world", hitpos);

        terrain_normal = ntransform(
            parms.terrain_path, "space:world",
            normalize(prim_normal(parms.terrain_path, hitprim, hituvw)));

        terrain_uvw = hituvw;
        terrain_prim = hitprim;

        // Project forces if necessary.
        if (parms.project_forces)
        {
            vector primtangent = normalize(
                cross(cross(normalize(force), terrain_normal), terrain_normal));
            force = dot(force, primtangent) * primtangent;
        }
    }
    else
    {
        terrain_prim = -1;
        terrain_uvw = {0, 0, 0};
        terrain_offset = {0, 0, 0};
    }
}


/// Simpler interface that doesn't record extra information about the
/// intersection, such as UV coordinates.
void
terrain_projectparticle(const TerrainProjectParms parms; const vector up;
                        vector P; vector terrain_normal)
{
    vector terrain_uvw, terrain_offset, force;
    int terrain_prim;

    terrain_projectparticle(
        parms, up, P, terrain_normal, terrain_uvw, terrain_prim, terrain_offset,
        force);
}

/// Returns the rotation needed to align the up vector with the terrain normal.
vector4
terrain_rotateupvector(const vector up; const vector terrain_normal)
{
    vector new_up = terrain_normal;

    // Terrain can be double-sided, so stay on the same side as the current up
    // vector.
    if (dot(up, terrain_normal) < 0)
        new_up *= -1;

    vector4 r = dihedral(up, new_up);
    return r;
}

#endif
