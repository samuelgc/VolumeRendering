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
 * NAME:    ragdoll.h
 */

#ifndef __ragdoll_h__
#define __ragdoll_h__

#include <math.h>

vector4
extractorient(const matrix m)
{
    vector r = cracktransform(XFORM_SRT, XFORM_XYZ, 1, {0, 0, 0}, m);
    return eulertoquaternion(radians(r), XFORM_XYZ);
}

vector
gettranslation(const matrix m)
{
    return set(m.wx, m.wy, m.wz);
}

int
equalzero(const float f; const float tol)
{
    return f >= -tol && f <= tol;
}

/// Find the closest joint that has a shape attached to it (for ragdolls,
/// intermediate or leaf joints remain fixed to their parent).
int
ragdoll_findparentshape(const int input; const int primnum;
                        const string layername; const int xform_id)
{
    // Find the closest parent that has a collision shape in this layer.
    int parent_shape_xform_id = xform_id;
    while (parent_shape_xform_id >= 0 &&
           !len(agentlayershapes(input, primnum, layername,
                                 parent_shape_xform_id)))
    {
        parent_shape_xform_id =
            agentrigparent(input, primnum, parent_shape_xform_id);
    }

    return parent_shape_xform_id;
}

/// Compute the rotation to align the y-axis with the default twist axis for
/// cone twist constraints (i.e. along the parent bone).
vector4
ragdoll_establishtwistaxis(const int primnum; const int xform_idx)
{
    vector twist_dir = {0, 1, 0};
    int i = xform_idx;

    while (i >= 0)
    {
        vector trans = gettranslation(agentlocaltransform(0, primnum, i));
        if (!equalzero(length2(trans), 0.00001))
        {
            twist_dir = normalize(trans);
            break;
        }
        else
            i = agentrigparent(0, primnum, i);
    }

    return dihedral({0, 1, 0}, twist_dir);
}

/// Helper function to set up a constraint between two shapes.
int
ragdoll_buildconstraint(const int primnum; const vector anchor;
                        const string agent_name; const string xform_names[];
                        const int xform_idx_0; const matrix xform0;
                        const int xform_idx_1; const matrix xform1;
                        const int is_pin)
{
    int point_0 = addpoint(0, anchor);
    string anchor_name_0 =
        sprintf("%s/%s", agent_name, xform_names[xform_idx_0]);
    setpointattrib(0, "name", point_0, anchor_name_0);
    setpointattrib(0, "anchor_type", point_0, "agent");
    setpointattrib(0, "local_P", point_0, anchor * invert(xform0));

    // For pin constraints, set the anchor's orientation so that the initial
    // relative orientation is maintained.
    vector4 r0 = extractorient(xform0);
    vector4 r1 = extractorient(xform1);
    vector4 offset = qmultiply(qinvert(r0), r1);
    setpointattrib(0, "local_orient", point_0, is_pin ? offset : {0, 0, 0, 1});

    int point_1 = addpoint(0, anchor);
    string anchor_name_1 =
        sprintf("%s/%s", agent_name, xform_names[xform_idx_1]);
    setpointattrib(0, "name", point_1, anchor_name_1);
    setpointattrib(0, "anchor_type", point_1, "agent");
    setpointattrib(0, "local_P", point_1, anchor * invert(xform1));
    setpointattrib(0, "local_orient", point_1, {0, 0, 0, 1});

    int prim = addprim(0, "polyline");
    addvertex(0, prim, point_0);
    addvertex(0, prim, point_1);

    // Record the source agent's name to make it easier to match up constraints
    // with their agent.
    setprimattrib(0, "name", prim, agent_name);

    // For pin constraints we want a rest length of zero, and cone twist
    // constraints don't do anything with that attribute (and therefore behave
    // as though the rest length was 0).
    setprimattrib(0, "restlength", prim, 0.0);

    return prim;
}

void
ragdoll_buildconstraintnetwork(const int ptnum; const int primnum;
                               const int pin_root_shapes;
                               const int pin_unconfigured_shapes)
{
    if (primintrinsic(0, "typename", primnum) != "PackedAgent")
        return;

    string required_attribs[] = {
        "name", "agent_jointgoalxforms", "agent_jointparents",
        "agent_jointlimits", "agent_jointconstrainedxforms"
    };
    foreach(string attrib; required_attribs)
    {
        if (!haspointattrib(0, attrib))
        {
            warning(
                "Point attribute '%s' not found - no constraints were created. "
                "This attribute can be set up by the Agent Configure Joints "
                "SOP.",
                attrib);

            return;
        }
    }

    string agent_name = point(0, "name", ptnum);
    matrix goal_xforms[] = point(0, "agent_jointgoalxforms", ptnum);
    vector joint_limits[] = point(0, "agent_jointlimits", ptnum);
    vector4 constrained_xforms[] =
        point(0, "agent_jointconstrainedxforms", ptnum);
    int joint_parents[] = point(0, "agent_jointparents", ptnum);

    // Create the primitive attributes for the constraint network.
    addpointattrib(0, "anchor_type", "");
    addpointattrib(0, "local_P", {0, 0, 0});
    addpointattrib(0, "local_orient", {0, 0, 0, 1});
    addprimattrib(0, "constraint_name", "");
    addprimattrib(0, "constraint_type", "");
    addprimattrib(0, "constrained_twist_axis", {0, 0, 0});
    addprimattrib(0, "goal_twist_axis", {0, 0, 0});
    addprimattrib(0, "constrained_up_axis", {0, 0, 0});
    addprimattrib(0, "goal_up_axis", {0, 0, 0});
    addprimattrib(0, "max_up_rotation", 0.0);
    addprimattrib(0, "max_twist", 0.0);
    addprimattrib(0, "max_out_rotation", 0.0);
    addprimattrib(0, "name", "");
    addprimattrib(0, "restlength", 0.0);

    matrix agent_xform = primintrinsic(0, "packedfulltransform", primnum);
    string layername = agentcollisionlayer(0, primnum);

    // TODO - avoid constructing duplicate constraints if multiple shapes are
    // bound to the same transform.
    int xform_ids[] = agentlayerbindings(0, primnum, layername, "all");
    string xform_names[] = agenttransformnames(0, primnum);
    int root_ids[] = {};

    foreach (int xform_id; xform_ids)
    {
        int parent_xform_id = joint_parents[xform_id];

        // Skip transforms that are root nodes or that don't have parent with
        // an attached shape.
        if (ragdoll_findparentshape(0, primnum, layername, parent_xform_id) < 0)
        {
            append(root_ids, xform_id);
            continue;
        }

        matrix current_xform =
            agentworldtransform(0, primnum, xform_id) * agent_xform;
        matrix parent_xform =
            agentworldtransform(0, primnum, parent_xform_id) * agent_xform;
        matrix goal_xform_world = goal_xforms[xform_id] * parent_xform;

        vector limits = joint_limits[xform_id];
        int has_cone = (limits != set(0, 0, 0));

        if (!has_cone && !pin_unconfigured_shapes)
            continue;

        vector anchor_pt =
            gettranslation(has_cone ? goal_xform_world : current_xform);
        int prim = ragdoll_buildconstraint(
            primnum, anchor_pt, agent_name, xform_names, xform_id,
            current_xform, parent_xform_id, parent_xform, !has_cone);

        if (has_cone)
        {
            vector twist_axis = {0, 1, 0};
            vector up_axis = {1, 0, 0};
            matrix3 goal_xform3 = matrix3(goal_xforms[xform_id]);

            // Find the orientation of the twist/up axis in the constrained
            // (child) object's space.
            vector4 constrained_xform = constrained_xforms[xform_id];
            vector current_twist_dir =
                normalize(qrotate(constrained_xform, twist_axis));
            vector current_up_dir =
                normalize(qrotate(constrained_xform, up_axis));

            // Set up the cone twist constraint.
            setprimattrib(0, "constraint_name", prim, "ConeTwist");
            setprimattrib(0, "constraint_type", prim, "rotation");

            setprimattrib(0, "constrained_twist_axis", prim, current_twist_dir);
            setprimattrib(0, "goal_twist_axis", prim,
                          twist_axis * goal_xform3);

            setprimattrib(0, "constrained_up_axis", prim, current_up_dir);
            setprimattrib(0, "goal_up_axis", prim, up_axis * goal_xform3);

            setprimattrib(0, "max_up_rotation", prim, limits.x);
            setprimattrib(0, "max_twist", prim, limits.y);
            setprimattrib(0, "max_out_rotation", prim, limits.z);
        }
        else if (pin_unconfigured_shapes)
        {
            // Create pin constraints if the joint limit was not set.
            setprimattrib(0, "constraint_name", prim, "Pin");
            setprimattrib(0, "constraint_type", prim, "all");
        }
    }

    if (pin_root_shapes)
    {
        // If there are multiple root shapes (i.e. there are multiple connected
        // components in the layer), join them with hard constraints so that
        // the ragdoll doesn't fall apart.
        for (int i = 1; i < len(root_ids); ++i)
        {
            int xform_0 = root_ids[i - 1];
            int xform_1 = root_ids[i];

            vector P = gettranslation(agentworldtransform(0, primnum, xform_0));
            P *= agent_xform;

            matrix child_xform =
                agentworldtransform(0, primnum, xform_0) * agent_xform;
            matrix parent_xform =
                agentworldtransform(0, primnum, xform_1) * agent_xform;

            int prim = ragdoll_buildconstraint(
                primnum, P, agent_name, xform_names, xform_0, child_xform,
                xform_1, parent_xform, 1);
            setprimattrib(0, "constraint_name", prim, "Pin");
            setprimattrib(0, "constraint_type", prim, "all");
        }
    }
}

/// Convert from the Agent Configure Joints node's parameters to the format
/// used by Bullet's cone twist constraints.
void
ragdoll_getjointinfo(const int primnum; const int xform_idx;
                     const vector offset; const vector rot;
                     const vector2 twist_range; const vector2 up_range;
                     const vector2 out_range; export float twist_angle;
                     export float up_angle; export float out_angle;
                     export matrix xform)
{
    vector twist_axis = {0, 1, 0};
    vector up_axis = {1, 0, 0};
    vector out_axis = {0, 0, 1};

    // Initially align the twist axis along the direction from parent to child.
    vector4 initial_alignment = ragdoll_establishtwistaxis(primnum, xform_idx);

    // Convert the rotation limit parameters into the format that Bullet
    // expects. The cone twist constraint uses the maximum rotation around the
    // up/out axis in either direction instead of a min/max rotation range, so
    // we need to figure out where the twist axis should be.
    vector4 establish_twist = set(0, 0, 0, 1);
    up_angle = 0;
    out_angle = 0;
    twist_angle = 0;

    if (up_range.y >= up_range.x)
    {
        // If the span of the cone exceeds 2 * M_PI, we proceed as if
        // up_range.y = up_range.x + 2 * M_PI.
        if(up_range.y - up_range.x > 2.0 * M_PI)
        {
            up_angle = M_PI;
            establish_twist = quaternion(up_range.x + M_PI, up_axis);
        }
        else
        {
            up_angle = 0.5 * (up_range.y - up_range.x);
            establish_twist = quaternion(0.5 * (up_range.x + up_range.y),
					 up_axis);
        }
    }

    if (out_range.y >= out_range.x)
    {
        // If the span of the cone exceeds 2 * M_PI, we proceed as if
        // out_range.y = out_range.x + 2 * M_PI.
        if(out_range.y - out_range.x > 2.0 * M_PI)
        {
            out_angle = M_PI;
            establish_twist = qmultiply(establish_twist,
                              quaternion(out_range.x + M_PI, out_axis));
        }
        else
        {   
            out_angle = 0.5 * (out_range.y - out_range.x);
            establish_twist = qmultiply(establish_twist,
                              quaternion(0.5 * (out_range.x + out_range.y),
					 out_axis));
        }
    }

    if (twist_range.y >= twist_range.x)
    {
        twist_angle = 0.5 * (twist_range.y - twist_range.x);

        establish_twist = qmultiply(
            establish_twist,
            quaternion(0.5 * (twist_range.x + twist_range.y), twist_axis));
    }

    matrix align_joints =
        qconvert(qmultiply(initial_alignment, establish_twist));

    // Apply the rotation and translation parameters.
    matrix local_xform =
        maketransform(XFORM_SRT, XFORM_XYZ, offset, rot, {1, 1, 1});

    xform = align_joints * local_xform;
}

/// Adds a line from 'origin_pt' to the given position.
void
ragdoll_drawline(const int origin_pt; const vector pt; const vector color)
{
    int prim = addprim(0, "polyline");
    addvertex(0, prim, origin_pt);
    addvertex(0, prim, addpoint(0, pt));

    setprimattrib(0, "Cd", prim, color);
}

/// Draw guide geometry for the cone (see SIM_ConRelConeTwist).
void
ragdoll_drawcone(const int origin_ptnum; const vector origin; float up_range;
                 float out_range; const matrix3 xform; const float guide_scale;
                 const vector color)
{
    int nsamples = 360;
    int prim = addprim(0, "polyline");
    setprimattrib(0, "Cd", prim, color);

    for (int i = 0; i <= nsamples; ++i)
    {
        float angle = 2.0 * M_PI * (i / (float)nsamples);

        float x = up_range * cos(angle);
        float z = out_range * sin(angle);
        vector axis = set(x, 0, z);
        float range = length(axis);

        vector4 rot = quaternion(range, normalize(axis));
        vector dir = normalize(qrotate(rot, {0, 1, 0}) * xform);
        vector p = origin + dir * guide_scale;

        int pt = addpoint(0, p);
        addvertex(0, prim, pt);
        
        // Draw links to origin.
        if (i % (nsamples / 4) == 0)
        {
            addvertex(0, prim, origin_ptnum);
            addvertex(0, prim, pt);
        }
    }
}

vector4
ragdoll_getchildorient(const int primnum; const int xform_idx;
                       const vector child_rot)
{
    return qmultiply(eulertoquaternion(radians(child_rot), XFORM_XYZ),
                     ragdoll_establishtwistaxis(primnum, xform_idx));
}

/// Draw guide geometry for a single joint.
void
ragdoll_drawjointguide(const int primnum; const int xform_idx;
                       const int parent_xform_idx; const vector offset;
                       const vector rot; const vector child_rot;
                       const vector2 twist_range; const vector2 up_range;
                       const vector2 out_range; const float scale;
                       const vector cone_guide_color;
                       const vector twist_limit_guide_color;
                       const vector twist_axis_color;
                       const vector up_axis_color)
{
    vector twist_axis = {0, 1, 0};
    vector up_axis = {1, 0, 0};
    float twist_limit_angle;
    float up_limit_angle;
    float out_limit_angle;
    matrix xform;

    ragdoll_getjointinfo(primnum, xform_idx, offset, rot, twist_range, up_range,
                         out_range, twist_limit_angle, up_limit_angle,
                         out_limit_angle, xform);

    matrix agent_xform = primintrinsic(0, "packedfulltransform", primnum);
    matrix parent_xform;
    if (parent_xform_idx < 0)
        parent_xform = agent_xform;
    else
    {
        parent_xform =
            agentworldtransform(0, primnum, parent_xform_idx) * agent_xform;
    }

    xform *= parent_xform;

    // Size the cones to initially be half of the bone length, and then
    // multiply by the user's guide scale.
    int children[] = agentrigchildren(0, primnum, xform_idx);
    matrix current_xform =
        agentworldtransform(0, primnum, xform_idx) * agent_xform;
    vector current_pos = gettranslation(current_xform);
    float guide_scale = 0;
    foreach(int child; children)
    {
        vector delta = gettranslation(agentworldtransform(0, primnum, child) *
                                      agent_xform) - current_pos;
        guide_scale = max(guide_scale, 0.5 * length(delta));
    }

    // Ensure that the guides are still visible for very small bones.
    guide_scale = scale * max(guide_scale, 0.09);

    addprimattrib(0, "Cd", {0, 0, 0}, "color");
    vector anchor_pos = gettranslation(xform);
    int anchor_pt = addpoint(0, anchor_pos);

    vector4 child_orient = ragdoll_getchildorient(primnum, xform_idx,
						  child_rot);

    // Draw the cone limit.
    ragdoll_drawcone(anchor_pt, anchor_pos, up_limit_angle, out_limit_angle,
                     matrix3(xform), guide_scale, cone_guide_color);

    // Draw the current orientation of the twist/up axis.
    vector twist_axis_child =
        normalize(qrotate(child_orient, twist_axis) * matrix3(current_xform));
    ragdoll_drawline(anchor_pt, anchor_pos + twist_axis_child * guide_scale,
                     twist_axis_color);

    vector up_axis_child =
        normalize(qrotate(child_orient, up_axis) * matrix3(current_xform));
    ragdoll_drawline(anchor_pt, anchor_pos + up_axis_child * guide_scale,
                     up_axis_color);

    // Draw twist guides.
    vector4 q_parent = quaternion(matrix3(xform));
    vector4 q_child = qmultiply(
        quaternion(matrix3(current_xform)), child_orient);

    // Compute current twist amount.
    vector4 q_diff = qmultiply(qinvert(q_parent), q_child);
    vector4 q_swing = dihedral(twist_axis, qrotate(q_diff, twist_axis));
    vector4 q_twist = qmultiply(qinvert(q_swing), q_diff);

    vector rotation_axis = qconvert(q_twist);
    float current_twist_angle = length(rotation_axis);
    if (dot(rotation_axis, twist_axis) < 0)
        current_twist_angle *= -1;

    // Rotate back by the current twist amount, and then draw the twist guides
    // in the child's space.
    float angles[];
    push(angles, -twist_limit_angle - current_twist_angle);
    push(angles, twist_limit_angle - current_twist_angle);
    foreach (float angle; angles)
    {
        vector limit_dir = qrotate(
            qmultiply(child_orient, quaternion(angle, twist_axis)), up_axis);
        limit_dir = normalize(limit_dir * matrix3(current_xform));

        ragdoll_drawline(anchor_pt,
                         anchor_pos + limit_dir * guide_scale,
                         twist_limit_guide_color);
    }
}

/// Compute the motor's maximum impulse, taking a transform group into account.
float
ragdoll_computemaximpulse(const string agent_geo; const int agent_ptnum;
                          const int agent_primnum; const int xform_idx;
                          const string impulse_attrib; const string state_table;
                          const int state_id; const string group_attrib)
{
    float max_impulse = point(agent_geo, impulse_attrib, agent_ptnum);

    string group_name = point(state_table, group_attrib, state_id);
    int group_idx =
        agentfindtransformgroup(agent_geo, agent_primnum, group_name);

    if (group_idx >= 0)
    {
        if (!agenttransformgroupmember(agent_geo, agent_primnum, group_idx,
                                       xform_idx))
        {
            max_impulse = 0;
        }
        else
        {
            max_impulse *= agenttransformgroupweight(agent_geo, agent_primnum,
                                                     group_idx, xform_idx);
        }
    }
    else if (len(group_name) > 0)
    {
        warning("Transform group '%s' not found.", group_name);
    }

    return max_impulse;
}

/// Compute the target orientation for a cone twist constraint's motor.
vector4
ragdoll_computemotortarget(const string agent_geo; const int agent_primnum;
                           const int xform_idx; const int parent_xform_idx;
                           const vector constrained_twist_axis;
                           const vector constrained_up_axis;
                           const vector goal_twist_axis;
                           const vector goal_up_axis)
{
    matrix3
    buildconexform(const vector twist_axis; const vector up_axis)
    {
        matrix3 rot = transpose(lookat({0, 0, 0}, {1, 0, 0}, {0, 0, 1}));
        rot *= lookat({0, 0, 0}, twist_axis, up_axis);
        return rot;
    }

    matrix agent_xform =
        primintrinsic(agent_geo, "packedfulltransform", agent_primnum);

    matrix joint_xform0 =
        agentworldtransform(agent_geo, agent_primnum, xform_idx) * agent_xform;
    matrix joint_xform1 =
        agentworldtransform(agent_geo, agent_primnum, parent_xform_idx) *
        agent_xform;

    matrix3 xform0 =
        buildconexform(constrained_twist_axis, constrained_up_axis) *
        matrix3(joint_xform0);
    matrix3 xform1 =
        buildconexform(goal_twist_axis, goal_up_axis) * matrix3(joint_xform1);

    return quaternion(xform0 * invert(xform1));
}

#endif
