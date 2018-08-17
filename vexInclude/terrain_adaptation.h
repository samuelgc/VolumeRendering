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
 * NAME:    terrain_adaptation.h
 */

#ifndef __terrain_adaptation__
#define __terrain_adaptation__

#include <agent_util.h>

#define FOOT_LOCK_TOLERANCE 0.01

/// Parameters used for terrain adaptation (version 1.0 and 2.0).
struct TerrainAdaptParms
{
    string terrain_path; /// Path to the terrain geometry.
    matrix agent_xform; /// Agent's overall transform ('packedfulltransform').
    matrix agent_xform_inv; /// Inverse of agent_xform.
    vector terrain_normal; /// The terrain normal (in world space).
    vector ref_dir; /// Direction the agent faces (in agent space).
    vector ref_up; /// Up vector (in agent space).

    float tilt_angle; /// The max angle to tilt per frame.
    float tilt_min; /// The minimum overall tilt angle.
    float tilt_max; /// The maximum overall tilt angle.
    float ankle_offset; /// Vertical offset for the ankle joint.
    float toe_offset; /// Vertical offset for the toe joint.
};

/// Parameters used for terrain adaptation (version 3.0).
struct TerrainAdaptParmsV3
{
    string terrain_path; /// Path to the terrain geometry.
    matrix agent_xform; /// Agent's overall transform ('packedfulltransform').
    matrix agent_xform_inv; /// Inverse of agent_xform.
    vector terrain_normal; /// The terrain normal (in world space).
    vector ref_dir; /// Direction the agent faces (in agent space).
    vector ref_up; /// Up vector (in agent space).

    float tilt_angle; /// The max angle to tilt per frame.
    float tilt_min; /// The minimum overall tilt angle.
    float tilt_max; /// The maximum overall tilt angle.

    int adjust_hips;
    float hip_offset;
    /// Limits how quickly the knee angle can change, once the upperleg ->
    /// ankle distance reaches this percentage of the maximum leg length.
    float knee_damping_threshold;
    /// Limit change in hip shift per frame.
    float hip_shift_per_frame;
};

struct TerrainAdaptLegData
{
    int upperleg_idx;
    int knee_idx;
    int ankle_idx;
    int toe_idx;
    vector ankle_offset;
    vector toe_offset;
    float ankle_down;
    float toe_down;

    int ankle_locked;
    int toe_locked;
    vector ankle_lock_pos;
    vector toe_lock_pos;

    vector knee_axis_local;

    void initFromArrays(const int primnum; const int i;
                        const string lowerlimb_joints[];
                        const vector knee_axes[];
                        const vector foot_offsets[]; const float foot_down[];
                        const int foot_locked[]; const vector foot_lock_pos[])
    {
        int j = i * 4;
        int k = i * 2;

        this.upperleg_idx = agentrigfind(0, primnum, lowerlimb_joints[j]);
        this.knee_idx = agentrigfind(0, primnum, lowerlimb_joints[j + 1]);
        this.ankle_idx = agentrigfind(0, primnum, lowerlimb_joints[j + 2]);
        this.toe_idx = agentrigfind(0, primnum, lowerlimb_joints[j + 3]);
        this.ankle_offset = foot_offsets[k];
        this.toe_offset = foot_offsets[k + 1];
        this.ankle_down = foot_down[k];
        this.toe_down = foot_down[k + 1];

        this.ankle_locked = foot_locked[k];
        this.toe_locked = foot_locked[k + 1];
        this.ankle_lock_pos = foot_lock_pos[k];
        this.toe_lock_pos = foot_lock_pos[k + 1];

        this.knee_axis_local = knee_axes[i];
    }

    void writeToArrays(const int i; vector knee_axes[]; int foot_locked[];
                       vector foot_lock_pos[])
    {
        knee_axes[i] = knee_axis_local;

        int k = i * 2;
        foot_locked[k] = ankle_locked;
        foot_locked[k + 1] = toe_locked;
        foot_lock_pos[k] = ankle_lock_pos;
        foot_lock_pos[k + 1] = toe_lock_pos;
    }

    float solveHipShift(const float target_length_2; const vector upperleg_pos;
                        const vector ankle_pos; const vector target_ankle_pos)
    {
        float adjusted_length_2 = length2(target_ankle_pos - upperleg_pos);

        // length2(target_ankle_pos - (upperleg_pos + delta)) = target_length ^ 2
        float a = 1;
        float b = -2 * (target_ankle_pos.y - upperleg_pos.y);
        float c = adjusted_length_2 - target_length_2;
        float delta_1, delta_2;

        if (solvequadratic(a, b, c, delta_1, delta_2))
        {
            // Pick solution that maintains the same leg direction.
            if (sign(upperleg_pos.y + delta_1 - target_ankle_pos.y) ==
                sign(upperleg_pos.y - ankle_pos.y))
            {
                return delta_1;
            }
            else
                return delta_2;
        }
        else
            return 0;
    }

    void computeTargetHipShift(const int primnum; const vector target_ankle_pos;
                               const int hip_idx; const vector ref_up;
                               matrix world_matrices[]; int valid_matrices[];
                               float target_shift; float shift_limit)
    {
        vector upperleg_pos = getMatrixPosition(getworldtransform(
            primnum, this.upperleg_idx, world_matrices, valid_matrices));
        vector knee_pos = getMatrixPosition(getworldtransform(
            primnum, this.knee_idx, world_matrices, valid_matrices));
        vector ankle_pos = getMatrixPosition(getworldtransform(
            primnum, this.ankle_idx, world_matrices, valid_matrices));

        float max_length =
            length(knee_pos - upperleg_pos) + length(ankle_pos - knee_pos);
        float max_length_2 = max_length * max_length;
        float orig_length_2 = length2(ankle_pos - upperleg_pos);

        // Solve for the vertical hip adjustment that would maintain the leg
        // length from the source animation.
        target_shift = this->solveHipShift(orig_length_2, upperleg_pos,
                                           ankle_pos, target_ankle_pos);

        // Solve for how far the hips can be shifted without stretching the leg.
        shift_limit = this->solveHipShift(max_length_2, upperleg_pos, ankle_pos,
                                          target_ankle_pos);
    }
};

/// Find the index in 'parent_joints' of the closest parent to the given joint.
int
TADAPTfindParent(const int primnum; const int child_xform_idx;
                 const int parent_joints[])
{
    int xform_idx = agentrigparent(0, primnum, child_xform_idx);

    while (xform_idx >= 0)
    {
        int i = find(parent_joints, xform_idx);
        if (i >= 0)
            return i;

        xform_idx = agentrigparent(0, primnum, xform_idx);
    }

    return -1;
}

void
TADAPTadjustToe(const int primnum; const int ankle_idx;
                const vector current_toe_pos; const vector target_toe_pos;
                const vector ankle_pos; matrix world_matrices[];
                int valid_matrices[])
{
    vector4 adjust_ankle =
        dihedral(current_toe_pos - ankle_pos, target_toe_pos - ankle_pos);

    matrix ankle_xform =
        getworldtransform(primnum, ankle_idx, world_matrices, valid_matrices);
    rotateAboutSelf(ankle_xform, adjust_ankle);
    world_matrices[ankle_idx] = ankle_xform;

    AGENTcomputeChildTransforms(primnum, ankle_idx, array(ankle_idx),
                                world_matrices, valid_matrices);
}

/// Adjusts the joint position by pushing it up to the terrain's surface. If
/// 'bidirectional' is enabled, it can also be pushed down if it is above the
/// terrain.
int
TADAPTfindTargetJointPos(const TerrainAdaptParms tadapt_parms;
                         const float joint_offset; const vector ray_dir_world;
                         const vector up_dir_world; const int bidirectional;
                         vector joint_pos)
{
    vector current_joint_pos_world = joint_pos * tadapt_parms.agent_xform;

    vector target_joint_pos_world;
    int found_intersection =
        lineIntersect(tadapt_parms.terrain_path, current_joint_pos_world,
                      ray_dir_world, 1, target_joint_pos_world);

    if (found_intersection >= 0)
    {
        vector delta = current_joint_pos_world - target_joint_pos_world;
        if (bidirectional ||
            projectedDistance(delta, up_dir_world) < joint_offset)
        {
            joint_pos = (target_joint_pos_world + joint_offset * up_dir_world) *
                        tadapt_parms.agent_xform_inv;

            return 1;
        }
    }

    return 0;
}

/// Adjusts the joint position, including an offset that is a position in local
/// space rather than a vertical offset in world space.
/// Returns 1 if the joint was pushed up.
int
TADAPTfindTargetJointPos(const int primnum;
                         const TerrainAdaptParmsV3 tadapt_parms;
                         const vector joint_offset; const vector ray_dir_world;
                         const int force_plant; matrix world_matrices[];
                         int valid_matrices[]; const int joint_idx;
                         vector joint_pos)
{
    matrix xform =
        getworldtransform(primnum, joint_idx, world_matrices, valid_matrices);
    vector delta = joint_offset * xform - getMatrixPosition(xform);

    vector current_pos_world = (joint_pos + delta) * tadapt_parms.agent_xform;

    // Search in both directions to find the closest intersection point, but
    // only use that point if the joint must be planted or it is below the
    // terrain.
    vector target_pos_world;
    int found_intersection =
        lineIntersect(tadapt_parms.terrain_path, current_pos_world,
                      ray_dir_world, /* bidirectional */ 1, target_pos_world);

    if (found_intersection >= 0)
    {
        vector adjustment = target_pos_world - current_pos_world;
        int pushed_up = dot(adjustment, ray_dir_world) > 0;
        if (pushed_up || force_plant)
        {
            joint_pos = target_pos_world * tadapt_parms.agent_xform_inv - delta;
            return pushed_up;
        }
    }

    return 0;
}

// For agentterrainadaptation 1.0 / 2.0.
vector
TADAPTcomputeIKPlane(const int primnum; const int upperleg_idx;
                     const int knee_idx; const int ankle_idx;
                     const vector target_ankle_pos; const vector ref_dir;
                     matrix world_matrices[]; int valid_matrices[])
{
    vector upperleg_pos = getMatrixPosition(getworldtransform(
        primnum, upperleg_idx, world_matrices, valid_matrices));
    vector knee_pos = getMatrixPosition(
        getworldtransform(primnum, knee_idx, world_matrices, valid_matrices));
    vector ankle_pos = getMatrixPosition(
        getworldtransform(primnum, ankle_idx, world_matrices, valid_matrices));

    // Determine the second direction vector that (along with target -
    // origin) defines the plane of the 2-bone IK solution.
    vector dir_1_adapted = normalize(target_ankle_pos - upperleg_pos);
    vector dir_1 = normalize(ankle_pos - upperleg_pos);
    vector dir_2 = normalize(knee_pos - ankle_pos);

    if (abs(dot(dir_1, dir_2)) >= 0.999f)
        dir_2 = ref_dir;

    vector4 adjust_leg = dihedral(dir_1, dir_1_adapted);
    return qrotate(adjust_leg, dir_2);
}

// For agentterrainadaptation 3.0
vector
TADAPTgetCurrentIKPlane(const int primnum; const int upperleg_idx;
                        const int knee_idx; const int ankle_idx;
                        const vector ref_dir; vector knee_axis_local;
                        matrix world_matrices[]; int valid_matrices[])
{
    vector upperleg_pos = getMatrixPosition(getworldtransform(
        primnum, upperleg_idx, world_matrices, valid_matrices));
    matrix knee_xform =
        getworldtransform(primnum, knee_idx, world_matrices, valid_matrices);
    vector ankle_pos = getMatrixPosition(
        getworldtransform(primnum, ankle_idx, world_matrices, valid_matrices));

    vector dir_1 = normalize(ankle_pos - upperleg_pos);
    vector knee_axis;

    if (length2(knee_axis_local) == 0)
    {
        vector dir_2 = normalize(getMatrixPosition(knee_xform) - ankle_pos);
        vector dir_3 = cross(dir_1, dir_2);

        if (length2(dir_3) < 0.005)
            return ref_dir;
        else
        {
            // Compute and store the normal to the plane of the IK solution, in
            // the space of the knee joint.
            knee_axis = normalize(dir_3);
            knee_axis_local =
                normalize(knee_axis * matrix3(invert(knee_xform)));
        }
    }
    else
        knee_axis = normalize(knee_axis_local * matrix3(knee_xform));

    // Determine the second direction vector that (along with target -
    // origin) defines the plane of the 2-bone IK solution.
    return normalize(cross(knee_axis, dir_1));
}

int
TADAPTadaptAnkleToTerrain(const int primnum; const TerrainAdaptParms parms;
                          const vector up; const int upperleg_idx;
                          const int knee_idx; const int ankle_idx;
                          matrix world_matrices[]; int valid_matrices[])
{
    vector ankle_pos = getMatrixPosition(
        getworldtransform(primnum, ankle_idx, world_matrices, valid_matrices));

    vector target_ankle_pos = ankle_pos;
    if (TADAPTfindTargetJointPos(parms, parms.ankle_offset, up, up, 0,
                                 target_ankle_pos))
    {
        vector dir_2 = TADAPTcomputeIKPlane(
            primnum, upperleg_idx, knee_idx, ankle_idx, target_ankle_pos,
            parms.ref_dir, world_matrices, valid_matrices);

        AGENTsolve2BoneIK(primnum, dir_2, target_ankle_pos, 1.0, world_matrices,
                          valid_matrices, upperleg_idx, knee_idx, ankle_idx);
        return 1;
    }

    return 0;
}

int
TADAPTadaptToeToTerrain(const int primnum; const TerrainAdaptParms parms;
                        const vector up; const int knee_idx;
                        const int ankle_idx; const int toe_idx;
                        matrix world_matrices[]; int valid_matrices[])
{
    vector knee_pos = getMatrixPosition(
        getworldtransform(primnum, knee_idx, world_matrices, valid_matrices));
    vector ankle_pos = getMatrixPosition(
        getworldtransform(primnum, ankle_idx, world_matrices, valid_matrices));
    vector toe_pos = getMatrixPosition(
        getworldtransform(primnum, toe_idx, world_matrices, valid_matrices));

    vector knee_pos_world = knee_pos * parms.agent_xform;
    vector ankle_pos_world = ankle_pos * parms.agent_xform;

    vector target_toe_pos = toe_pos;
    if (TADAPTfindTargetJointPos(parms, parms.toe_offset,
                                 normalize(knee_pos_world - ankle_pos_world),
                                 up, 0, target_toe_pos))
    {
        TADAPTadjustToe(primnum, ankle_idx, toe_pos, target_toe_pos, ankle_pos,
                        world_matrices, valid_matrices);
        return 1;
    }

    return 0;
}

/// Adapt a leg of the agent to the terrain.
void
TADAPTadjustLeg(const int primnum; const TerrainAdaptParms parms;
                const vector up; const int upperleg_idx; const int knee_idx;
                const int ankle_idx; const int toe_idx)
{
    matrix world_matrices[];
    int valid_matrices[];

    // Adjust the leg if necessary.
    TADAPTadaptAnkleToTerrain(primnum, parms, up, upperleg_idx, knee_idx,
                              ankle_idx, world_matrices, valid_matrices);

    // Adjust the toe if necessary.
    if (toe_idx >= 0)
    {
        TADAPTadaptToeToTerrain(primnum, parms, up, knee_idx, ankle_idx,
                                toe_idx, world_matrices, valid_matrices);
    }
}

void
TADAPTcomputeTargetAnklePosition(const int primnum;
                                 const TerrainAdaptParmsV3 tadapt_parms;
                                 const vector up; const TerrainAdaptLegData leg;
                                 const vector inital_ankle_pos;
                                 matrix world_matrices[]; int valid_matrices[];
                                 vector target_ankle_pos)
{
    // Adjust the initial ankle_position to the terrain.
    int force_plant_ankle =
        (!leg.ankle_locked && leg.ankle_down > FOOT_LOCK_TOLERANCE);

    vector adapted_ankle_pos = inital_ankle_pos;
    int ankle_pushed_up = TADAPTfindTargetJointPos(
        primnum, tadapt_parms, leg.ankle_offset, up, force_plant_ankle,
        world_matrices, valid_matrices, leg.ankle_idx, adapted_ankle_pos);

    // Adjust the target ankle position if the ankle or toe is locked.
    if (leg.ankle_locked)
    {
        // When locked, the channel value specifies how to blend between the
        // locked position and the current position.
        target_ankle_pos = leg.ankle_lock_pos * tadapt_parms.agent_xform_inv;
        target_ankle_pos =
            lerp(adapted_ankle_pos, target_ankle_pos, leg.ankle_down);
    }
    else if (!ankle_pushed_up && leg.ankle_down > FOOT_LOCK_TOLERANCE)
    {
        // If we're blending into the foot lock and the foot is above the
        // ground, we want to start blending down toward the terrain. We don't
        // actually choose the locked position until the channel value reaches
        // 1, though.
        target_ankle_pos =
            lerp(inital_ankle_pos, adapted_ankle_pos, leg.ankle_down);
    }
    else
        target_ankle_pos = adapted_ankle_pos;
}

void
TADAPTcomputeTargetToePosition(
    const int primnum; const TerrainAdaptParmsV3 tadapt_parms; const vector up;
    const TerrainAdaptLegData leg; const vector initial_toe_pos;
    const vector toe_handle_delta; matrix world_matrices[];
    int valid_matrices[]; vector target_toe_handle_pos; vector target_ankle_pos)
{
    // Adjust the target toe position for the terrain.
    vector target_toe_pos = initial_toe_pos;
    int toe_pushed_up = 0;
    {
        int force_plant_toe =
            (!leg.toe_locked && leg.toe_down > FOOT_LOCK_TOLERANCE);

        toe_pushed_up = TADAPTfindTargetJointPos(
            primnum, tadapt_parms, leg.toe_offset, up, force_plant_toe,
            world_matrices, valid_matrices, leg.toe_idx, target_toe_pos);
    }

    // Lock the toe if necessary.
    if (leg.toe_locked)
    {
        // Adjust the toe bone so that toe + toe_offset is at the target
        // position.
        target_toe_handle_pos =
            lerp(target_toe_pos + toe_handle_delta,
                 leg.toe_lock_pos * tadapt_parms.agent_xform_inv, leg.toe_down);
    }
    else if (!toe_pushed_up && leg.toe_down > FOOT_LOCK_TOLERANCE)
    {
        // Start moving down to the terrain, if necessary.
        target_toe_handle_pos =
            lerp(initial_toe_pos + toe_handle_delta,
                 target_toe_pos + toe_handle_delta, leg.toe_down);
    }
    else
    {
        target_toe_handle_pos = target_toe_pos + toe_handle_delta;
    }
}

void
TADAPTcomputeTargetFootPosition(const int primnum;
                                const TerrainAdaptParmsV3 tadapt_parms;
                                const vector up; TerrainAdaptLegData leg;
                                vector target_ankle_pos;
                                vector target_toe_handle_pos;
                                matrix world_matrices[]; int valid_matrices[])
{
    // Check if the foot lock is finishing.
    if (leg.toe_locked && leg.toe_down < FOOT_LOCK_TOLERANCE)
        leg.toe_locked = 0;
    if (leg.ankle_locked && leg.ankle_down < FOOT_LOCK_TOLERANCE)
        leg.ankle_locked = 0;

    vector initial_ankle_pos = getMatrixPosition(getworldtransform(
        primnum, leg.ankle_idx, world_matrices, valid_matrices));
    matrix initial_toe_xform =
        getworldtransform(primnum, leg.toe_idx, world_matrices, valid_matrices);
    vector initial_toe_pos = getMatrixPosition(initial_toe_xform);

    vector toe_to_ankle = initial_ankle_pos - initial_toe_pos;
    vector toe_handle_delta =
        leg.toe_offset * initial_toe_xform - initial_toe_pos;
    vector toe_handle_to_ankle = toe_to_ankle - toe_handle_delta;

    // If the ankle is planted and the toe is free, determine the toe position
    // after the target ankle position. Otherwise, we push up / plant the toe
    // and then determine the ankle position.
    int compute_toe_first =
        leg.toe_locked || leg.ankle_down < FOOT_LOCK_TOLERANCE;
    if (leg.toe_idx >= 0 && compute_toe_first)
    {
        TADAPTcomputeTargetToePosition(primnum, tadapt_parms, up, leg,
                                       initial_toe_pos, toe_handle_delta,
                                       world_matrices, valid_matrices,
                                       target_toe_handle_pos, target_ankle_pos);

        initial_ankle_pos = target_toe_handle_pos + toe_handle_to_ankle;
    }

    TADAPTcomputeTargetAnklePosition(primnum, tadapt_parms, up, leg,
                                     initial_ankle_pos, world_matrices,
                                     valid_matrices, target_ankle_pos);

    if (leg.toe_idx >= 0)
    {
        if (!compute_toe_first)
        {
            initial_toe_pos = target_ankle_pos - toe_to_ankle;

            TADAPTcomputeTargetToePosition(
                primnum, tadapt_parms, up, leg, initial_toe_pos,
                toe_handle_delta, world_matrices, valid_matrices,
                target_toe_handle_pos, target_ankle_pos);
        }
        else
        {
            vector foot_dir =
                normalize(target_ankle_pos - target_toe_handle_pos);
            target_ankle_pos =
                target_toe_handle_pos + foot_dir * length(toe_handle_to_ankle);
        }
    }
}

void
TADAPTadjustLegAndToe(const int primnum; const TerrainAdaptParmsV3 tadapt_parms;
                      const vector target_ankle_pos;
                      const vector target_toe_handle_pos;
                      TerrainAdaptLegData leg; matrix world_matrices[];
                      int valid_matrices[];)
{
    // 2-bone IK to adjust for the target ankle position.
    {
        vector dir_2 = TADAPTgetCurrentIKPlane(
            primnum, leg.upperleg_idx, leg.knee_idx, leg.ankle_idx,
            tadapt_parms.ref_dir, leg.knee_axis_local, world_matrices,
            valid_matrices);

        AGENTsolve2BoneIK(primnum, dir_2, target_ankle_pos,
                          tadapt_parms.knee_damping_threshold, world_matrices,
                          valid_matrices, leg.upperleg_idx, leg.knee_idx,
                          leg.ankle_idx);
    }

    // Adjust the toe position.
    if (leg.toe_idx >= 0)
    {
        matrix ankle_xform = getworldtransform(primnum, leg.ankle_idx,
                                               world_matrices, valid_matrices);
        matrix toe_xform = getworldtransform(primnum, leg.toe_idx,
                                             world_matrices, valid_matrices);

        vector toe_pos = getMatrixPosition(toe_xform);
        vector toe_handle_pos = leg.toe_offset * toe_xform;

        vector4 rotate_ankle =
            dihedral(toe_handle_pos - target_ankle_pos,
                     target_toe_handle_pos - target_ankle_pos);
        rotateAboutSelf(ankle_xform, rotate_ankle);

        world_matrices[leg.ankle_idx] = ankle_xform;
        int modified_xforms[] = array(leg.ankle_idx);

        AGENTcomputeChildTransforms(primnum, leg.ankle_idx, modified_xforms,
                                    world_matrices, valid_matrices);
    }

    // Check if the foot plant is starting.
    if (!leg.toe_locked && leg.toe_down > (1 - FOOT_LOCK_TOLERANCE))
    {
        leg.toe_locked = 1;

        matrix toe_xform = getworldtransform(primnum, leg.toe_idx,
                                             world_matrices, valid_matrices);
        leg.toe_lock_pos =
            leg.toe_offset * toe_xform * tadapt_parms.agent_xform;
    }

    if (!leg.ankle_locked && leg.ankle_down > (1 - FOOT_LOCK_TOLERANCE))
    {
        leg.ankle_locked = 1;

        vector ankle_pos = getMatrixPosition(getworldtransform(
            primnum, leg.ankle_idx, world_matrices, valid_matrices));
        leg.ankle_lock_pos = ankle_pos * tadapt_parms.agent_xform;
    }
}

void
TADAPTadjustHips(const int primnum; const TerrainAdaptParmsV3 tadapt_parms;
                 const int hip_indices[]; const float hip_shifts[];
                 matrix world_matrices[]; int valid_matrices[])
{
    // First, figure out the order in which to adjust the hips.
    // TODO - consider precomputing this in the Agent Prep SOP.
    int num_hips = len(hip_indices);
    int root = -1;
    int childs[];

    resize(childs, num_hips);
    for (int i = 0; i < num_hips; ++i)
        childs[i] = -1;

    foreach (int i; int hip_joint_idx; hip_indices)
    {
        int parent = TADAPTfindParent(primnum, hip_joint_idx, hip_indices);

        if (parent >= 0)
            childs[parent] = i;
        else
            root = i;
    }

    int adjustment_order[] = array(root);
    for (int i = 0; i < len(adjustment_order); ++i)
    {
        int child = childs[adjustment_order[i]];
        if (child >= 0)
            append(adjustment_order, child);
    }
    adjustment_order = reverse(adjustment_order);

    // Adjust the hips and spine.
    vector hip_positions_before[];
    vector hip_positions_after[];
    int modified_joints[];

    foreach (int hip_idx; adjustment_order)
    {
        int hip_joint_idx = hip_indices[hip_idx];

        matrix hips_xform = getworldtransform(primnum, hip_joint_idx,
                                              world_matrices, valid_matrices);
        vector hip_pos = getMatrixPosition(hips_xform);
        hip_positions_before[hip_idx] = hip_pos;

        hip_pos += hip_shifts[hip_idx] * tadapt_parms.ref_up;

        hip_positions_after[hip_idx] = hip_pos;
        setMatrixPosition(hips_xform, hip_pos);

        append(modified_joints, hip_joint_idx);
        world_matrices[hip_joint_idx] = hips_xform;

        // Rotate spine.
        int child_idx = childs[hip_idx];
        if (child_idx >= 0)
        {
            // Find the spine joint to rotate towards the child hips.
            int hips_children[] = agentrigchildren(0, primnum, hip_joint_idx);
            int spine_idx = TADAPTfindParent(primnum, hip_indices[child_idx],
                                             hips_children);
            if (spine_idx >= 0)
            {
                int spine_joint = hips_children[spine_idx];
                matrix spine_xform =
                    agentlocaltransform(0, primnum, spine_joint) * hips_xform;

                vector v1 = hip_positions_before[child_idx] -
                            hip_positions_before[hip_idx];
                vector v2 = hip_positions_after[child_idx] -
                            hip_positions_after[hip_idx];

                vector4 q = dihedral(v1, v2);
                rotateAboutSelf(spine_xform, q);

                append(modified_joints, spine_joint);
                world_matrices[spine_joint] = spine_xform;
            }
        }
    }

    AGENTcomputeChildTransforms(primnum, hip_indices[root], modified_joints,
                                world_matrices, valid_matrices);
}

/// Adapt an agent to terrain and/or perform foot locking, given a list of
/// lower limbs and torso(s).
void
TADAPTadaptAndLockLowerLimbs(
    const int primnum; const TerrainAdaptParmsV3 tadapt_parms;
    const string hip_joints[]; const string lowerback_joints[];
    const string lowerlimb_joints[]; vector tilt_dirs[]; vector knee_axes[];
    float hip_shifts[]; const vector foot_offsets[]; const float foot_down[];
    int foot_locked[]; vector foot_lock_pos[])
{
    vector terrain_normal_local = normalize(
        tadapt_parms.terrain_normal * matrix3(tadapt_parms.agent_xform_inv));
    vector up_world =
        normalize(tadapt_parms.ref_up * matrix3(tadapt_parms.agent_xform));

    int hip_indices[];
    foreach (string joint; hip_joints)
        append(hip_indices, agentrigfind(0, primnum, joint));

    matrix world_matrices[];
    int valid_matrices[];

    vector target_ankle_positions[];
    vector target_toe_handle_positions[];
    float target_shifts[];
    int leg_hip_joints[]; // Hip joint associated with each leg.

    // Track maximum allowed shift per hip joint.
    float shift_limits[];
    resize(shift_limits, len(hip_joints));
    for (int i = 0; i < len(shift_limits); ++i)
        shift_limits[i] = 1e5;

    // Compute the target positions for the ankles and toes.
    int num_lowerlimbs = len(lowerlimb_joints) / 4;
    for (int i = 0; i < num_lowerlimbs; ++i)
    {
        TerrainAdaptLegData leg;
        leg->initFromArrays(primnum, i, lowerlimb_joints, knee_axes,
                            foot_offsets, foot_down, foot_locked,
                            foot_lock_pos);

        TADAPTcomputeTargetFootPosition(
            primnum, tadapt_parms, up_world, leg, target_ankle_positions[i],
            target_toe_handle_positions[i], world_matrices, valid_matrices);

        if (tadapt_parms.adjust_hips)
        {
            int hip_joint_idx = -1;
            int hip_num =
                TADAPTfindParent(primnum, leg.upperleg_idx, hip_indices);

            if (hip_num < 0)
            {
                warning("Could not find hip joint for %s",
                        lowerlimb_joints[i * 4]);
            }
            else
                hip_joint_idx = hip_indices[hip_num];

            // For each leg, compute how far it would like to shift the hips
            // along with the limits to avoid over-stretching the leg.
            float leg_shift_limit;
            leg->computeTargetHipShift(
                primnum, target_ankle_positions[i], hip_joint_idx,
                tadapt_parms.ref_up, world_matrices, valid_matrices,
                target_shifts[i], leg_shift_limit);

            shift_limits[hip_num] = min(shift_limits[hip_num], leg_shift_limit);
            leg_hip_joints[i] = hip_joint_idx;
        }

        // Note: can't cache the leg structs in an array until bug 76016 is
        // fixed.
        leg->writeToArrays(i, knee_axes, foot_locked, foot_lock_pos);
    }

    if (tadapt_parms.adjust_hips)
    {
        // Compute the desired shift amount for each hip.
        foreach (int i; int hip_joint_idx; hip_indices)
        {
            float prev_hip_shift = hip_shifts[i];
            float min_target_shift = 1e5;
            float avg_target_shift = 0.0;
            int num_legs = 0;

            foreach (int j; float target_shift; target_shifts)
            {
                if (leg_hip_joints[j] == hip_joint_idx)
                {
                    min_target_shift = min(min_target_shift, target_shift);
                    avg_target_shift += target_shift;
                    ++num_legs;
                }
            }

            avg_target_shift /= num_legs;

            // Calculate the hip offset, preferring the average target offset
            // unless that's too close to the maximum offset. This is from
            // equation 7.17 in Johansen's paper ("Automated Semi-Procedural
            // Animation for Character Locomotion").
            float min_to_avg = avg_target_shift - min_target_shift;
            float min_to_max = shift_limits[i] - min_target_shift;
            hip_shifts[i] = min_target_shift + (min_to_avg * min_to_max) /
                                                   (min_to_avg + min_to_max);

            // Add in the global hip offset.
            hip_shifts[i] += tadapt_parms.hip_offset;

            // Optionally limit the hip adjustment rate.
            if (tadapt_parms.hip_shift_per_frame >= 0)
            {
                float delta = hip_shifts[i] - prev_hip_shift;
                delta = min(abs(delta),
                            tadapt_parms.hip_shift_per_frame) * sign(delta);
                hip_shifts[i] = prev_hip_shift + delta;
            }
        }

        // Adjust the skeleton.
        TADAPTadjustHips(primnum, tadapt_parms, hip_indices, hip_shifts,
                         world_matrices, valid_matrices);
    }

    // Tilt the agent back depending on the terrain.
    if (tadapt_parms.tilt_angle > 0)
    {
        int n = len(lowerback_joints);
        if (len(tilt_dirs) < n)
        {
            for (int i = len(tilt_dirs); i < n; ++i)
                append(tilt_dirs, tadapt_parms.ref_up);
        }

        for (int i = 0; i < n; ++i)
        {
            tilt_dirs[i] = AGENTtiltBack(
                primnum, tadapt_parms.ref_up, tadapt_parms.ref_dir,
                terrain_normal_local, tilt_dirs[i], tadapt_parms.tilt_angle,
                tadapt_parms.tilt_min, tadapt_parms.tilt_max,
                agentrigfind(0, primnum, lowerback_joints[i]));
        }
    }

    // Solve 2-bone IK for legs and adjust toes.
    for (int i = 0; i < num_lowerlimbs; ++i)
    {
        TerrainAdaptLegData leg;
        leg->initFromArrays(primnum, i, lowerlimb_joints, knee_axes,
                            foot_offsets, foot_down, foot_locked,
                            foot_lock_pos);

        TADAPTadjustLegAndToe(primnum, tadapt_parms, target_ankle_positions[i],
                              target_toe_handle_positions[i], leg,
                              world_matrices, valid_matrices);

        leg->writeToArrays(i, knee_axes, foot_locked, foot_lock_pos);
    }
}

float []
TADAPTsampleFootDownChannels(const int primnum; const string lowerlimb_joints[];
                             const string foot_channels[])
{
    int sample_clip(const int primnum; const string clipname;
                      const float cliptime; const string transformgroup;
                      const string channel; const string joint; float sample)
    {
        // Don't add warnings for missing channels if e.g. there is no toe
        // joint.
        int xform_idx = agentrigfind(0, primnum, joint);
        if (xform_idx < 0)
            return 0;

        // Skip layered clips that don't affect the lower body.
        int in_group =
            agenttransformgroupmember(0, primnum, transformgroup, xform_idx);
        if (!in_group)
            return 0;

        int channel_index = agentclipchannel(0, primnum, clipname, channel);
        if (channel_index < 0)
        {
            warning(
                "Foot plant channel '%s' for joint '%s' not found in clip '%s'",
                channel, joint, clipname);
            return 0;
        }

        sample = agentclipsample(0, primnum, clipname, cliptime, channel_index);
        return 1;
    }

    float sample_clips(const int primnum; const string clipnames[];
                       const float cliptimes[]; const float clipweights[];
                       const string transformgroups[];
                       const string channel_name; const string joint_name)
    {
        float sample = 0.0;
        sample_clip(primnum, clipnames[0], cliptimes[0], transformgroups[0],
                    channel_name, joint_name, sample);

        float sample_b = 0.0;
        if (len(clipnames) > 1 &&
            sample_clip(primnum, clipnames[1], cliptimes[1], transformgroups[1],
                        channel_name, joint_name, sample_b))
        {
            sample = lerp(sample, sample_b, clipweights[1]);
        }

        return sample;
    }

    float foot_down[];
    string clipnames[] = agentclipnames(0, primnum);
    float cliptimes[] = agentcliptimes(0, primnum);
    float clipweights[] = agentclipweights(0, primnum);
    string transformgroups[] = agentcliptransformgroups(0, primnum);

    int n = len(lowerlimb_joints) / 4;
    for (int i = 0; i < n; ++i)
    {
        float sample = sample_clips(primnum, clipnames, cliptimes, clipweights,
                                    transformgroups, foot_channels[i * 2],
                                    lowerlimb_joints[i * 4 + 2]);
        append(foot_down, sample);

        sample = sample_clips(primnum, clipnames, cliptimes, clipweights,
                              transformgroups, foot_channels[i * 2 + 1],
                              lowerlimb_joints[i * 4 + 3]);
        append(foot_down, sample);
    }

    return foot_down;
}

/// Builds the terrain adaptation guide geometry for the given agent. The agent
/// is removed and replaced with several points indicating the joint positions
/// and whether the feet are locked (using the 'P', 'Cd', and 'pscale' point
/// attributes).
void
TADAPTdrawGuides(const int primnum; const string color_ramp_parm;
                 const float locked_scale; const string lower_limbs[];
                 const vector foot_offsets[]; const float foot_down[])
{
    void add_point(const int primnum; const string color_ramp_parm;
                   const float locked_scale; const string joint_name;
                   const matrix agent_xform; const vector offset;
                   const float planted)
    {
        int xform_idx = agentrigfind(0, primnum, joint_name);
        if (xform_idx >= 0)
        {
            matrix xform =
                agentworldtransform(0, primnum, xform_idx) * agent_xform;
            int pt = addpoint(0, offset * xform);
            vector Cd = chramp(color_ramp_parm, planted);
            setpointattrib(0, "Cd", pt, Cd);

            vector scales =
                cracktransform(XFORM_SRT, XFORM_XYZ, 2, {0, 0, 0}, agent_xform);
            float guide_scale = avg(scales);
            if (planted >= 1 - FOOT_LOCK_TOLERANCE)
                guide_scale *= locked_scale;

            setpointattrib(0, "pscale", pt, guide_scale);
        }
    }

    int n = len(lower_limbs) / 4;
    matrix agent_xform = primintrinsic(0, "packedfulltransform", primnum);

    for (int i = 0; i < n; ++i)
    {
        string upper_leg = lower_limbs[i * 4];
        string knee = lower_limbs[i * 4 + 1];
        string ankle = lower_limbs[i * 4 + 2];
        string toe = lower_limbs[i * 4 + 3];

        add_point(primnum, color_ramp_parm, locked_scale, upper_leg,
                  agent_xform, {0, 0, 0}, 0);
        add_point(primnum, color_ramp_parm, locked_scale, knee, agent_xform,
                  {0, 0, 0}, 0);
        add_point(primnum, color_ramp_parm, locked_scale, ankle, agent_xform,
                  foot_offsets[i * 2], foot_down[i * 2]);
        add_point(primnum, color_ramp_parm, locked_scale, toe, agent_xform,
                  foot_offsets[i * 2 + 1], foot_down[i * 2 + 1]);
    }

    removeprim(0, primnum, 1);
}

void
TADAPTadaptAgentLowerLimbs(const int primnum;
                           const TerrainAdaptParms tadapt_parms;
                           const string lowerback_joints[];
                           const string upperleg_joints[];
                           const string knee_joints[];
                           const string ankle_joints[];
                           const string toe_joints[]; vector tilt_dir)
{
    vector terrain_normal_local = normalize(
        tadapt_parms.terrain_normal * matrix3(tadapt_parms.agent_xform_inv));
    vector up_world =
        normalize(tadapt_parms.ref_up * matrix3(tadapt_parms.agent_xform));

    // Tilt the agent back depending on the terrain.
    int n = len(lowerback_joints);
    for (int i = 0; i < n; ++i)
    {
        tilt_dir = AGENTtiltBack(primnum, tadapt_parms.ref_up,
                                 tadapt_parms.ref_dir, terrain_normal_local,
                                 tilt_dir, tadapt_parms.tilt_angle,
                                 tadapt_parms.tilt_min, tadapt_parms.tilt_max,
                                 agentrigfind(0, primnum, lowerback_joints[i]));
    }

    // Adjust lower limbs.
    n = len(upperleg_joints);
    for (int i = 0; i < n; ++i)
    {
        TADAPTadjustLeg(primnum, tadapt_parms, up_world,
                        agentrigfind(0, primnum, upperleg_joints[i]),
                        agentrigfind(0, primnum, knee_joints[i]),
                        agentrigfind(0, primnum, ankle_joints[i]),
                        agentrigfind(0, primnum, toe_joints[i]));
    }
}

///
/// Description: Adapt agent to terrain. Function will adjust left/right leg.
///
/// Parameters: 
///     * primnum - the primitive number of the agent
///     * agent_input - the input corresponding to the agent
///     * terrain_input - the input corresponding to the terrain
///     * tilt_angle - the angle to tilt per frame
///     * tilt_min - the minimum tilt angle
///     * tilt_max - the maximum tilt angle
///     * ankle_offset - the vertical offset to apply when raying the terrain
///     * toe_offset - TODO
///
void TADAPTadaptAgent(
    int primnum; 
    int ptnum;
    int agent_input; 
    string terrain_input_path; 
    float tilt_angle; 
    float tilt_min;
    float tilt_max;
    float ankle_offset;
    float toe_offset)
{
    string rootname = point(0, "agentrig_root", ptnum);
    string lowerbackname = point(0, "agentrig_lowerback", ptnum);
    string l_upperlegname = point(0, "agentrig_leftupperleg", ptnum);
    string l_kneename = point(0, "agentrig_leftknee", ptnum);
    string l_anklename = point(0, "agentrig_leftankle", ptnum);
    string l_toename = point(0, "agentrig_lefttoe", ptnum);
    string r_upperlegname = point(0, "agentrig_rightupperleg", ptnum);
    string r_kneename = point(0, "agentrig_rightknee", ptnum);
    string r_anklename = point(0, "agentrig_rightankle", ptnum);
    string r_toename = point(0, "agentrig_righttoe", ptnum);

    vector tiltdir = point(0, "agentterrainadaptation_tiltdir", ptnum);

    // Read in point attributes
    vector P = point(0, 'P', ptnum);
    vector up = normalize(point(0, 'up', ptnum));
    vector tnorm = normalize(point(0, 'terrainnormal', ptnum));
    vector4 qorient = point(0, 'orient', ptnum);
    vector4 qorient_inv = qinvert(qorient);

    matrix agent_xform = primintrinsic(0, "packedfulltransform", primnum);
    matrix agent_xform_inv = invert(agent_xform);

    // Default agent orientation along the positive z-axis.
    // TODO - don't have this hardcoded to z.
    vector orient = {0,0,1};

    // Rotate the terrain normal to be relative to the agent.
    tnorm = qrotate(qorient_inv, tnorm);

    // Tilt the agent back depending on the terrain.
    tiltdir = AGENTtiltBack(
        primnum,
        up,
        orient, 
        tnorm, 
        tiltdir,
        tilt_angle,
        tilt_min,
        tilt_max,
        agentrigfind(0, primnum, lowerbackname));
    setpointattrib(geoself(), "agentterrainadaptation_tiltdir", ptnum, tiltdir);

    TerrainAdaptParms parms;
    parms.terrain_path = terrain_input_path;
    parms.agent_xform = agent_xform;
    parms.agent_xform_inv = agent_xform_inv;
    parms.terrain_normal = tnorm;
    parms.ref_dir = orient;
    parms.ref_up = {0, 1, 0};
    parms.ankle_offset = ankle_offset;
    parms.toe_offset = toe_offset;

    // Adjust the left leg.
    TADAPTadjustLeg(
        primnum,
        parms,
        up,
        agentrigfind(0, primnum, l_upperlegname), 
        agentrigfind(0, primnum, l_kneename), 
        agentrigfind(0, primnum, l_anklename),
        agentrigfind(0, primnum, l_toename));

    // Adjust the right leg.
    TADAPTadjustLeg(
        primnum,
        parms,
        up,
        agentrigfind(0, primnum, r_upperlegname),
        agentrigfind(0, primnum, r_kneename), 
        agentrigfind(0, primnum, r_anklename),
        agentrigfind(0, primnum, r_toename));
}

#endif
