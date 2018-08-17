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
 * NAME:    locomotion.h
 */

#ifndef __locomotion_h__
#define __locomotion_h__

#include <agent_util.h>
#include <crowd_clips.h>
#include <crowd_cliplayers.h>

int
locomotion_findstateid(const string state_table; const string state)
{
    return findattribval(state_table, "point", "state", state);
}

/// Returns whether agents in the state should be active, static, or ignored by
/// the Bullet solver.
string
locomotion_ragdollstate(const string state_table; const int state_id)
{
    return point(state_table, "ragdoll", state_id);
}

/// Returns the gait speed for an in-place animation clip. For a locomotive
/// clip, -1 is returned.
float
locomotion_gaitspeed(const string state_table; const int state_id;
                     const string clip_properties; const int clip_id)
{
    float gaitspeed = point(state_table, "gaitspeed", state_id);

    // For in-place clips, gait speed can be overridden in the clip properties.
    if (gaitspeed >= 0 && clip_id >= 0)
    {
        float custom_gaitspeed = crowdclip_gaitspeed(clip_properties, clip_id);
        if (custom_gaitspeed >= 0)
            gaitspeed = custom_gaitspeed;
    }

    return gaitspeed;
}

/// Returns the allowed variance for the gait speed.
float
locomotion_gaitspeedvariance(const string state_table; const int state_id)
{
    return point(state_table, "speedvariance", state_id);
}

/// Returns whether the clip should be retimed based on the particle speed.
int
locomotion_shouldretimeclip(const string state_table; const int state_id)
{
    return point(state_table, "retime", state_id);
}

/// Returns whether the particle speed should be limited to fall within the
/// gait speed range.
int
locomotion_limitparticlespeed(const string state_table; const int state_id)
{
    return point(state_table, "limitparticlespeed", state_id);
}

/// Returns whether the last frame matches the first frame of the clip, and
/// should be skipped when incrementing the clip time. This is useful for
/// correctly looping locomotive clips (see comment in
/// locomotion_wrapcliptime).
int
locomotion_lastframematchesfirst(const string state_table; const int state_id)
{
    return point(state_table, "lastframematchesfirst", state_id);
}

/// Returns the clip name for a state if no randomization is enabled.
string
locomotion_clipname(const string state_table; const int state_id)
{
    return point(state_table, "clipname", state_id);
}

/// Returns whether the clips should be randomized for a state.
int
locomotion_randomizeclips(const string state_table; const int state_id)
{
    return point(state_table, "randomclipname", state_id);
}

/// Returns the seed used for randomizing clip names.
float
locomotion_clipnameseed(const string state_table; const int state_id)
{
    return point(state_table, "randomclipnameseed", state_id);
}

/// Returns the list of clip patterns used for random clip selection.
string []
locomotion_randomclippatterns(const string state_table; const int state_id)
{
    return point(state_table, "randomclippatterns", state_id);
}

/// Returns the list of clip pattern weights used for random clip selection.
float []
locomotion_randomclipweights(const string state_table; const int state_id)
{
    return point(state_table, "randomclipweights", state_id);
}

/// Blends between the speed limits for two states and handles states with no
/// speed limits.
float
locomotion_blendspeedlimit(const float limit_a; const float weight_a;
                           const float limit_b; const float weight_b)
{
    if (limit_a < 0)
        return limit_b;
    else if (limit_b < 0)
        return limit_a;
    else
        return limit_a * weight_a + limit_b * weight_b;
}

/// Returns a list of the clip names and weights that should be used when
/// randomly selecting a clip from the given list of clip name patterns and
/// weights.
void
locomotion_getrandomclipchoices(const string clip_properties;
                                const string agentname; const int primnum;
                                const string patterns[];
                                const float pattern_weights[];
                                string clip_names[]; float clip_weights[])
{
    foreach (int i; string pattern; patterns)
    {
        int num_matches = 0;
        foreach (string clipname;
                 crowdclip_catalog(clip_properties, agentname, primnum))
        {
            if (match(pattern, clipname))
            {
                append(clip_names, clipname);
                ++num_matches;
            }
        }

        // Split weight among the matched clips.
        if (num_matches > 0)
        {
            float w = pattern_weights[i] / num_matches;
            for (int j = 0; j < num_matches; ++j)
                append(clip_weights, w);
        }
    }
}

/// Returns a list of the clip names and weights that should be used when
/// randomly selecting a clip for the given state.
void
locomotion_getrandomclipchoices(const string clip_properties;
                                const string state_table; const int state_id;
                                const string agentname; const int primnum;
                                string clip_names[]; float clip_weights[])
{
    string patterns[] = locomotion_randomclippatterns(state_table, state_id);
    float pattern_weights[] =
        locomotion_randomclipweights(state_table, state_id);

    locomotion_getrandomclipchoices(clip_properties, agentname, primnum,
                                    patterns, pattern_weights, clip_names,
                                    clip_weights);
}

/// Select a clip for the agent based on the given state, not taking the clip
/// transition graph into account. See crowd_transitions.h for a more complex
/// version.
string
locomotion_chooseclip(const string clip_properties; const string state_table;
                      const int state_id; const string agentname; const int id;
                      const int primnum)
{
    if (locomotion_randomizeclips(state_table, state_id) == 1)
    {
        string clip_names[];
        float clip_weights[];
        locomotion_getrandomclipchoices(clip_properties, state_table, state_id,
                                        agentname, primnum, clip_names,
                                        clip_weights);

        float seed = locomotion_clipnameseed(state_table, state_id);
        float u = rand(set(id, seed));
        return clip_names[sample_discrete(clip_weights, u)];
    }
    else
        return locomotion_clipname(state_table, state_id);
}

int
locomotion_wrapcliptime(const string state_table; const int state_id;
                        const vector2 loop_range; const float start_time;
                        const float sample_rate; const int enable_looping;
                        export float clip_time)
{
    float tolerance = 0.001;
    float last_sample = loop_range[1];

    if (locomotion_lastframematchesfirst(state_table, state_id))
    {
        // If a locomotive clip is looping, assume that the last frame matches
        // the first (bug #67058). The last frame is used when sampling the
        // locomotion channel to get the correct translation when moving from
        // the last pose back to the original pose, but is skipped when
        // advancing the clip time to avoid having a double frame.
        float sample_t = 1.0 / sample_rate;
        last_sample -= sample_t;
    }

    // If the time is earlier than both the start time and the beginning of the
    // loop range, it's invalid.
    if ((clip_time - min(start_time, loop_range[0])) < -tolerance)
    {
        clip_time = loop_range[0] + (clip_time % last_sample);
        return 1;
    }

    if ((clip_time - last_sample) >= -tolerance)
    {
        // If looping is disabled, clamp clip time to the last sample.
        if (!enable_looping)
        {
            clip_time = last_sample;
            return 0;
        }

        if (clip_time < last_sample)
            clip_time = loop_range[0];
        else
            clip_time = loop_range[0] + (clip_time % last_sample);

        return 1;
    }
    else
        return 0;
}

/// Choose an initial clip time to start playback from, taking the clip
/// properties into account.
float
locomotion_choosecliptime(const int primnum; const string state_table;
                          const int state_id; const string clip_properties;
                          const int clip_id; const string clip_name;
                          const int id; const float current_clip_time;
                          const int override_current_clip_time;
                          const float clip_time_offset;
                          const int randomize_clip_time;
                          const float random_time_offset; const float seed)
{
    float start_time = crowdclip_starttime(clip_properties, clip_id);
    vector2 loop_range =
        crowdclip_looprange(clip_properties, clip_id, primnum, clip_name);
    float sample_rate =
        crowdclip_samplerate(clip_properties, clip_id, primnum, clip_name);

    float cliptime;
    if (override_current_clip_time)
    {
        cliptime = start_time + clip_time_offset;

        if (randomize_clip_time)
            cliptime += rand(set(id, seed)) * random_time_offset;
    }
    else
        cliptime = max(start_time, current_clip_time);

    // Ensure the clip time is valid (e.g. not in the last sample of a
    // looping locomotive clip.
    locomotion_wrapcliptime(state_table, state_id, loop_range, start_time,
                            sample_rate, /* enable_looping */ 1, cliptime);
    return cliptime;
}

matrix
locomotion_clipsampleworld(const string clip_properties; const int clip_id;
                           const int primnum; const string clip_name;
                           const float t; const int xform_idx)
{
    string source_clip = clip_name;
    if (clip_id >= 0)
        source_clip = crowdclip_sourcename(clip_properties, clip_id);

    return agentclipsampleworld(0, primnum, source_clip, t, xform_idx);
}

vector4
locomotion_clipsampleorient(const string clip_properties; const int clip_id;
                            const int primnum; const string clip_name;
                            const float t)
{
    int locomotion_idx = agentrigfind(0, primnum, "__locomotion__");
    matrix end_xform = locomotion_clipsampleworld(
        clip_properties, clip_id, primnum, clip_name, t, locomotion_idx);
    return quaternion(matrix3(end_xform));
}

/// Advance the animation clip.
void
locomotion_updatecliptime(
    const int primnum; const string clip_properties; const string state_table;
    const string state; const string agentname; const string clip_name;
    const vector v; const float pscale; const float timestep;
    const float locomotion_clip_retime; export float clip_time;
    export float locomotion_clip_time; export float clip_retime;
    int clip_looped; vector4 locomotion_frame)
{
    int state_id = locomotion_findstateid(state_table, state);
    int clip_id = crowdclip_findclipid(clip_properties, agentname, clip_name);

    float gait_speed =
        locomotion_gaitspeed(state_table, state_id, clip_properties, clip_id);
    float sample_rate =
        crowdclip_samplerate(clip_properties, clip_id, primnum, clip_name);
    float speed = length(v);
    int enable_looping = crowdclip_enablelooping(clip_properties, clip_id);
    vector2 loop_range =
        crowdclip_looprange(clip_properties, clip_id, primnum, clip_name);
    float start_time = crowdclip_starttime(clip_properties, clip_id);
    int is_locomotive = gait_speed < 0;

    // Retime in-place animation.
    if (gait_speed > 0 && locomotion_shouldretimeclip(state_table, state_id))
    {
        float allowed_variance =
            locomotion_gaitspeedvariance(state_table, state_id);

        gait_speed *= pscale;
        clip_retime *= clamp(speed / gait_speed, 1 - allowed_variance,
                             1 + allowed_variance);
    }

    float clip_timestep = clip_retime * timestep;
    float locomotion_timestep = locomotion_clip_retime * timestep;

    clip_time += clip_timestep;
    if (is_locomotive)
        locomotion_clip_time += locomotion_timestep;

    int wrapped =
        locomotion_wrapcliptime(state_table, state_id, loop_range, start_time,
                                sample_rate, enable_looping, clip_time);
    // This odd expression is needed to work around bug #72853.
    clip_looped = wrapped ? 1 : clip_looped;

    // Update the locomotion clip time.
    if (is_locomotive)
    {
        wrapped = locomotion_wrapcliptime(state_table, state_id, loop_range,
                                          start_time, sample_rate,
                                          enable_looping, locomotion_clip_time);

        // If the clip looped, update the rotation that is applied when
        // sampling the locomotion channel.
        if (wrapped)
        {
            int locomotion_idx = agentrigfind(0, primnum, "__locomotion__");
            matrix start_xform = locomotion_clipsampleworld(
                clip_properties, clip_id, primnum, clip_name, loop_range[0],
                locomotion_idx);
            matrix end_xform = locomotion_clipsampleworld(
                clip_properties, clip_id, primnum, clip_name,
                loop_range[1] - 1.0 / sample_rate, locomotion_idx);

            vector4 loop_orient =
                quaternion(matrix3((end_xform * invert(start_xform))));
            locomotion_frame = qmultiply(locomotion_frame, loop_orient);
        }
    }
}

void
locomotion_updateclips(const int primnum; const float pscale; const vector v;
                       const string agentname; const string clip_properties;
                       const string state_table; const string state;
                       const string next_state; const float timestep;
                       const int in_state_transition;
                       const int cliplayer_nextstate;
                       const float locomotion_clip_retimes[];
                       const string clip_names[]; const int clip_layerids[];
                       float clip_times[]; float locomotion_clip_times[];
                       float clip_retimes[]; int clip_looped[];
                       vector4 locomotion_frames[])
{
    for (int i = 0; i < len(clip_names); ++i)
    {
        // Only use the crowd state's retiming settings for the base animation
        // clips.
        string state_name;
        if (i == 0)
        {
            state_name = state;
        }
        else if (i == 1 && (crowdcliplayer_isbaselayer(clip_layerids[i],
                                                       cliplayer_nextstate)))
        {
            // If there aren't clip layers, during a state transition clip 1
            // belongs to layer 0 but is associated with the next state.
            state_name = in_state_transition ? next_state : state;
        }

        locomotion_updatecliptime(
            primnum, clip_properties, state_table, state_name, agentname,
            clip_names[i], v, pscale, timestep, locomotion_clip_retimes[i],
            clip_times[i], locomotion_clip_times[i], clip_retimes[i],
            clip_looped[i], locomotion_frames[i]);
    }
}

/// Limits the rotation from a to b based on the max turn rate.
vector
locomotion_limitturn(const vector a; const vector b; const int project_in_plane;
                     const vector up; const float max_turn_rate;
                     const float stiffness; const float damping;
                     const float max_ang_accel; const float timestep;
                     float current_turn_rate)
{
    // Compute the signed angle between the vectors.
    vector axis = cross(a, b);
    float dp = dot(a, b);
    float delta_angle = atan2(length(axis), dp);

    if (dot(axis, up) < 0)
    {
        delta_angle *= -1;
        axis *= -1;
    }

    // If the vectors are antiparallel and the forces have been projected into
    // the agent's plane, make sure to rotate around the up axis instead of
    // some arbitrary axis.
    if (project_in_plane)
        axis = up;

    // Determine angular acceleration as a function of the angle from the
    // target and the current turn rate.
    // This is similar to equation 1 in "Behavioural Dynamics of Human
    // Locomotion" (Warren & Fajen, 2004).
    if (stiffness >= 0)
    {
        delta_angle = agent_integratespring(
            stiffness, damping, radians(max_ang_accel), timestep, -delta_angle,
            radians(current_turn_rate));
    }

    // Limit the max turn rate.
    if (max_turn_rate >= 0)
    {
        float max_delta_angle = radians(max_turn_rate) * timestep;
        delta_angle = min(max_delta_angle, abs(delta_angle)) * sign(delta_angle);
    }

    current_turn_rate = degrees(delta_angle) / timestep;

    // Rotate 'a' around 'axis' towards 'b'.
    vector4 turn = quaternion(delta_angle, normalize(axis));
    return qrotate(turn, normalize(a) * length(b));
}

vector
locomotion_limitturn(const vector a; const vector b; const int project_in_plane;
                     const vector up; const float max_turn_rate;
                     const float timestep)
{
    float turn_rate;
    return locomotion_limitturn(a, b, project_in_plane, up, max_turn_rate, -1,
                                -1, -1, timestep, turn_rate);
}

void
locomotion_limittilt(const vector4 orient; const vector ref_dir;
                     const vector ref_up; const vector target_up;
                     const float max_tilt_rate; const float stiffness;
                     const float damping; const float max_ang_accel;
                     const float timestep; float tilt_rate; vector heading;
                     export vector limited_up)
{
    vector current_up = qrotate(orient, ref_up);

    limited_up = locomotion_limitturn(
        current_up, target_up, 0, {0, 0, 0}, max_tilt_rate, stiffness, damping,
        max_ang_accel, timestep, tilt_rate);

    vector4 delta_r = dihedral(current_up, limited_up);
    heading = qrotate(delta_r, heading);
}

vector4
locomotion_lookat(const vector ref_dir; const vector dir; const vector ref_up;
                  const vector up)
{
    matrix3 rot = transpose(lookat({0, 0, 0}, ref_dir, ref_up));
    rot *= lookat({0, 0, 0}, dir, up);

    return quaternion(rot);
}

/// Determine a target velocity, heading, and external forces for an agent,
/// depending on whether its clip is in-place or locomotive animation.
void
locomotion_getgoalxform(
    const int primnum; const float pscale; const vector v; const vector up;
    const vector ref_dir; const vector ref_up; const vector current_force;
    const vector current_heading; const string clip_properties;
    const string state_table; const string state; const string agentname;
    const string clip_name; const float locomotion_clip_time;
    const float locomotion_clip_retime; const float timestep;
    const float orient_change_threshold; const float max_turn_rate;
    const float turn_stiffness; const float turn_damping;
    const float turn_accel_max; const int constrain_v; const int project_force;
    const float sim_influence; const float locomotion_strength;
    const float inplace_airresist; export vector4 locomotion_frame;
    export vector heading; export vector force; export vector targetv;
    export float airresist; export float target_speed_min;
    export float target_speed_max; float turn_rate)
{
    int state_id = locomotion_findstateid(state_table, state);
    int clip_id = crowdclip_findclipid(clip_properties, agentname, clip_name);
    float gait_speed =
        locomotion_gaitspeed(state_table, state_id, clip_properties, clip_id);
    float speed = length(v);
    float tolerance = 0.001;

    force = current_force;
    if (project_force)
        force -= project(force, up);

    if (gait_speed < 0)
    {
        // Determine a target velocity and orientation from the locomotion
        // channel.
        int locomotion_idx = agentrigfind(0, primnum, "__locomotion__");

        // Find the current translation/rotation.
        matrix cog_xform_1 = locomotion_clipsampleworld(
            clip_properties, clip_id, primnum, clip_name, locomotion_clip_time,
            locomotion_idx);
        vector4 rot_1 = quaternion(matrix3(cog_xform_1));
        vector trans_1 = getMatrixPosition(cog_xform_1);

        // Find the next time to sample the locomotion channel at, wrapping
        // around if necessary.
        float t_2 = locomotion_clip_time + timestep * locomotion_clip_retime;
        float sample_rate =
            crowdclip_samplerate(clip_properties, clip_id, primnum, clip_name);

        int enable_looping = crowdclip_enablelooping(clip_properties, clip_id);
        vector2 loop_range =
            crowdclip_looprange(clip_properties, clip_id, primnum, clip_name);
        float start_time = crowdclip_starttime(clip_properties, clip_id);

        int wrapped = locomotion_wrapcliptime(state_table, state_id, loop_range,
                                              start_time, sample_rate,
                                              enable_looping, t_2);

        matrix cog_xform_2 = locomotion_clipsampleworld(
            clip_properties, clip_id, primnum, clip_name, t_2, locomotion_idx);
        vector trans_2 = getMatrixPosition(cog_xform_2);
        vector4 rot_2 = quaternion(matrix3(cog_xform_2));

        // If the clip time wraps at some point during the timestep, some extra
        // work is necessary to figure out the correct translation and
        // rotation.
        if (wrapped)
        {
            matrix cog_xform_start = locomotion_clipsampleworld(
                clip_properties, clip_id, primnum, clip_name, loop_range[0],
                locomotion_idx);
            vector4 rot_start = quaternion(matrix3(cog_xform_start));

            trans_2 -= getMatrixPosition(cog_xform_start);
            rot_2 = qmultiply(qinvert(rot_start), rot_2);

            float last_sample = loop_range[1];
            if (locomotion_lastframematchesfirst(state_table, state_id))
                last_sample -= 1.0 / sample_rate;

            matrix cog_xform_end = locomotion_clipsampleworld(
                clip_properties, clip_id, primnum, clip_name, last_sample,
                locomotion_idx);
            vector4 rot_end = quaternion(matrix3(cog_xform_end));

            trans_2 =
                qrotate(rot_end, trans_2) + getMatrixPosition(cog_xform_end);
            rot_2 = qmultiply(rot_end, rot_2);
        }

        // Determine how much the agent should rotate by according to the
        // locomotion channel.
        vector4 delta_rot = qmultiply(qinvert(rot_1), rot_2);
        vector4 delta_rot_world = qmultiply(
            locomotion_frame, qmultiply(delta_rot, qinvert(locomotion_frame)));
        vector locomotion_heading =
            qrotate(delta_rot_world, normalize(current_heading));

        // If the agent is below the orientation change threshold, or the
        // external forces are zero, use the current heading as the target
        // direction.
        vector target_heading = normalize(force);
        if (length(v) < orient_change_threshold || length2(force) < 0.00001)
            target_heading = locomotion_heading;

        // Otherwise, allow external forces to influence the agent, limited
        // by the max turn rate and sim influence parameters.
        target_heading = locomotion_limitturn(
            locomotion_heading, target_heading, project_force, up,
            max_turn_rate * sim_influence, turn_stiffness, turn_damping,
            turn_accel_max, timestep, turn_rate);

        vector4 orient = locomotion_lookat(ref_dir, target_heading, ref_up, up);
        heading = qrotate(orient, ref_dir);

        // Take the orientation changes from external forces into account
        // when determining the target velocity.
        vector4 locomotion_orient = locomotion_lookat(
            ref_dir, locomotion_heading, ref_up, up);
        locomotion_frame = qmultiply(
            locomotion_frame, qmultiply(qinvert(locomotion_orient), orient));

        vector trans = pscale * (trans_2 - trans_1);
        targetv = qrotate(locomotion_frame, trans) / timestep;

        force *= sim_influence;
        airresist = locomotion_strength;

        // Lock speed to match locomotion.
        target_speed_min = length(trans) / timestep;
        target_speed_max = target_speed_min;
    }
    else
    {
        // For in-place animation, prevent orientation changes if we have a
        // small velocity.
        if (speed >= max(tolerance, orient_change_threshold))
            heading = -cross(cross(normalize(v), up), up);
        else
            heading = current_heading;

        // The velocity constraint already enforces the angular acceleration
        // and speed limits, so we don't need to limit the orientation change
        // in that situation.
        if (!constrain_v)
        {
            heading = locomotion_limitturn(
                current_heading, heading, 0, up, max_turn_rate, turn_stiffness,
                turn_damping, turn_accel_max, timestep, turn_rate);
        }

        // External forces are unchanged, and there is an optional drag force.
        targetv = {0, 0, 0};
        airresist = inplace_airresist;

        // Limit particle speed to allowed variance.
        if (locomotion_limitparticlespeed(state_table, state_id))
        {
            float allowed_variance =
                locomotion_gaitspeedvariance(state_table, state_id);

            target_speed_min = (1 - allowed_variance) * gait_speed * pscale;
            target_speed_max = (1 + allowed_variance) * gait_speed * pscale;
        }
        else
        {
            target_speed_min = -1;
            target_speed_max = -1;
        }
    }

    if (constrain_v)
    {
        // If v is initially zero, use the initial heading.
        vector v_dir = speed > tolerance ? v : heading;

        // Don't allow turning while under the speed threshold.
        float turn_limit = speed >= orient_change_threshold ? max_turn_rate : 0;

        vector predicted_v = v + force * timestep;

        float v_turn_rate = turn_rate;
        predicted_v = locomotion_limitturn(
            v_dir, predicted_v, project_force, up, turn_limit, turn_stiffness,
            turn_damping, turn_accel_max, timestep, v_turn_rate);
        // Don't write out the turn rate again for locomotive clips.
        if (gait_speed >= 0)
            turn_rate = v_turn_rate;

        force = (predicted_v - v) / timestep;
    }
}

/// Determine a target velocity, heading, and external forces for an agent, and
/// blend between locomotive and in-place animation when transitioning between
/// states.
void
locomotion_steerparticle(
    const int primnum; const float pscale; const vector v; const vector up;
    const vector ref_dir; const vector ref_up; const string clip_properties;
    const string state_table; const string state; const string next_state;
    const int in_state_transition; const int cliplayer_nextstate;
    const string agentname; const string clip_names[];
    const int clip_layerids[]; const float locomotion_clip_times[];
    const float locomotion_clip_retimes[]; const float clip_weights[];
    const float timestep; const float orient_change_threshold;
    const float max_turn_rate; const float max_tilt_rate;
    const float turn_stiffness; const float turn_damping;
    const float turn_accel_max; const float tilt_stiffness;
    const float tilt_damping; const float tilt_accel_max; const int constrain_v;
    const int project_force; const float sim_influence;
    const float locomotion_strength; const float inplace_airresist;
    export vector4 locomotion_frames[]; export vector4 orient;
    export vector heading; export vector force; export vector targetv;
    export float airresist; export float speed_min; export float speed_max;
    float turn_rate; float tilt_rate)
{
    // Apply max tilt rate.
    vector up_limited;
    locomotion_limittilt(orient, ref_dir, ref_up, up, max_tilt_rate,
                         tilt_stiffness, tilt_damping, tilt_accel_max, timestep,
                         tilt_rate, heading, up_limited);

    // Update facing direction, targetv, etc
    vector current_force, current_heading;
    locomotion_getgoalxform(
        primnum, pscale, v, up_limited, ref_dir, ref_up, force, heading,
        clip_properties, state_table, state, agentname, clip_names[0],
        locomotion_clip_times[0], locomotion_clip_retimes[0], timestep,
        orient_change_threshold, max_turn_rate, turn_stiffness, turn_damping,
        turn_accel_max, constrain_v, project_force, sim_influence,
        locomotion_strength, inplace_airresist, locomotion_frames[0],
        current_heading, current_force, targetv, airresist, speed_min,
        speed_max, turn_rate);

    // Only include the base animation clips for locomotion.
    if (len(clip_names) > 1 &&
        crowdcliplayer_isbaselayer(clip_layerids[1], cliplayer_nextstate))
    {
        vector next_heading, next_targetv, next_force;
        float next_airresist, next_speed_min, next_speed_max;

        string state_b = in_state_transition ? next_state : state;
        locomotion_getgoalxform(
            primnum, pscale, v, up_limited, ref_dir, ref_up, force, heading,
            clip_properties, state_table, state_b, agentname, clip_names[1],
            locomotion_clip_times[1], locomotion_clip_retimes[1], timestep,
            orient_change_threshold, max_turn_rate, turn_stiffness,
            turn_damping, turn_accel_max, constrain_v, project_force,
            sim_influence, locomotion_strength, inplace_airresist,
            locomotion_frames[1], next_heading, next_force, next_targetv,
            next_airresist, next_speed_min, next_speed_max, turn_rate);

        heading =
            current_heading * clip_weights[0] + next_heading * clip_weights[1];
        force = current_force * clip_weights[0] + next_force * clip_weights[1];
        targetv = targetv * clip_weights[0] + next_targetv * clip_weights[1];
        airresist =
            airresist * clip_weights[0] + next_airresist * clip_weights[1];

        speed_min = locomotion_blendspeedlimit(speed_min, clip_weights[0],
                                               next_speed_min, clip_weights[1]);
        speed_max = locomotion_blendspeedlimit(speed_max, clip_weights[0],
                                               next_speed_max, clip_weights[1]);
    }
    else
    {
        heading = current_heading;
        force = current_force;
    }

    // Update orientation.
    vector4 prev_orient = orient;
    orient = locomotion_lookat(ref_dir, heading, ref_up, up_limited);
    if (dot(prev_orient, orient) < 0)
        orient *= -1;
}

#endif
