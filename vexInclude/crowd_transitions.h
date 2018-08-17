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

#ifndef __crowd_transitions_h__
#define __crowd_transitions_h__

#include <locomotion.h>

/// Possible states that a crowd transition can be in.
#define TRANSITION_DELAY 0
#define TRANSITION_NEW_STATE 1
#define TRANSITION_NEW_CLIP 2
#define TRANSITION_ADVANCE 3

int
crowdtransition_findclippoint(const string clip_graph; const string agentname;
                              const string clipname)
{
    // TODO - it would be nice to have an array version of findattribval(),
    // similar to neighbours().

    int n = findattribvalcount(clip_graph, "point", "clipname", clipname);
    for (int i = 0; i < n; ++i)
    {
        int pt = findattribval(clip_graph, "point", "clipname", clipname, i);
        if (agentname == point(clip_graph, "agentname", pt))
            return pt;
    }

    return -1;
}

int
crowdtransition_findpath(const string clip_graph; const string path_geo;
                         const string agentname; const int clip_id_a;
                         const string clipname_b)
{
    int clip_id_b =
        crowdtransition_findclippoint(clip_graph, agentname, clipname_b);
    if (clip_id_b < 0)
        return -1;

    int n = findattribvalcount(path_geo, "prim", "end_clipname_id", clip_id_b);
    for (int i = 0; i < n; ++i)
    {
        int primnum =
            findattribval(path_geo, "prim", "end_clipname_id", clip_id_b, i);

        if (clip_id_a == prim(path_geo, "start_clipname_id", primnum))
            return primnum;
    }

    return -1;
}

/// Select a path in the clip graph to get to the next state. If the
/// destination state allows randomizing from a set of clips, only the
/// reachable clips are considered.
int
crowdtransition_selectpath(const string clip_graph; const string path_geo;
                           const string clip_properties;
                           const string state_table; const int primnum;
                           const string agentname; const int agent_id;
                           const int nextstate_id; const int clip_id_a)
{
    if (locomotion_randomizeclips(state_table, nextstate_id) == 1)
    {
        string clip_names[];
        float clip_weights[];
        locomotion_getrandomclipchoices(clip_properties, state_table,
                                        nextstate_id, agentname, primnum,
                                        clip_names, clip_weights);

        // Find the random clip choices that are actually reachable.
        int eligible_paths[];
        float path_weights[];
        foreach (int i; string clipname; clip_names)
        {
            int path_prim = crowdtransition_findpath(
                clip_graph, path_geo, agentname, clip_id_a, clipname);

            if (path_prim >= 0)
            {
                append(eligible_paths, path_prim);
                append(path_weights, clip_weights[i]);
            }
        }

        if (!len(eligible_paths))
            return -1;

        float seed = locomotion_clipnameseed(state_table, nextstate_id);
        float u = rand(set(agent_id, seed));

        return eligible_paths[sample_discrete(path_weights, u)];
    }
    else
    {
        string clipname = locomotion_clipname(state_table, nextstate_id);
        return crowdtransition_findpath(clip_graph, path_geo, agentname,
                                        clip_id_a, clipname);
    }
}

/// Inspect the clip transition graph to figure out how to transition from clip
/// A to clip B.
/// Returns TRANSITION_CHANGE_STATE or TRANSITION_CHANGE_CLIP if the transition
/// is allowed to start.
int
crowdtransition_find(const string clip_graph; const string path_geo;
                     const string state_table; const string clip_properties;
                     const int primnum; const string agentname;
                     const int agent_id; const int nextstate_id;
                     const string clip_a; const float clip_time_a;
                     export string clip_b; export float clip_time_b;
                     export float blend_duration; int path_prim;
                     int path_edge)
{
    float tolerance = 0.00001;

    if (npoints(clip_graph) == 0)
    {
        warning("Clip transition graph does not contain any points.");
        return TRANSITION_DELAY;
    }

    // If we're at the start of the transition, choose a path to a clip in the
    // new state.
    if (path_prim < 0)
    {
        int pt_a = crowdtransition_findclippoint(clip_graph, agentname, clip_a);
        if (pt_a < 0)
        {
            warning("Could not find clip '%s' in clip transition graph.",
                    clip_a);
            return TRANSITION_DELAY;
        }

        path_prim = crowdtransition_selectpath(
            clip_graph, path_geo, clip_properties, state_table, primnum,
            agentname, agent_id, nextstate_id, pt_a);

        if (path_prim < 0)
        {
            string nextstate = point(state_table, "state", nextstate_id);
            warning(
                "Could not find any transitions from clip '%s' to a clip of "
                "state '%s'",
                clip_a, nextstate);
            return TRANSITION_DELAY;
        }

        path_edge = 0;
    }

    // Find the pair of clips for our current edge in the path.
    int pt_a = point(path_geo, "clipname_id",
                     primpoint(path_geo, path_prim, path_edge));
    int pt_b = point(path_geo, "clipname_id",
                     primpoint(path_geo, path_prim, path_edge + 1));

    clip_b = point(clip_graph, "clipname", pt_b);

    // Find the edge in the clip transition graph.
    int hedge = pointhedge(clip_graph, pt_a, pt_b);
    if (hedge < 0)
        return TRANSITION_DELAY;

    int transition_prim = hedge_prim(clip_graph, hedge);
    if (transition_prim < 0)
        return TRANSITION_DELAY;

    // Check if the clip time is within a region where the transition can start.
    vector2 transition_regions[] =
        prim(clip_graph, "transition_regions", transition_prim);

    foreach(int i; vector2 transition_region; transition_regions)
    {
        if ((clip_time_a - transition_region[0]) >= -tolerance &&
            (clip_time_a - transition_region[1]) <= tolerance)
        {
            float blend_durations[] =
                prim(clip_graph, "blend_durations", transition_prim);
            vector2 sync_points[] =
                prim(clip_graph, "sync_points", transition_prim);

            // Use the sync points to select the initial clip time for clip B.
            clip_time_b = sync_points[i][1] + (clip_time_a - sync_points[i][0]);
            blend_duration = blend_durations[i];

            // Ensure that the new clip time is valid.
            int clip_b_id =
                crowdclip_findclipid(clip_properties, agentname, clip_b);
            float start_time = crowdclip_starttime(clip_properties, clip_b_id);
            vector2 loop_range = crowdclip_looprange(clip_properties, clip_b_id,
                                                     primnum, clip_b);
            float sample_rate = crowdclip_samplerate(clip_properties, clip_b_id,
                                                     primnum, clip_b);
            locomotion_wrapcliptime(state_table, nextstate_id, loop_range,
                                    start_time, sample_rate,
                                    /* enable_looping */ 1, clip_time_b);

            // Advance along the path. If we reach the end, also trigger the
            // state change.
            path_edge += 1;

            if (path_edge == primvertexcount(path_geo, path_prim) - 1)
            {
                path_prim = -1;
                path_edge = -1;
                return TRANSITION_NEW_STATE;
            }
            else
                return TRANSITION_NEW_CLIP;
        }
    }

    return TRANSITION_DELAY;
}

#endif
