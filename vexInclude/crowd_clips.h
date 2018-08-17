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
 * NAME:    crowd_clips.h
 *
 * DESCRIPTION: Functions for dealing with the trimmed and renamed clips
 * defined in the clip properties geometry. See the Agent Clip Properties SOP's
 * documentation.
 *
 */

#ifndef __crowd_clips_h__
#define __crowd_clips_h__

int
crowdclip_findclipid(const string geometry; const string agentname;
                     const string clipname)
{
    if (npoints(geometry) == 0)
        return -1;

    // TODO - it would be nice to have an array version of findattribval(),
    // similar to neighbours().

    int n = findattribvalcount(geometry, "point", "clipname", clipname);
    for (int i = 0; i < n; ++i)
    {
        int pt = findattribval(geometry, "point", "clipname", clipname, i);
        if (agentname == point(geometry, "agentname", pt))
            return pt;
    }

    return -1;
}

int
crowdclip_enablelooping(const string clip_properties; const int clip_id)
{
    int success = 0;
    int enable =
        pointattrib(clip_properties, "enable_looping", clip_id, success);
    return success ? enable : 1;
}

float
crowdclip_blendtimebefore(const string clip_properties; const int clip_id)
{
    return point(clip_properties, "blend_duration_before", clip_id);
}

float
crowdclip_blendtimeafter(const string clip_properties; const int clip_id)
{
    return point(clip_properties, "blend_duration_after", clip_id);
}

vector2
crowdclip_looprange(const string clip_properties; const int clip_id)
{
    return point(clip_properties, "loop_range", clip_id);
}

vector2
crowdclip_looprange(const string clip_properties; const int clip_id;
                    const int primnum; const string clipname)
{
    if (clip_id >= 0)
        return crowdclip_looprange(clip_properties, clip_id);
    else
        return set(0, agentcliplength(0, primnum, clipname));
}

float
crowdclip_starttime(const string clip_properties; const int clip_id)
{
    return point(clip_properties, "start_time", clip_id);
}

float
crowdclip_gaitspeed(const string clip_properties; const int clip_id)
{
    return point(clip_properties, "gait_speed", clip_id);
}

float
crowdclip_samplerate(const string clip_properties; const int clip_id;
                     const int opinput; const int primnum;
                     const string clipname)
{
    if (clip_id >= 0)
        return point(clip_properties, "sample_rate", clip_id);
    else
        return agentclipsamplerate(opinput, primnum, clipname);
}

float
crowdclip_samplerate(const string clip_properties; const int clip_id;
                     const int primnum; const string clipname)
{
    return crowdclip_samplerate(clip_properties, clip_id, 0, primnum, clipname);
}

/// Name of the underlying animation clip.
string
crowdclip_sourcename(const string clip_properties; const int clip_id)
{
    return point(clip_properties, "clipname_source", clip_id);
}

/// Name of the virtual clip, which may not necessarily be the name of the
/// underlying animation clip.
string
crowdclip_name(const string clip_properties; const int clip_id)
{
    return point(clip_properties, "clipname", clip_id);
}

/// Extension of the agentclipcatalog() VEX function that includes any virtual
/// clips from the clip properties.
string[]
crowdclip_catalog(const string clip_properties; const string agentname;
                  const int opinput; const int primnum)
{
    string clipnames[] = agentclipcatalog(opinput, primnum);

    for (int i = 0, n = npoints(clip_properties); i < n; ++i)
    {
        string clip = crowdclip_name(clip_properties, i);
        string source_clip = crowdclip_sourcename(clip_properties, i);
        if (clip != source_clip)
            append(clipnames, clip);
    }

    return clipnames;
}

string[]
crowdclip_catalog(const string clip_properties; const string agentname;
                  const int primnum)
{
    return crowdclip_catalog(clip_properties, agentname, 0, primnum);
}

/// Blend clip i with itself as it loops around.
void
crowdclip_blend(const string clip_properties; const int primnum;
                const string agentname; const int i; const int clip_looped[];
                string clip_names[]; float clip_times[]; float clip_weights[];
                string clip_transformgroups[]; int clip_layerids[])
{
    string clipname = clip_names[i];
    int clip_id = crowdclip_findclipid(clip_properties, agentname, clipname);
    if (clip_id < 0)
        return;

    // Replace the virtual clip with its source clip name.
    clip_names[i] = crowdclip_sourcename(clip_properties, clip_id);

    float time_before = crowdclip_blendtimebefore(clip_properties, clip_id);
    float time_after = crowdclip_blendtimeafter(clip_properties, clip_id);
    if (time_before <= 0 && time_after <= 0)
        return;

    float sample_rate =
        crowdclip_samplerate(clip_properties, clip_id, primnum, clipname);
    vector2 loop_range = crowdclip_looprange(clip_properties, clip_id);
    float start_sample = loop_range[0];
    float end_sample = loop_range[1] - 1 / sample_rate;

    float t_after = start_sample + time_after;
    float t_before = end_sample - time_before;

    int j = len(clip_names);
    float w = clip_weights[i];
    int need_blend = 0;

    if (clip_times[i] < t_after && clip_times[i] >= start_sample &&
        clip_looped[i])
    {
        // After looping around, blend out of the clip's last sample.
        clip_weights[i] = fit(clip_times[i], start_sample, t_after, 0.5 * w, w);
        clip_weights[j] = w - clip_weights[i];
        clip_times[j] = end_sample;
        need_blend = 1;
    }
    else if (clip_times[i] >= t_before)
    {
        // Blend towards the clip's first sample.
        clip_weights[i] = fit(clip_times[i], t_before, end_sample, w, 0.5 * w);
        clip_weights[j] = w - clip_weights[i];
        clip_times[j] = start_sample;
        need_blend = 1;
    }

    if (need_blend)
    {
        clip_names[j] = clip_names[i];
        clip_transformgroups[j] = clip_transformgroups[i];
        clip_layerids[j] = clip_layerids[i];
    }
}

/// Helper functions for unit conversion in SOPs that allow parameters to
/// specify either samples or frames.
struct UnitConvertParms
{
    int unit_type;
    float fps;

    float sample_rate;
    float start_time;

    void setclip(const string clip_properties; const string agent_name;
                 const int opinput; const int primnum; const string clip_name)
    {
        int clip_id =
            crowdclip_findclipid(clip_properties, agent_name, clip_name);

        sample_rate = crowdclip_samplerate(clip_properties, clip_id, opinput,
                                           primnum, clip_name);
        start_time = crowdclip_starttime(clip_properties, clip_id);
    }
};

float
crowdclip_timetoseconds(const float val; const UnitConvertParms parms)
{
    if (parms.unit_type == 0) // Frames
    {
        return parms.start_time + (val - 1) / parms.fps;
    }
    else // Samples
    {
        return parms.start_time + val / parms.sample_rate;
    }
}

float
crowdclip_durationtoseconds(const float val; const UnitConvertParms parms)
{
    return val / (parms.unit_type == 0 ? parms.fps : parms.sample_rate);
}

int
crowdclip_durationtosamples(const float val; const UnitConvertParms parms)
{
    return int(crowdclip_durationtoseconds(val, parms) * parms.sample_rate);
}

float
crowdclip_secondstoduration(const float seconds; const UnitConvertParms parms)
{
    return seconds * (parms.unit_type == 0 ? parms.fps : parms.sample_rate);
}

float
crowdclip_secondstotime(const float seconds; const UnitConvertParms parms)
{
    if (parms.unit_type == 0) // Frames
    {
        return 1 + parms.fps * (seconds - parms.start_time);
    }
    else // Samples
    {
        return parms.sample_rate * (seconds - parms.start_time);
    }
}

// Insert an animation clip with appropriate default values at the specified
// index in the clip list.
void
crowdclip_insert(const int i; string clip_names[]; float clip_times[];
                 float clip_weights[]; float clip_retimes[]; int clip_looped[];
                 string clip_transformgroups[]; int clip_layerids[];
                 int clip_layerblendstates[]; float clip_layerblendtimes[];
                 float locomotion_clip_times[]; vector4 locomotion_frames[];
                 float locomotion_clip_retimes[])
{
    insert(clip_names, i, "");
    insert(clip_times, i, 0.0);
    insert(clip_weights, i, 0.0);
    insert(clip_retimes, i, 0.0);
    insert(clip_looped, i, 0);

    insert(clip_transformgroups, i, "");
    insert(clip_layerids, i, 0);
    insert(clip_layerblendstates, i, 0);
    insert(clip_layerblendtimes, i, 0);

    insert(locomotion_clip_times, i, 0.0);
    insert(locomotion_frames, i, {0, 0, 0, 1});
    insert(locomotion_clip_retimes, i, 0.0);
}

// Remove an animation clip from the specified index in the clip list.
void
crowdclip_remove(const int i; string clip_names[]; float clip_times[];
                 float clip_weights[]; float clip_retimes[]; int clip_looped[];
                 string clip_transformgroups[]; int clip_layerids[];
                 int clip_layerblendstates[]; float clip_layerblendtimes[];
                 float locomotion_clip_times[]; vector4 locomotion_frames[];
                 float locomotion_clip_retimes[])
{
    removeindex(clip_names, i);
    removeindex(clip_times, i);
    removeindex(clip_weights, i);
    removeindex(clip_retimes, i);
    removeindex(clip_looped, i);

    if (len(clip_transformgroups))
    {
        removeindex(clip_transformgroups, i);
        removeindex(clip_layerids, i);
        removeindex(clip_layerblendstates, i);
        removeindex(clip_layerblendtimes, i);
    }

    removeindex(locomotion_clip_times, i);
    removeindex(locomotion_frames, i);
    removeindex(locomotion_clip_retimes, i);
}

#endif
