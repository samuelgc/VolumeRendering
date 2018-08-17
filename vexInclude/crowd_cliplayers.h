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
 * NAME:    crowd_cliplayers.h
 *
 */

#ifndef __crowd_cliplayers_h__
#define __crowd_cliplayers_h__

#define BLEND_INTERPOLATE 0
#define BLEND_ADDITIVE 1

void
crowdcliplayer_insert(const int i; int layer_modes[]; int layer_parents[];
                      float layer_weights[])
{
    insert(layer_modes, i, BLEND_INTERPOLATE);
    insert(layer_parents, i, -1);
    insert(layer_weights, i, 0.0);
}

void
crowdcliplayer_remove(const int i; int layer_modes[]; int layer_parents[];
                      float layer_weights[])
{
    removeindex(layer_modes, i);
    removeindex(layer_parents, i);
    removeindex(layer_weights, i);
}

/// Returns the index of the layer for blending between two states.
int
crowdcliplayer_blendindex(const int next_state_start)
{
    return next_state_start - 1;
}

/// Update the index where layers from the next state will start.
void
crowdcliplayer_updatenextstate(const int last_layer_idx; int next_state_start)
{
    // Reserve room for the layer that blends between the two states.
    next_state_start = last_layer_idx + 2;
}

/// Returns whether the specified layer is the base layer of a state.
int
crowdcliplayer_isbaselayer(const int layer_idx; const int next_state_start)
{
    return layer_idx == 0 || layer_idx == next_state_start;
}

#endif
