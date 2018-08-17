/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *  	Side Effects Software Inc
 *  	123 Front Street West
 *  	Toronto, Ontario
 *  	Canada   M5V 3E7
 *  	416-504-9876
 *
 * NAME:    fuzzy_util.h
 */
#ifndef __fuzzy_util__
#define __fuzzy_util__

#include <math.h>

struct FuzzySet
{
    float samples[];
    float confidence;
    float minimum;
    float maximum;
}

///
/// Description: Samples from a finite 1D signal (linearly interpolating between discrete values)
/// 
/// Parameters:
///     * signal - signal to sample from
///     * x - location to sample
///     * minimum - start of signal
///     * maximum - end of signal
float sample_lerp(const float signal[]; const float x; const float minimum; const float maximum)
{
    if (x < minimum || x > maximum)
        return 0;

    float   array_ind = fit(x, minimum, maximum, 0, len(signal) - 1);
    int     ind_low   = int(floor(array_ind));
    int     ind_high  = int(ceil(array_ind));
     
    return fit(array_ind, ind_low, ind_high, signal[ind_low], signal[ind_high]);
}

///
/// Description: Creates a FuzzySet from a ramp (and metadata)
/// 
/// Parameters:
///     * ramp_basis - Ramp basis
///     * ramp_values - Ramp values
///     * ramp_positions - Ramp positions
///     * nsamples - Number of samples to take for the LUT
///	* confidence - degree of membership in the Fuzzy Set given a crisp input (from fuzzify())
///     * minimum - crisp minimum associated with the Fuzzy Set
///     * maximum - crisp maximum associated with the Fuzzy Set
FuzzySet fuzzy_create_set(const string ramp_basis[]; const float ramp_values[]; const float ramp_positions[]; const int nsamples; const float confidence; const float minimum; const float maximum)
{
    // Store the sampled membership function
    float samples[];
    resize(samples, nsamples);
    for(int i = 0; i < nsamples; ++i)
    {	
	float x = fit(i, 0, nsamples - 1, 0, 1);
	samples[i] = spline(ramp_basis, x, ramp_values, ramp_positions);
    }

    // Store the membership value
    return FuzzySet(samples, confidence, minimum, maximum);
}

///
/// Description: Creates a mirrored FuzzySet from a ramp (and metadata)
/// 
/// Parameters:
///     * ramp_basis - Ramp basis
///     * ramp_values - Ramp values
///     * ramp_positions - Ramp positions
///     * nsamples - Number of samples to take for the LUT
///	* confidence - degree of membership in the Fuzzy Set given a crisp input (from fuzzify())
///     * minimum - crisp minimum associated with the Fuzzy Set
///     * maximum - crisp maximum associated with the Fuzzy Set
FuzzySet fuzzy_create_set_mirror(const string ramp_basis[]; const float ramp_values[]; const float ramp_positions[]; const int nsamples; const float confidence; const float minimum; const float maximum)
{
    // Store the sampled membership function
    float samples[];
    resize(samples, nsamples);
    for(int i = 0; i < nsamples; ++i)
    {	
	float x = fit(i, 0, nsamples - 1, 0, 1);
	samples[i] = spline(ramp_basis, 1 - x, ramp_values, ramp_positions);
    }

    // Store the membership value
    return FuzzySet(samples, confidence, minimum, maximum);
}

///
/// Description: Creates a FuzzySet that spans range A from a FuzzySet that spans range B, and returns 
///		 the FuzzySet masked against range B
/// 
/// Parameters:
///     * ramp_basis - Ramp basis
///     * ramp_values - Ramp values
///     * ramp_positions - Ramp positions
///     * nsamples - Number of samples to take for the LUT
///	* confidence - degree of membership in the Fuzzy Set given a crisp input (from fuzzify())
///     * minimum - crisp minimum associated with the Fuzzy Set
///     * maximum - crisp maximum associated with the Fuzzy Set
FuzzySet fuzzy_create_set_tform(const string ramp_basis[]; const float ramp_values[]; const float ramp_positions[]; const int nsamples; const float confidence; 
				const float min_old; const float max_old; const float min_new; const float max_new)
{
    // Edge case when the output range is invalid
    if (max_new - min_new < M_TOLERANCE)
	return fuzzy_create_set({ "linear" }, { 0 }, { 0 }, nsamples, confidence, 0, 1);
    
    float positions_tform[];
    float samples[];

    resize(positions_tform, len(ramp_positions));
    resize(samples, nsamples);

    // Re-map the basis to the new min/max 
    for (int i = 0; i < len(ramp_positions); ++i)
    {
	positions_tform[i] = fit(ramp_positions[i], 0, 1, min_new, max_new); 
    }

    for(int i = 0; i < nsamples; ++i)
    {	
	float x = fit(i, 0, nsamples - 1, min_old, max_old);
	samples[i] = spline(ramp_basis, x, ramp_values, positions_tform);
    }

    // Store the membership value
    return FuzzySet(samples, confidence, min_old, max_old);
}

///
/// Description: Aggregates two membership functions that cover two distinct crisp ranges using "max" aggregation (maximum value of each sample) 
///		 The functions are zero-padded and linearly interpolated so that the aggregate covers the union of the ranges. The smallest delta 
///		 between the two functions is preserved.
/// 
/// Parameters:
///     * aggregated_membership - Aggregated membership function (modified by function)
///     * min_a	- crisp minimum associated with the aggregated membership function (modified by function)
///     * max_a - crisp maximum associated with the aggregated membership function (modified by function)
///	* membership - membership function to aggregate into aggregated_membership
///     * minimum - crisp minimum associated with the membership function to aggregate
///     * maximum - crisp maximum associated with the membership function to aggregate
void fuzzy_max_aggregate(float aggregated_membership[]; float min_a; float max_a; const float membership[]; const float min_b; const float max_b; const float confidence)
{
    float   d_x = -1, min_new = min_a, max_new = max_a;
    float   membership_up[]; 
    int	    array_len = -1;

    if (len(aggregated_membership) > 0)
	d_x = (max_a - min_a) / len(aggregated_membership);
    if (len(membership) > 0)
    {
	float d_b = (max_b - min_b) / len(membership);
	if (d_x == -1)
	    d_x = d_b;
	else
	    d_x = min(d_x, d_b);
    }

    // If d_x == -1, then there is no valid input
    if (d_x == -1)
	return;

    min_new = min(min_a, min_b);
    max_new = max(max_a, max_b);
    array_len = int(ceil((max_new - min_new) / d_x));
    resize(membership_up, array_len);
    
    for(int i = 0; i < array_len; ++i)
    {
	float x = fit(i, 0, array_len - 1, min_new, max_new);
	float dom = min(confidence, sample_lerp(membership, x, min_b, max_b));

	membership_up[i] = sample_lerp(aggregated_membership, x, min_a, max_a);
	membership_up[i] = max(membership_up[i], dom);
    }

    // Set the associated metadata
    aggregated_membership = membership_up;
    min_a = min_new;
    max_a = max_new;
}

///
/// Description: Aggregates two membership functions that cover two distinct crisp ranges using "sum" aggregation (sum the values at each sample) 
///		 The functions are zero-padded and linearly interpolated so that the aggregate covers the union of the ranges. The smallest delta 
///		 between the two functions is preserved.
/// 
/// Parameters:
///     * aggregated_membership - Aggregated membership function (modified by function)
///     * min_a	- crisp minimum associated with the aggregated membership function (modified by function)
///     * max_a - crisp maximum associated with the aggregated membership function (modified by function)
///	* membership - membership function to aggregate into aggregated_membership
///     * minimum - crisp minimum associated with the membership function to aggregate
///     * maximum - crisp maximum associated with the membership function to aggregate
void fuzzy_sum_aggregate(float aggregated_membership[]; float min_a; float max_a; const float membership[]; const float min_b; const float max_b; const float confidence)
{
    float   d_x = -1, min_new = min_a, max_new = max_a;
    float   membership_up[]; 
    int	    array_len = -1;

    if (len(aggregated_membership) > 0)
	d_x = (max_a - min_a) / len(aggregated_membership);
    if (len(membership) > 0)
    {
	float d_b = (max_b - min_b) / len(membership);
	if (d_x == -1)
	    d_x = d_b;
	else
	    d_x = min(d_x, d_b);
    }

    // If d_x == -1, then there is no valid input
    if (d_x == -1)
	return;

    min_new = min(min_a, min_b);
    max_new = max(max_a, max_b);
    array_len = int(ceil((max_new - min_new) / d_x));
    resize(membership_up, array_len);
    
    for(int i = 0; i < array_len; ++i)
    {
	float x = fit(i, 0, array_len - 1, min_new, max_new);
	float dom = min(confidence, sample_lerp(membership, x, min_b, max_b));

	membership_up[i] = sample_lerp(aggregated_membership, x, min_a, max_a);
	membership_up[i] += dom;
    }

    // Set the associated metadata
    aggregated_membership = membership_up;
    min_a = min_new;
    max_a = max_new;
}

///
/// Description: Aggregates two membership functions that cover the same crisp range (using "max" aggregation) 
///
/// Parameters:
///     * aggregated_membership - Aggregated membership function (modified by function)
///	* membership - membership function to aggregate into aggregated_membership
///	* confidence - confidence value of membership to clip membership with
void fuzzy_max_aggregate(float aggregated_membership[]; const float membership[]; const float confidence)
{
    // Up sample whichever signal has a smaller sampling rate
    int	  array_len = max(len(aggregated_membership), len(membership));
    float membership_up[] = resample_linear(membership, array_len);
    aggregated_membership = resample_linear(aggregated_membership, array_len);

    for(int i = 0; i < array_len; ++i)
    {
	float dom = min(confidence, membership_up[i]);
	aggregated_membership[i] = max(aggregated_membership[i], dom);
    }
}

///
/// Description: Aggregates two membership functions that cover the same crisp range (using "sum" aggregation) 
///
/// Parameters:
///     * aggregated_membership - Aggregated membership function (modified by function)
///	* membership - membership function to aggregate into aggregated_membership
///	* confidence - confidence value of membership to clip membership with
void fuzzy_sum_aggregate(float aggregated_membership[]; const float membership[]; const float confidence)
{
    // Up sample whichever signal has a smaller sampling rate
    int	  array_len = max(len(aggregated_membership), len(membership));
    float membership_up[] = resample_linear(membership, array_len);
    aggregated_membership = resample_linear(aggregated_membership, array_len);

    for(int i = 0; i < array_len; ++i)
    {
	float dom = min(confidence, membership_up[i]);
	aggregated_membership[i] += dom;
    }
}

///
/// Description: Fuzzifies an array of inputs that uniformly covers the membership function into a fuzzy value. 
///
/// Parameters:
///     * ramp_basis - ramp basis
///	* ramp_values - ramp values
///	* ramp_positions - ramp positions
///	* crisp_range - input array where each element is a uniformly spaced sample of the membership function
///	* min_value - crisp minimum associated with the membership function
///	* max_value - crisp maximum associated with the membership function
float
fuzzify_range(const string ramp_basis[]; const float ramp_values[];
              const float ramp_positions[]; const string imp_ramp_basis[];
              const float imp_ramp_values[]; const float imp_ramp_positions[];
              const float crisp_range[]; const float min_value;
              const float max_value)
{
    float weight = 0, total = 0;
    int samples = len(crisp_range);

    for (int i = 0; i < samples; ++i)
    {
        float x = fit(i, 0, samples - 1, 0, 1);
        float w =
            spline(imp_ramp_basis, x, imp_ramp_values, imp_ramp_positions);
        weight += crisp_range[i] * w;
        total += w;
    }

    if (total == 0)
        return 0;

    float crisp_val = weight / total;
    return spline(ramp_basis, fit(crisp_val, min_value, max_value, 0, 1),
                  ramp_values, ramp_positions);
}
#endif
