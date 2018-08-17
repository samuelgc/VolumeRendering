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
 * NAME:	shaderlayer.h ( VOP Library Functions, VEX )
 *
 * COMMENTS:	Support functions for Shader Layers.
 *
 */

#ifndef __shaderlayer__
#define __shaderlayer__

struct ShaderExports
{
    string	names_f[];
    float	values_f[];
    string	names_v[];
    vector	values_v[];
    string	names_v4[];
    vector4	values_v4[];
}

struct ShaderLayer
{
    bsdf	    F;
    vector	    Of;
    vector	    Ce;
    vector	    P;
    vector	    N;
    float	    layeralpha;
    float	    thickness;
    vector	    absorption;
    float	    masks[];
    ShaderExports   exports;
}

void init_layerexports(export ShaderExports exports)
{
    exports.names_f = {};
    exports.values_f = {};
    exports.names_v = {};
    exports.values_v = {};
    exports.names_v4 = {};
    exports.values_v4 = {};
}

void init_layer(export ShaderLayer layer)
{
    layer.F = bsdf();
    layer.Of = {1,1,1};
    layer.Ce = {0,0,0};
    layer.P = P;
    layer.N = normalize(N);
    layer.layeralpha = 1.0;
    layer.masks = {};
    init_layerexports(layer.exports);
}

// Color and alpha expression for each composite operation.
#define COMP_AOVERB_COLOR(A, Aa, B, Ba)		Aa*A+(1-Aa)*Ba*B
#define COMP_AOVERB_ALPHA(Aa, Ba)		Aa+(1-Aa)*Ba

#define COMP_AINSIDEB_COLOR(A, Aa, B, Ba)	Aa*A*Ba
#define COMP_AINSIDEB_ALPHA(Aa, Ba)		Aa*Ba

#define COMP_AOUTSIDEB_COLOR(A, Aa, B, Ba)	Aa*A*(1-Ba)
#define COMP_AOUTSIDEB_ALPHA(Aa, Ba)		Aa*(1-Ba)

#define COMP_AATOPB_COLOR(A, Aa, B, Ba)		Aa*A*Ba + Ba*B*(1-Aa)
#define COMP_AATOPB_ALPHA(Aa, Ba)		Ba

#define COMP_AXORB_COLOR(A, Aa, B, Ba)		Aa*A*(1-Ba) + Ba*B*(1-Aa)
#define COMP_AXORB_ALPHA(Aa, Ba)		Aa+Ba-2*(Aa*Ba)

#define CREATE_ARRAY_COMPOP_TYPE(NAME, VALUE_TYPE, COMPOP) \
void NAME(const string names_a[]; const VALUE_TYPE vals_a[]; const float Aa; \
	  const string names_b[]; const VALUE_TYPE vals_b[]; const float Ba; \
	  string out_names[]; VALUE_TYPE out_values[]) \
{ \
    VALUE_TYPE val_a, val_b; \
    int done_in_b[]; \
 \
    /* loop over exports in A */ \
    foreach(int index_a; string name_a; names_a) \
    { \
	int index_b = index_a; \
	if (names_b[index_b] != name_a) \
	    index_b = find(names_b, name_a); \
 \
	/* if this export is also in B, comp with B's value */ \
	/* if not, comp with black */ \
	int found = index_b >= 0; \
	val_b = select(found, vals_b[index_b], 0.0); \
 \
	append(out_names, name_a); \
	append(out_values, COMP_AOVERB_COLOR(vals_a[index_a], Aa, val_b, Ba)); \
 \
	/* record whether we've handled this B export */ \
	done_in_b[index_b] = found; \
    } \
 \
    /* loop over exports in B */ \
    foreach(int index_b; string name_b; names_b) \
    { \
	/* Skip this iteam if we've already handled it when we looped through */ \
	/* A's exports. */ \
	if(done_in_b[index_b]) \
	    continue; \
 \
	/* If we get to here, we know that this export doesn't exist */ \
	/* in A, so just comp with black (just multiply with alpha) */ \
	val_a = 0.0; \
	val_b = vals_b[index_b]; \
	append(out_names, name_b); \
	append(out_values, COMPOP(val_a, Aa, val_b, Ba)); \
    } \
}

#define CREATE_ARRAY_COMPOP(VALUE_TYPE) \
CREATE_ARRAY_COMPOP_TYPE(comp_export_arrays_aoverb, VALUE_TYPE, COMP_AOVERB_COLOR) \
CREATE_ARRAY_COMPOP_TYPE(comp_export_arrays_ainsideb, VALUE_TYPE, COMP_AINSIDEB_COLOR) \
CREATE_ARRAY_COMPOP_TYPE(comp_export_arrays_aoutsideb, VALUE_TYPE, COMP_AOUTSIDEB_COLOR) \
CREATE_ARRAY_COMPOP_TYPE(comp_export_arrays_aatopb, VALUE_TYPE, COMP_AATOPB_COLOR) \
CREATE_ARRAY_COMPOP_TYPE(comp_export_arrays_axorb, VALUE_TYPE, COMP_AXORB_COLOR)

CREATE_ARRAY_COMPOP(float)
CREATE_ARRAY_COMPOP(vector)
CREATE_ARRAY_COMPOP(vector4)

#define CREATE_EXPORT_COMPOP(COMPOP) \
void composite_exports_##COMPOP(const ShaderLayer A; const float Aa; const ShaderLayer B; const float Ba; ShaderLayer out) \
{ \
    out.exports.names_f = {}; \
    out.exports.values_f = {}; \
    out.exports.names_v = {}; \
    out.exports.values_v = {}; \
    out.exports.names_v4 = {}; \
    out.exports.values_v4 = {}; \
 \
    comp_export_arrays_##COMPOP(A.exports.names_f, A.exports.values_f, Aa, \
				B.exports.names_f, B.exports.values_f, Ba, \
				out.exports.names_f, out.exports.values_f); \
    comp_export_arrays_##COMPOP(A.exports.names_v, A.exports.values_v, Aa, \
				B.exports.names_v, B.exports.values_v, Ba, \
				out.exports.names_v, out.exports.values_v); \
    comp_export_arrays_##COMPOP(A.exports.names_v4, A.exports.values_v4, Aa, \
				B.exports.names_v4, B.exports.values_v4, Ba, \
				out.exports.names_v4, out.exports.values_v4); \
}

CREATE_EXPORT_COMPOP(aoverb)
CREATE_EXPORT_COMPOP(ainsideb)
CREATE_EXPORT_COMPOP(aoutsideb)
CREATE_EXPORT_COMPOP(aatopb)
CREATE_EXPORT_COMPOP(axorb)

// Define a composite operation applied to all layer members.
// The top variant lets the caller specify alpha values
// while the bottom variant users the intrinsic alpha values of each layer
#define CREATE_COMPOP(NAME, COLOROP, ALPHAOP) \
ShaderLayer composite_##NAME(const ShaderLayer A; const float Aa; \
	         const ShaderLayer B; const float Ba) \
{ \
    ShaderLayer C; \
    C.masks = {}; \
    C.F = COLOROP(A.F, Aa, B.F, Ba); \
    C.Of = COLOROP(A.Of, Aa, B.Of, Ba); \
    C.Ce = COLOROP(A.Ce, Aa, B.Ce, Ba); \
    C.P = COLOROP(A.P, Aa, B.P, Ba); \
    C.N = normalize(COLOROP(A.N, Aa, B.N, Ba)); \
    C.layeralpha = ALPHAOP(Aa, Ba); \
    composite_exports_##NAME(A, Aa, B, Ba, C); \
 \
    int nummasks_a = len(A.masks); \
    int nummasks_b = len(B.masks); \
    int maxmasks = max(nummasks_a, nummasks_b); \
 \
    for(int i=0; i<maxmasks; i++) \
    { \
	float mask_a = select(i<nummasks_a, A.masks[i], 0.0); \
	float mask_b = select(i<nummasks_b, B.masks[i], 0.0); \
	C.masks[i] = COLOROP(mask_a, Aa, mask_b, Ba); \
    } \
    return C; \
} \
ShaderLayer composite_##NAME(const ShaderLayer A; const ShaderLayer B) \
{ \
    return composite_##NAME(A, A.layeralpha, B, B.layeralpha); \
}

CREATE_COMPOP(aoverb, COMP_AOVERB_COLOR, COMP_AOVERB_ALPHA)
CREATE_COMPOP(ainsideb, COMP_AINSIDEB_COLOR, COMP_AINSIDEB_ALPHA)
CREATE_COMPOP(aoutsideb, COMP_AOUTSIDEB_COLOR, COMP_AOUTSIDEB_ALPHA)
CREATE_COMPOP(aatopb, COMP_AATOPB_COLOR, COMP_AATOPB_ALPHA)
CREATE_COMPOP(axorb, COMP_AXORB_COLOR, COMP_AXORB_ALPHA)

#define CREATE_SETLAYEREXPORT(TYPE, POSTFIX) \
void set_layer_export(ShaderLayer layer; string name; TYPE value) \
{ \
    int index = find(layer.exports.names_##POSTFIX, name); \
    if(index < 0) \
	index = len(layer.exports.names_##POSTFIX); \
    layer.exports.names_##POSTFIX[index] = name; \
    layer.exports.values_##POSTFIX[index] = value; \
}

CREATE_SETLAYEREXPORT(float, f)
CREATE_SETLAYEREXPORT(vector, v)
CREATE_SETLAYEREXPORT(vector4, v4)

#endif
