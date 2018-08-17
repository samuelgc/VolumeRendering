/*
 * PROPRIETARY INFORMATION.  This software is proprietary to
 * Side Effects Software Inc., and is not to be reproduced,
 * transmitted, or disclosed in any way without written permission.
 *
 * Produced by:
 *	Side Effects Software Inc
 *	123 Front Street West, Suite 1401
 *	Toronto, Ontario
 *	Canada   M5J 2M2
 *	416-504-9876
 *
 * NAME:	assert.h (VEX Library, C++)
 *
 * COMMENTS:
 */

#ifndef __vex_assert__
#define __vex_assert__

// VEX Assertions.  The assert_enabled() function checks the state of the
// HOUDINI_VEX_ASSERT environment variable, so assertions will only be printed
// if the variable is set.

#define assert(EXPR)	\
    if (assert_enabled()) { \
	if (!(EXPR)) print_once(sprintf('VEX Assertion Failed %s:%d - (%s)\n', \
		__FILE__, __LINE__, #EXPR)); \
    }

#endif
