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
 * NAME:	pop.h (VEX)
 *
 * COMMENTS:	Contains useful POP definitions
 */

#ifndef __pop_h__
#define __pop_h__

#define	PSTATE_PRIMARY	0x01	// Particle not birthed from another particle
#define PSTATE_DYING	0x02	// Particle will die before the next frame
#define PSTATE_STOPPED	0x04	// Particle is stopped (motion doesn't occur)
#define PSTATE_COLLIDE	0x08	// Particle collided with something
#define PSTATE_STUCK	0x10	// Particle is stuck to geometry
#define PSTATE_ISRBD	0x20	// Particle is part of a rigid body simulation
#define PSTATE_ISACTIVE	0x40	// Particle is an RBD active object
#define PSTATE_ISOVERRIDE	0x80	// Particle motion is done by a CHOP

#define PSTATE_DONTMOVE	(PSTATE_DYING|PSTATE_STOPPED|PSTATE_STUCK)

#endif
