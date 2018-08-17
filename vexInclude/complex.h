#ifndef __complex_h
#define __complex_h

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
 * NAME:	complex.h ( Complex number functions, VEX )
 *
 * COMMENTS:	Complex number functions, a subset of the C99 complex.h
 *		functions, with a few additions for assuming only real or
 *		imaginary input for thesecond parameter (e.g. cmultimag).
 */

// Complex number defintion.
struct complex
{
    float real, imag;
}

// Complex number constructor function.
complex CMPLX(float real; float imag)
{
    return complex(real, imag);
}

// Return the real component of a complex number.
float creal(const complex c)
{
    return c.real;
}

// Return the imaginary component of a complex number.
float cimag(const complex c)
{
    return c.imag;
}

// Multiply two complex numbers.
complex cmult(const complex a; const complex b)
{
    return complex(a.real * b.real - a.imag * b.imag,
		   a.real * b.imag + a.imag * b.real);
}

// Multiply a complex number by a real scalar.
complex cmultreal(const complex a; const float b)
{
    return complex(b * a.real, b * a.imag);
}

// Multiply a complex number by an imaginary scalar.
complex cmultimag(const complex a; const float b)
{
    return complex(-a.imag * b, a.real * b);
}

// Add two complex numbers.
complex cadd(const complex a; const complex b)
{
    return complex(a.real + b.real, a.imag + b.imag);
}

// Subtract the second complex number from the first.
complex csub(const complex a; const complex b)
{
    return complex(a.real - b.real, a.imag - b.imag);
}

// Compute the complex conjugate of a complex number.
complex conj(const complex a)
{
    return complex(a.real, -a.imag);
}

// Compute the absolute value of a complex number.
float cabs(const complex a)
{
    return sqrt(a.real * a.real + a.imag * a.imag);
}

// Compute the complex exponential of a complex number.
complex cexp(const complex a)
{
    float rexp = exp(a.real);
    return complex(rexp * cos(a.imag), rexp * sin(a.imag));
}

// Compute the complex exponential of a real scalar,
complex cexpreal(const float b)
{
    return complex(exp(b), 0);
}

// Compute the complex exponential of an imaginary scalar,
// i.e. Euler's formula: e^(i x) = cos(x) + i sin(x)
complex cexpimag(const float x)
{
    return complex(cos(x), sin(x));
}

#endif
