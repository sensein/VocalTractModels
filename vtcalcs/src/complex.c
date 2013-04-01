

/***************************************************************************
*									   *
*	File : complex.c						   *
*	Nore : Basic complex operations.				   *
*									   *
***************************************************************************/

#include	<math.h>
#include	<complex.h>

cmplx complex ( float x, float y )		/* (x, jy)  */
{
	cmplx	a;
	a.re = x;
	a.im = y;
	return ( a );
}

cmplx cmplx_add ( cmplx a, cmplx b )		/* a + b */
{
	return( complex( a.re+b.re, a.im+b.im ) );
}

cmplx cmplx_sub ( cmplx a, cmplx b )		/*  a - b  */
{
	return( complex( a.re-b.re, a.im-b.im ) );
}

cmplx cmplx_scl ( float x, cmplx a )		/*  x*a  */
{
	return( complex( x*a.re, x*a.im ) );
}

cmplx cmplx_mul ( cmplx a, cmplx b )		/*  a*b  */
{
	return( complex( a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re ) );
}

cmplx cmplx_div ( cmplx a, cmplx b )		/*  a/b  */
{
	float	d;
	d = b.re*b.re + b.im*b.im;
	return( complex( (a.re*b.re + a.im*b.im)/d, (a.im*b.re - a.re*b.im)/d ) );
}

cmplx cmplx_add_mul ( cmplx a, cmplx b, cmplx c )	/* a + b*c */
{
	return( complex( a.re+b.re*c.re-b.im*c.im, a.im+b.re*c.im+b.im*c.re ) );
}

cmplx cmplx_sub_mul ( cmplx a, cmplx b, cmplx c )	/* a - b*c */
{
	return( complex( a.re-b.re*c.re+b.im*c.im, a.im-b.re*c.im-b.im*c.re ) );
}

cmplx cmplx_inv ( cmplx a )			/* 1/a */
{
	float	d;
	d = a.re*a.re + a.im*a.im;
	return( complex( a.re/d, -a.im/d) );
}

float abs_cmplx ( cmplx a )			/* |a| */
{
	return( sqrt(a.re*a.re + a.im*a.im) );
}

float pha_cmplx ( cmplx a )			/* phase angle */
{
	return( atan2(a.im, a.re) );
}

