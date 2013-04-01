

/*************************************************************************
*	File : complex.h						 *
*	Note : prototypes for comlex operations.			 *
*************************************************************************/

typedef struct { float re; float im; } cmplx;

cmplx	complex ( float x, float y );			/*  x + jy	*/
cmplx	cmplx_add ( cmplx a, cmplx b );			/*  a + b	*/
cmplx	cmplx_sub ( cmplx a, cmplx b );			/*  a - b	*/
cmplx	cmplx_scl ( float x, cmplx a );			/*  xa		*/
cmplx	cmplx_mul ( cmplx a, cmplx b );			/*  ab		*/
cmplx	cmplx_div ( cmplx a, cmplx b );			/*  a/b		*/
cmplx	cmplx_add_mul ( cmplx a, cmplx b, cmplx c );	/*  a + bc	*/
cmplx	cmplx_sub_mul ( cmplx a, cmplx b, cmplx c );	/*  a - bc	*/
cmplx	cmplx_inv ( cmplx a );				/*   1/a	*/
float	abs_cmplx ( cmplx a );				/*   |a|	*/
float	pha_cmplx ( cmplx a );				/* phase angle  */
