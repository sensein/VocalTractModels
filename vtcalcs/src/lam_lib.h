

/********
*	File :	lam-lib.h
*	Note :	The definition of structures and prototypes.
********/

typedef	struct{ int x, y;} int2D;
typedef	struct{ float x, y;} float2D;

void	read_model_spec( void );
void	convert_scale( void );
void	semi_polar( void );
void	lam( float *para);
void	sagittal_to_area( int *ns, area_function *af );
void	appro_area_function (int ns1, area_function *A1,
			     int ns2, area_function *A2);

void	plot_semi_polar( int nvp );
void	plot_lam( int nvp );
void	deplot_lam( void );
