/****
*	File : mydir.h
*	Note : directory specifications.
***/

/* For windows */
#ifdef WINDOWS_COMPILE
	char	*NTAFpath    = "data\\area\\Nt.are";
/*	char	*VTAFdir     = "data\\area\\"; */

/*	char	*SIGdir      = "data\\sig\\"; 
	char	*tempSIGpath = "data\\sig\\tempsig.sig";
	char	*tempGLTpath = "data\\sig\\tempglt.sig";
*/

	char	*MODELSPECpath = "data\\specdata\\pb1_spec.dat";
	char	*RADIMPpath    = "data\\specdata\\rad_imp.dat";
#else
	char	*NTAFpath    = "data//area//Nt.are";
/*	char	*VTAFdir     = "data//area//";

	char	*SIGdir      = "data//sig//";
	char	*tempSIGpath = "data//sig//tempsig.sig";
	char	*tempGLTpath = "data//sig//tempglt.sig";
*/
	char	*MODELSPECpath = "data//specdata//pb1_spec.dat";
	char	*RADIMPpath    = "data//specdata//rad_imp.dat";
#endif
