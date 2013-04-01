

/*****
*	File :	menu.h
*	Note :	various declarations.
*****/

#define MAX_FOOTNOTE_LEN  78
#define	MAX_MENU_ITEM_LEN 50
#define	RESP_OVERHEAD 2
#define	BkSp 8
#define CR 13
#define Esc 27
#define Ins 82
#define UP_ARROW 72
#define LEFT_ARROW 75
#define RIGHT_ARROW 77
#define DOWN_ARROW 80

typedef	struct	{int flag; int p; int p_old; int n; } menu_state;

	/* note:							*/
	/*	flag  = 0 at an initial entry, replaced by 1		*/
	/*	p     = an item selected ( 0, 1,..., or n-1)		*/
	/*	p_old = the old p					*/
	/*	n     = number of items in menu				*/

typedef	struct	{char mssg[MAX_MENU_ITEM_LEN];
		 char footnote[MAX_FOOTNOTE_LEN];
		 int  pos; int len;
		 int left; int top; int right; int bottom;} menu_items;

	/* note:							*/
	/*	mssg  = item message ( maximum 50 characters )		*/
	/*	footnote = If NULL, no display.				*/
	/*	pos   = position (in pixels) of the first letter	*/
	/*	len   = the length of the string in pixels		*/
	/*	(left, top, ...) = item box coordinates in pixels	*/


void	get_str    ( int y_pix, int x_pix, char *mssg,
		     int max_len, char *response );
void	menu_main  ( menu_state *state, menu_items *item );
void	menu_horiz ( int x_pix, int y_pix,
		     menu_state *state, menu_items *item, void **image_buf );
void	menu_vert  ( int x_pix, int y_pix,
		     menu_state *state, menu_items *item, void **image_buf );
