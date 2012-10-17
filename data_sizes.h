/******************************************************************************/
/*                                                                            */
/*                          sizes of data structures                          */
/*                                                                            */
/******************************************************************************/

#if DIM == 2

#define MAXVERT         NS*NS
#define MAXNODE         MAXVERT
#define MAXLINK         6*MAXNODE
#define MAXELEM         2*(NS-1)*(NS-1)
#define MAXFACE         (NS-1)*(3*NS-1)
#define MAXEDGE         (NS-1)*(3*NS-1)
#define MAXBPOINT       60*(NS-1)

#if DATA_S & F_LINK_TO_FACES
#define MAXFLINK        4*MAXFACE
#else
#define MAXFLINK        1
#endif

#if DATA_S & N_LINK_TO_FACES
#define MAXNFLINK       12*MAXNODE         
#else
#define MAXNFLINK       1         
#endif

#if DATA_S & N_LINK_TO_ELEMENTS 
#define MAXNELINK       6*MAXNODE
#else
#define MAXNELINK       1
#endif

#if DATA_S & F_LINK_TO_NODES
#define MAXFNLINK       4*MAXFACE         
#else
#define MAXFNLINK       1         
#endif

#if DATA_S & SPECIAL_NODES_AND_FACES
#define MAXSNODE        4*NS
#define MAXSFACE        4*NS
#else
#define MAXSNODE        1
#define MAXSFACE        1
#endif

#if E_DATA & E_E_NEIGHBOURS
#define MAXELINK        12*MAXELEM
#else
#define MAXELINK        1
#endif

/*----------------------------------------------------------------------------*/

#else  /*  if DIM == 3  */

#define MAXSNODE        1
#define MAXSFACE        1

#define MAXVERT         NS*NS*NS
#define MAXNODE         MAXVERT
#define MAXLINK         14*MAXNODE
#define MAXELEM         6*(NS-1)*(NS-1)*(NS-1)
#define MAXFACE         6*(NS-1)*(NS-1)*(2*NS-1)
#define MAXEDGE         (NS-1)*(7*NS*NS-5*NS+1) 
#define MAXBPOINT       4*(NS-1)*(3*NS-1)

#if DATA_S & F_LINK_TO_FACES
#define MAXFLINK        6*MAXFACE
#else
#define MAXFLINK        1
#endif

#if DATA_S & N_LINK_TO_FACES
#define MAXNFLINK       60*MAXNODE         
#else
#define MAXNFLINK       1         
#endif

#if DATA_S & N_LINK_TO_ELEMENTS 
#define MAXNELINK       24*MAXNODE
#else
#define MAXNELINK       1
#endif

#if DATA_S & F_LINK_TO_NODES
#define MAXFNLINK       5*MAXFACE         
#else
#define MAXFNLINK       1         
#endif

#if E_DATA & E_E_NEIGHBOURS
#define MAXELINK        92*MAXELEM
#else
#define MAXELINK        1
#endif

#endif

/*----------------------------------------------------------------------------*/

#define MAXLEVEL        31

#if INCLUDE_UMFPACK == YES
#define MAX_ROW         (10*MAXNODE)
#define MAX_ENT         (10*(MAXNODE+MAXLINK))
#else
#define MAX_ROW         1
#define MAX_ENT         1
#endif

#if DATA_STR & LG_DATA
#define MAXNLGLINK      10000
#define MAXLGNLINK      10000
#define MAXLGFLINK      100000
#define MAXLGLGLINK     10000
#define MAXFLGLINK      100000
#define MAXLGDATA       10000
#else
#define MAXNLGLINK      1
#define MAXLGNLINK      1
#define MAXLGFLINK      1
#define MAXLGLGLINK     1
#define MAXFLGLINK      1
#define MAXLGDATA       1
#endif
