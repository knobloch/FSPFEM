/******************************************************************************/
/*                                                                            */
/*                     specifications of data structures                      */
/*                                                                            */
/******************************************************************************/

#if (U_SPACE == P1_NC || U_SPACE == Q1ROT) && (P_SPACE == P0)

#define DATA_S          10
#define N_DATA           0
#define E_DATA         161
#define T_FOR_BC        2 

#if DATA_STR & KORN_MATRIX
#define F_DATA          40
#else
#define F_DATA          33
#endif

#elif (U_SPACE == P1_MOD) && (P_SPACE == P0 || P_SPACE == P1_DISC)

#define DATA_S          42
#define N_DATA           0
#define F_DATA         104
#define T_FOR_BC        2
#if P_SPACE == P0
#define E_DATA          33
#elif P_SPACE == P1_DISC
#define E_DATA          68
#endif

#elif (U_SPACE == P1C) && (P_SPACE == P0)

#define DATA_S          8
#define N_DATA          49
#define F_DATA          0
#define E_DATA          (2205|E_E_FNEIGHBOURS)
#define T_FOR_BC        2
#define USE_DIV_STAB    YES

#elif ((U_SPACE == P1C_FBUB) || (U_SPACE == P1C_NEW_FBUB)) && (P_SPACE == P0)

#define DATA_S          15
#define N_DATA          49
#define F_DATA          255
#define E_DATA          29
#if DATA_STR & REDUCED
#define T_FOR_BC      34
#else
#define T_FOR_BC      2
#endif

#elif (U_SPACE == P1C_ELBUB) && (P_SPACE == P1C)

#define DATA_S          8
#define N_DATA          27
#define F_DATA          16
#define E_DATA          11
#define T_FOR_BC        2

#elif (U_SPACE == P2C) && (P_SPACE == P0)

#define DATA_S          15
#define N_DATA          593
#define F_DATA          289
#define E_DATA          169
#define T_FOR_BC        2

#elif ((U_SPACE == P2C) && (P_SPACE == P1C || P_SPACE == P1C_P2L)) || ((U_SPACE == IP2C) && (P_SPACE == IP1C)) && (MOVING_BOUNDARY == NO)

#if P_SPACE == P1C || P_SPACE == IP1C
#define DATA_S          15
#elif P_SPACE == P1C_P2L
#define DATA_S          207
#endif
#define T_FOR_BC        2

#if DATA_STR & KORN_MATRIX
#if U_SPACE == P2C
#define F_DATA         1576
#elif U_SPACE == IP2C
#define F_DATA         5672
#endif
#else
#if U_SPACE == P2C
#define F_DATA         1329
#elif U_SPACE == IP2C
#define F_DATA         5425
#endif
#endif

#if (LOW_AND_HIGH == YES) && (U_SPACE_LOW == P1_NC) && (P_SPACE_LOW == P0)

#define E_DATA          161
#if DATA_STR & KORN_MATRIX
#define N_DATA          958
#else
#define N_DATA          891
#endif

#elif (LOW_AND_HIGH == YES) && (U_SPACE_LOW == P1C) && (P_SPACE_LOW == P0)

#define E_DATA          (2205|E_E_FNEIGHBOURS)
#if DATA_STR & KORN_MATRIX
#define N_DATA          958
#else
#define N_DATA          891
#endif
#define USE_DIV_STAB_LOW YES

#else

#define E_DATA          0
#if DATA_STR & KORN_MATRIX
#define N_DATA          446
#else
#define N_DATA          379
#endif

#endif

#elif ((U_SPACE == P2C) && (P_SPACE == P1C)) || ((U_SPACE == IP2C) && (P_SPACE == IP1C)) && (MOVING_BOUNDARY == YES)

#define DATA_S          15
#define T_FOR_BC        2
#define PROBLEM         1 
#define E_DATA			0
//
#if DATA_STR & KORN_MATRIX
#if DATA_STR & BD_DATA
#define N_DATA         2559 //958
#if U_SPACE == P2C
#define F_DATA         10025 //1576
#elif U_SPACE == IP2C
#define F_DATA         14121 //5672
#endif
#else
#define N_DATA         511 //958
#if U_SPACE == P2C
#define F_DATA         1833 //1576
#elif U_SPACE == IP2C
#define F_DATA         5929 //5672
#endif
#endif
#else
#define N_DATA          123 //891
#if U_SPACE == P2C
#define F_DATA         305 //1329
#elif U_SPACE == IP2C
#define F_DATA         4401 //5425
#endif
#endif

#elif (U_SPACE == P2C) && (P_SPACE == P1_DISC)

#define DATA_S          15
#define N_DATA          337
#define F_DATA          1313
#define E_DATA          12420
#define T_FOR_BC        2

#elif (U_SPACE == P2C_ELBUB) && (P_SPACE == P0)

#define DATA_S          15
#define N_DATA          593
#define F_DATA          289
#define E_DATA          4027
#define T_FOR_BC        2

#elif ((U_SPACE == P2C_ELBUB) && (P_SPACE == P1C || P_SPACE == P1C_P2L)) || ((U_SPACE == IP2C_ELBUB) && (P_SPACE == IP1C))

#if P_SPACE == P1C || P_SPACE == IP1C
#define DATA_S          15
#elif P_SPACE == P1C_P2L
#define DATA_S          207
#endif
#define E_DATA         3866
#define T_FOR_BC        2
#define N_DATA          891
#if U_SPACE == P2C_ELBUB
#define F_DATA         1329
#elif U_SPACE == IP2C_ELBUB
#define F_DATA         5425
#endif

#elif ((U_SPACE == P2C_ELBUB) && (P_SPACE == P1_DISC)) || ((U_SPACE == IP2C_ELBUB) && (P_SPACE == IP1_DISC))

#define DATA_S          15
#define T_FOR_BC        2
#if U_SPACE == P2C_ELBUB
#define F_DATA          289
#elif U_SPACE == IP2C_ELBUB
#define F_DATA         4385
#endif

#if (LOW_AND_HIGH == YES) && (U_SPACE_LOW == P1_NC) && (P_SPACE_LOW == P0)
#define N_DATA          593
#define E_DATA          16311
#elif (LOW_AND_HIGH == YES) && (U_SPACE_LOW == P1C) && (P_SPACE_LOW == P0)
#define N_DATA          625
#define E_DATA          (16319|E_E_FNEIGHBOURS)
#define USE_DIV_STAB_LOW YES
#else
#define N_DATA          593
#define E_DATA          16310
#endif

#endif

#if P_SPACE == NONE

#if U_SPACE==P1_NC

#define DATA_S          10
#define N_DATA          0
#define F_DATA          17
#define E_DATA          0
#define T_FOR_BC        2

#elif U_SPACE==P1_MOD

#define DATA_S          10
#define N_DATA          0
#define F_DATA          56
#define E_DATA          0
#define T_FOR_BC        2

#elif (U_SPACE == P1C) || (U_SPACE == Q1C) || (U_SPACE == MINI_L_DIV_FR)

#if U_STRUCTURE == SCALAR
#define N_DATA          9
#elif U_STRUCTURE == VECTOR
#define N_DATA          17
#endif
#if PERIODIC_BC == YES
#define DATA_S          136
#else
#define DATA_S          520
#endif
#define F_DATA          16
#define E_DATA          (2|E_TAU|E_TAU_SOLD)
#define T_FOR_BC        2

#elif (U_SPACE == P1C_ELBUB) && (U_STRUCTURE == SCALAR)

#define N_DATA          9
#define DATA_S          8
#define F_DATA          0 
#define E_DATA          3585
#define T_FOR_BC        2

#elif (U_SPACE == Q2C) && (U_STRUCTURE == SCALAR) && (PERIODIC_BC == NO)

#define DATA_S          15
#define N_DATA          9
#define F_DATA          17
#define E_DATA          3856
#define T_FOR_BC        2

#elif (U_SPACE == GP1C) || (U_SPACE == GP1X3C) || (U_SPACE == GP1C_ELBUB) || (U_SPACE == GQ1C) || (U_SPACE == GQ1X4C) || (U_SPACE == GQ1C_ELBUB) || (U_SPACE == GP2C) || (U_SPACE == GP2X3C) || (U_SPACE == GP2C_3ELBUB) || (U_SPACE == GP2C_6ELBUB) || (U_SPACE == GQ2C) || (U_SPACE == GQ2X4C) || (U_SPACE == GQ2B3C) || (U_SPACE == GQ2C_2ELBUB) || (U_SPACE == GQ2C_3ELBUB)

#if U_STRUCTURE == SCALAR

#define N_DATA          28672
#define F_DATA          114688
#define E_DATA          1032192
#define DATA_S          575
#define T_FOR_BC        2
#define N_OF_NODE_FUNC  1
#if (U_SPACE == GP1C) || (U_SPACE == GP1X3C) || (U_SPACE == GP1C_ELBUB) || (U_SPACE == GQ1C) || (U_SPACE == GQ1C_ELBUB)
#define N_OF_FACE_FUNC  0
#elif (U_SPACE == GP2C) || (U_SPACE == GP2X3C) || (U_SPACE == GP2C_3ELBUB) || (U_SPACE == GP2C_6ELBUB) || (U_SPACE == GQ1X4C) || (U_SPACE == GQ2C) || (U_SPACE == GQ2C_2ELBUB)  || (U_SPACE == GQ2C_3ELBUB)
#define N_OF_FACE_FUNC  1
#elif U_SPACE == GQ2B3C
#define N_OF_FACE_FUNC  2
#define DIRECTION_OF_EDGES   NO
#elif U_SPACE == GQ2X4C
#define N_OF_FACE_FUNC  3
#endif
#if (U_SPACE == GP1C) || (U_SPACE == GQ1C) || (U_SPACE == GP2C)
#define N_OF_ELEM_FUNC  0
#elif (U_SPACE == GP1X3C) || (U_SPACE == GP1C_ELBUB) || (U_SPACE == GQ1X4C) || (U_SPACE == GQ1C_ELBUB) || (U_SPACE == GQ2C)
#define N_OF_ELEM_FUNC  1
#elif U_SPACE == GQ2B3C
#define N_OF_ELEM_FUNC  4
#elif (U_SPACE == GP2C_3ELBUB) || (U_SPACE == GQ2C_2ELBUB)
#define N_OF_ELEM_FUNC  3
#elif (U_SPACE == GP2X3C) || (U_SPACE == GQ2C_3ELBUB)
#define N_OF_ELEM_FUNC  4
#elif U_SPACE == GP2C_6ELBUB
#define N_OF_ELEM_FUNC  6
#elif U_SPACE == GQ2X4C
#define N_OF_ELEM_FUNC  9
#endif

#elif U_STRUCTURE == VECTOR

#endif

#elif U_SPACE == P2C

#define DATA_S          15
#define N_DATA          73 
#define F_DATA          273
#define E_DATA          0
#define T_FOR_BC        2

#elif U_SPACE == Q1ROT

#define DATA_S          15
#define E_DATA          0
#define T_FOR_BC        2

#if DATA_STR & KORN_MATRIX
#define N_DATA          446
#define F_DATA         1576
#else
#define N_DATA          379
#define F_DATA         1329
#endif

#endif

#endif  /*  P_SPACE == NONE  */
