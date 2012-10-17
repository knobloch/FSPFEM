#if SC_EXAMPLE == 1                /*  Mizukami, Hughes, CMAME 50 (1985)  */
                                   /*  Example 4.1                        */
#define VALUE_OF_EPS  1.e-7
#define PHI           (PI/4.)
#define THETA         (PI/4.)
#define NV_POINTS     11           
#define CUBE_TYPE     2
#include "bc-cd/bc899.h"

#elif  SC_EXAMPLE == 2             /*  Mizukami, Hughes, CMAME 50 (1985)  */
                                   /*  Example 4.1                        */
#define VALUE_OF_EPS  1.e-7
#define THETA         (PI/3.)
#define NV_POINTS     11          
#define CUBE_TYPE     2
#include "bc-cd/bc899.h"

#elif  SC_EXAMPLE == 22

#define VALUE_OF_EPS  1.e-7
#define THETA         (PI/3.)
#define NV_POINTS     21          
#define CUBE_TYPE     3
#include "bc-cd/bc899a.h"

#elif  SC_EXAMPLE == 3             /*  Mizukami, Hughes, CMAME 50 (1985)  */
                                   /*  Example 4.2                        */
#define VALUE_OF_EPS  1.e-7
#define THETA         (PI/4.)
#define NV_POINTS     21         
#define CUBE_TYPE     2
#include "bc-cd/bc898.h"

#elif  SC_EXAMPLE == 4             /*  Mizukami, Hughes, CMAME 50 (1985)  */
                                   /*  Example 4.2                        */
#define VALUE_OF_EPS  1.e-7
#define THETA         (PI/3.)
#define NV_POINTS     21        
#define CUBE_TYPE     2
#include "bc-cd/bc898.h"

#elif  SC_EXAMPLE == 42            /*  Mizukami, Hughes, CMAME 50 (1985)  */
                                   /*  Example 4.2                        */
#define VALUE_OF_EPS  1.e-7
#define THETA         (PI/3.)
#define NV_POINTS     11        
#define CUBE_TYPE     2
#include "bc-cd/bc898.h"

#elif  SC_EXAMPLE == 5             /*  Mizukami, Hughes, CMAME 50 (1985)  */
                                   /*  Example 4.3                        */
#define VALUE_OF_EPS  1.e-7
#define NV_POINTS     11       
#define CUBE_TYPE     2
#include "bc-cd/bc897.h"

#elif  SC_EXAMPLE == 55            /*  Mizukami, Hughes, CMAME 50 (1985)  */
                                   /*  Example 4.3                        */
#define VALUE_OF_EPS  1.e-8
#define NV_POINTS     33        
#define CUBE_TYPE     2
#include "bc-cd/bc897.h"
//#include "bc-cd/bc897r.h"

#elif  SC_EXAMPLE == 52            /*  Mizukami, Hughes, CMAME 50 (1985)  */
                                   /*  Example 4.3                        */
#define VALUE_OF_EPS  1.e-7
#define NV_POINTS     6        
#define CUBE_TYPE     2
#include "bc-cd/bc897.h"

#elif  SC_EXAMPLE == 6             /*  Mizukami, Hughes, CMAME 50 (1985)  */
                                   /*  Example 4.4                        */
#define VALUE_OF_EPS  1.e-7
#define NV_POINTS     11      
#define CUBE_TYPE     2
#include "bc-cd/bc896.h"

#elif  SC_EXAMPLE == 66
#define VALUE_OF_EPS  1.e-6
#define NV_POINTS     65
#define CUBE_TYPE     1
#include "bc-cd/bc840.h"

#elif  SC_EXAMPLE == 67
#define VALUE_OF_EPS  1.e-8
#define NV_POINTS     65
#define CUBE_TYPE     2
#include "bc-cd/bc841.h"

#elif  SC_EXAMPLE == 68
#define VALUE_OF_EPS  1.e-8
#define NV_POINTS     33
#define CUBE_TYPE     2
#include "bc-cd/bc842.h"

#elif  SC_EXAMPLE == 7             /*                                     */
                                   /*  donut                              */
#define VALUE_OF_EPS  1.e-6
#define NV_POINTS     33      
#define CUBE_TYPE     1
#include "bc-cd/bc999.h"

#elif  SC_EXAMPLE == 72            /*                                     */
                                   /*  donut                              */
#define VALUE_OF_EPS  1.e-6
#define NV_POINTS     17      
#define CUBE_TYPE     1
#include "bc-cd/bc999.h"

#elif  SC_EXAMPLE == 8         /*  Hughes, Mallet, Mizukami, CMAME 54 (1986)  */
                               /*  Example 4.1                                */
#define VALUE_OF_EPS  1.e-8
#define THETA         (PI/3.)
#define NV_POINTS     33        
#define CUBE_TYPE     2
#include "bc-cd/bc898a.h"

#elif  SC_EXAMPLE == 80
#define VALUE_OF_EPS  1.e-8
#define THETA         (PI/3.)
#define NV_POINTS     33        
#define CUBE_TYPE     2
#include "bc-cd/bc898n.h"

#elif  SC_EXAMPLE == 81
#define VALUE_OF_EPS  1.e-8
#define THETA         (PI/2.1)
#define NV_POINTS     33        
#define CUBE_TYPE     1
#include "bc-cd/bc830.h"

#elif  SC_EXAMPLE == 82        /*  Hughes, Mallet, Mizukami, CMAME 54 (1986)  */
                               /*  Example 4.1                                */
#define VALUE_OF_EPS  1.e-7
#define THETA         (PI/3.)
#define NV_POINTS     11        
#define CUBE_TYPE     2
#include "bc-cd/bc898a.h"

#elif  SC_EXAMPLE == 88        /*  Hughes, Mallet, Mizukami, CMAME 54 (1986)  */
                               /*  Example 4.1                                */
#define VALUE_OF_EPS  1.e-8
#define THETA         (PI/3.)
#define NV_POINTS     65        
#define CUBE_TYPE     2
#include "bc-cd/bc898a.h"

#elif  SC_EXAMPLE == 9             /*  Knopp, Lube, Rapin, CMAME 191 (2002)  */
                                   /*  Example 6.2                           */
#define VALUE_OF_EPS  1.e-8
#define NV_POINTS     33        
#define CUBE_TYPE     1
#include "bc-cd/bc888.h"

#elif  SC_EXAMPLE == 91
#define VALUE_OF_EPS  1.e-8
#define NV_POINTS     33        
#define CUBE_TYPE     1
#include "bc-cd/bc888a.h"

#elif  SC_EXAMPLE == 92            /*  Knopp, Lube, Rapin, CMAME 191 (2002)  */
                                   /*  Example 6.2                           */
#define VALUE_OF_EPS  1.e-4
#define NV_POINTS     33        
#define CUBE_TYPE     1
#include "bc-cd/bc888.h"

#elif  SC_EXAMPLE == 93            /*  Knopp, Lube, Rapin, CMAME 191 (2002)  */
                                   /*  Example 6.2                           */
#define VALUE_OF_EPS  TNU  
#define NV_POINTS     65        
#define CUBE_TYPE     2
#include "bc-cd/bc888.h"

#elif  SC_EXAMPLE == 10            /*  Knopp, Lube, Rapin, CMAME 191 (2002)  */
                                   /*  Example 6.2                           */
#define VALUE_OF_EPS  1.e-8
#define NV_POINTS     65        
#define CUBE_TYPE     1
#include "bc-cd/bc888.h"

#elif  SC_EXAMPLE == 19            /*  John, Maubach, Tobiska, NM 78 (1985)  */
                                   /*  Example 2                             */
#define VALUE_OF_EPS  1.e-8
#define NV_POINTS     65        
#define CUBE_TYPE     1
#include "bc-cd/bc892.h"

#elif  SC_EXAMPLE == 20            /*  non-constant rhs                      */
                                   /*                                        */
#define VALUE_OF_EPS  1.e-8
#define NV_POINTS     33        
#define CUBE_TYPE     2
#include "bc-cd/bc894.h"

#elif  SC_EXAMPLE == 25            /*  Hemker                                */
                                   /*                                        */
#define VALUE_OF_EPS  1.e-6
#define NV_POINTS     65
#define CUBE_TYPE     1
#include "bc-cd/bc893.h"

#elif  SC_EXAMPLE == 70            /*  Burman, Hansbo, CMAME (2004)       */
                                   /*  Example 3.2                        */
#define VALUE_OF_EPS  1.e-5
#define NV_POINTS     21           
#define CUBE_TYPE     2
#include "bc-cd/bc895.h"

#elif  SC_EXAMPLE == 700
                       
#define VALUE_OF_EPS  1.e-8
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     33
#define CUBE_TYPE     1
#include "bc-cd/bc7.h"

#elif  SC_EXAMPLE == 701
                       
#define VALUE_OF_EPS  1.e-6
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     21
#define CUBE_TYPE     1
#include "bc-cd/bc7a.h"

#elif  SC_EXAMPLE == 702
                       
#define VALUE_OF_EPS  1.e-6
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     41
#define CUBE_TYPE     1
#include "bc-cd/bc7a.h"

#elif  SC_EXAMPLE == 703
                       
#define VALUE_OF_EPS  1.e-6
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     81
#define CUBE_TYPE     1
#include "bc-cd/bc7a.h"

#elif  SC_EXAMPLE == 704
                       
#define VALUE_OF_EPS  1.e-6
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     161
#define CUBE_TYPE     1
#include "bc-cd/bc7a.h"

#elif  SC_EXAMPLE == 800
                       
#define VALUE_OF_EPS  1.e-8
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     33
#define CUBE_TYPE     1
#include "bc-cd/bc8.h"

#elif  SC_EXAMPLE == 801
                       
#define VALUE_OF_EPS  1.e-8
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     21
#define CUBE_TYPE     1
#include "bc-cd/bc8a.h"

#elif  SC_EXAMPLE == 802
                       
#define VALUE_OF_EPS  1.e-8
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     41
#define CUBE_TYPE     1
#include "bc-cd/bc8a.h"

#elif  SC_EXAMPLE == 803
                       
#define VALUE_OF_EPS  1.e-8
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     81
#define CUBE_TYPE     1
#include "bc-cd/bc8a.h"

#elif  SC_EXAMPLE == 804
                       
#define VALUE_OF_EPS  1.e-8
double  SNU=VALUE_OF_EPS;
#define NV_POINTS     161
#define CUBE_TYPE     1
#include "bc-cd/bc8a.h"

#elif  SC_EXAMPLE == 805   /* Matthies, Skrzypacz, Tobiska, Preprint 44, 2007 */
                       
#define VALUE_OF_EPS  1.e-7
#define NV_POINTS     33
#define CUBE_TYPE     1
#include "bc-cd/bc801a.h"

#elif  SC_EXAMPLE == 810         /*  reaction-diffusion eq. with f = 1  */
                       
#define VALUE_OF_EPS  1.e-7
#define NV_POINTS     11 
#define CUBE_TYPE     2
#include "bc-cd/bc810.h"

#elif  SC_EXAMPLE == 811         /*  Burman, Ern, CMAME 191 (2002), 3833-3855 */
                       
#define VALUE_OF_EPS  1.e-7
#define NV_POINTS     21 
#define CUBE_TYPE     2
#include "bc-cd/bc811.h"

#elif  SC_EXAMPLE == 820
                       
#define VALUE_OF_EPS  1.e-7
#define NV_POINTS     21 
#define CUBE_TYPE     2
#include "bc-cd/bc820.h"

#elif  SC_EXAMPLE == 888          /*  rotating profile  */

#define VALUE_OF_EPS  1.e-6
#define NV_POINTS     65          
#include "bc-cd/bc888.h"

#elif  SC_EXAMPLE == 0

#else

#define VALUE_OF_EPS  1.e-7
#define NV_POINTS     11           /*  value of NV  */
#include "bc-cd/bc896.h"

#endif
