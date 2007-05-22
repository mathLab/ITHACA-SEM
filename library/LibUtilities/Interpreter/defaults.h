/*****************************************************************************
 * DEFAULTS.H:  default parameter initializations for initial.y.
 * All variables are now in a single table, everything is double precision.
 *
 * $Id: defaults.h,v 6.9 2005/06/14 04:15:09 hmb Exp $
 *****************************************************************************/

static struct {
  char*  name;
  double cval;
} consts[] = {

  /* -- Mathematical constants. */

  "E"           ,   2.71828182845904523536 ,
  "DEG"         ,  57.29577951308232087721 ,
  "PI"          ,   3.14159265358979323844 ,
  "TWOPI"       ,   6.28318530717958647688 ,
  "EULER"       ,   0.57721566490153286061 ,
  "GOLDEN"      ,   1.61803398874989484820 ,  
  
  0             ,   0.0
};




