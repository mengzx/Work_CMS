#include "menus_base.h"
#include "TString.h"
#include "iostream.h"
#include "TLegend.h"
//#include "vectors.h"
#include <stdio.h>
#include <stdlib.h>

menus_base::menus_base(){

  debug_=-1;
  drawOverflow_=false;
  ratioPlotErr_=0.1;
  mcscale_=1.;
}

