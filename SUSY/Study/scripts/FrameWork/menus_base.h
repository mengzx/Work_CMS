#ifndef menus_base_h
#define menus_base_h

#include "TString.h"
#include "TLegend.h"
#include <stdio.h>
#include <stdlib.h>

class menus_base{

 public:
  menus_base();
  ~menus_base(){;}

  bool drawOverflow_;
  int debug_;
  double mcscale_;
  double ratioPlotErr_;

  //for printTables
  int lowHTEdge;
  int highHTEdge;
  int lowATEdge;
  int highATEdge;

}; //end of menus
#endif
