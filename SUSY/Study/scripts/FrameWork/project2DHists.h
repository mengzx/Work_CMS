#ifndef project2DHists_h
#define project2DHists_h

#include <vector>
#include <iostream>
#include <fstream>
#include "TH2D.h"
#include "TString.h"
#include <stdio.h>

class project2DHists{
 public:
  project2DHists();
  ~project2DHists(){;}

  TH1D* projectX( TH2D* inh, double reqlowe, double requpe );
  TH1D* projectY( TH2D* inh, double reqlowe, double requpe );

};

#endif


