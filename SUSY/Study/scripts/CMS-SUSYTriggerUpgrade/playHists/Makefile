ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
CXXFLAGS      =-I$(ROOTSYS)/include -O -Wall -fPIC

SOURCEDIR	=  ../../FrameWork/

HEADERDIR	= ./

MAKE	=	make

object	=	menus.o vectors.o basicPlots.o printTables.o

objectFramWork	=	tdrstyle.o playHist2D.o playHist1D.o project2DHists.o menus_base.o rootdict.o

all: $(object) $(objectFramWork) main.exe
.PHONY : all

$(object): %.o: %.C
	g++ -c $? -o $@ -I $(ROOTCFLAGS) -I . -I $(SOURCEDIR)

$(objectFramWork): %.o : $(SOURCEDIR)%.C
	g++ -c $? -o $@ -I $(ROOTCFLAGS) -I . -I $(SOURCEDIR)

main.exe:main.C
	g++ -O3 -Wall -I $(ROOTCFLAGS) $^ -o $@ $(object) $(objectFramWork) -I . -I $(SOURCEDIR)  -I $(HEADERDIR) $(ROOTLIBS)

echo:
	echo "ROOTCFLAGS=$(ROOTCFLAGS)"
	echo "ROOTLIBS=$(ROOTLIBS)"
	echo "ROOTSYS=$(ROOTSYS)"
	echo "PATH=$(PATH)"
	echo "LD_LIBRARY_PATH=$(LD_LIBRARY_PATH)"

.PHONY : clean
clean: 
	rm -rf main.exe *.o *~ 
deadcor.exe:



