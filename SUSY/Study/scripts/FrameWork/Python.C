#include "menus_base.h"
#include "playHist1D.h"
#include "playHist2D.h"
#include "project2DHists.h"
#include "tdrstyle.h"

#include <Python.h>
#include <boost/python.hpp> 
#include <boost/python/suite/indexing/vector_indexing_suite.hpp> 

using namespace boost::python;
BOOST_PYTHON_MODULE(libFrameWork) {

  class_< menus_base >( "menus_base", boost::python::no_init ).;
  class_< playHist1D, bases<menus_base> >( "playHist1D", boost::python::no_init )
    .def( "getHist1D", &playHist1D::getHist1D )
    .def( "cloneHist1D", &playHist1D::cloneHist1D )
    .def( "getifileidir1D", &playHist1D::getifileidir1D )
    .def( "getHistInvFandvDir1D", &playHist1D::getHistInvFandvDir1D )
    .def( "addHistForDiffFoldersAndFiles1D", &playHist1D::addHistForDiffFoldersAndFiles1D )
    .def( "getifileidirih1D", &playHist1D::getifileidirih1D )
    .def( "getHistInvFvDirvH1D", &playHist1D::getHistInvFvDirvH1D )
    .def( "addHistForDiffFoldersFilesHists1D", &playHist1D::addHistForDiffFoldersFilesHists1D )
    .def( "formatHist", &playHist1D::formatHist )
    .def( "getOverflowbin", &playHist1D::getOverflowbin )
    .def( "getOverflowbinErr", &playHist1D::getOverflowbinErr )
    .def( "MaxHist", &playHist1D::MaxHist )
    .def( "SortHists", &playHist1D::SortHists )
    .def( "SortHists_index", &playHist1D::SortHists_index )
    .def( "invSortHists", &playHist1D::invSortHists )
    .def( "invSortHists_index", &playHist1D::invSortHists_index )
    .def( "MaxHist_index", &playHist1D::MaxHist_index )
    .def( "CumulativeH", &playHist1D::CumulativeH )
    //    .def( "CumulativeH", &playHist1D::CumulativeH )
    .def( "getRatioPlot", &playHist1D::getRatioPlot )
    .def( "getRatioErr", &playHist1D::getRatioErr )
    .def( "getFirstBinHasContent", &playHist1D::getFirstBinHasContent )
    .def( "getLastBinHasContent", &playHist1D::getLastBinHasContent )
    ;

  class_< playHist2D, bases<menus_base> >( "playHist2D", boost::python::no_init )
    .def( "getHist2D", &playHist2D::getHist2D )
    .def( "cloneHist2D", &playHist2D::cloneHist2D )
    .def( "getifileidir2D", &playHist2D::getifileidir2D )
    .def( "getHistInvFandvDir2D", &playHist2D::getHistInvFandvDir2D )
    .def( "addHistForDiffFoldersAndFiles2D", &playHist2D::addHistForDiffFoldersAndFiles2D )
    .def( "getifileidirih2D", &playHist2D::getifileidirih2D )
    .def( "getHistInvFvDirvH2D", &playHist2D::getHistInvFvDirvH2D )
    .def( "addHistForDiffFoldersFilesHists2D", &playHist2D::addHistForDiffFoldersFilesHists2D )
    .def( "addHistForDiffFoldersAndFiles_SubtrackHists2D", &playHist2D::addHistForDiffFoldersAndFiles_SubtrackHists2D )
    .def( "ReFillHist_AlphaTVSHT", &playHist2D::ReFillHist_AlphaTVSHT )
    .def( "ReFillHist_low", &playHist2D::ReFillHist_low )
    .def( "ReFillHist_high", &playHist2D::ReFillHist_high )
    .def( "formatHist", &playHist2D::formatHist )
    .def( "Lines", &playHist2D::Lines )
    .def( "CumulativeH", &playHist2D::CumulativeH )
    ;
  class_< project2DHists >( "project2DHists", boost::python::no_init )
    .def( "projectX", &project2DHists::projectX )
    .def( "projectY", &project2DHists::projectY )
    ;

  class_< tdrstyle >( "tdrstyle", boost::python::no_init )
    .def( "setTDRStyle", &tdrstyle::setTDRStyle )

}
