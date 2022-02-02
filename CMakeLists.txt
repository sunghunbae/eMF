CMAKE_MINIMUM_REQUIRED (VERSION 3.2)

# we will use FindGSL which is new in cmake version 3.2
# https://cmake.org/cmake/help/latest/module/FindGSL.html

PROJECT(eMF)

# chang the installation location (default: ~/bin)
SET(CMAKE_INSTALL_PREFIX ~/bin )

SET(USE_STATIC_LIBS off )
SET(CMAKE_VERBOSE_MAKEFILE off)
SET(CMAKE_CXX_FLAGS "-O3 -Wparentheses")

IF   ( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  LINK_DIRECTORIES ( /usr/lib /lib )
ELSE ( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  LINK_DIRECTORIES ( /usr/lib64 /lib64 )
ENDIF( CMAKE_SIZEOF_VOID_P EQUAL 8 )

FIND_PACKAGE (GSL REQUIRED)

INCLUDE_DIRECTORIES( ${GSL_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( ${eMF_SOURCE_DIR}/src )

SET ( emf
	src/emf.h
	src/func.cxx
	src/chisq.cxx
	src/aniso.cxx
	src/grid.cxx
	src/mle.cxx
	src/gmin.cxx
	src/brent.cxx
	src/conjgr.cxx
	src/simplex.cxx
	src/powell.cxx
	src/levmar.cxx
	src/anneal.cxx
	src/fitmodel.cxx
	src/rjmap.cxx
	src/file.cxx
	src/data.cxx
	src/main.cxx
	src/print.cxx
	)

ADD_EXECUTABLE ( eMF  ${emf} )
TARGET_LINK_LIBRARIES ( eMF ${GSL_LIBRARIES} m )
SET_TARGET_PROPERTIES ( eMF PROPERTIES OUTPUT_NAME eMF )
INSTALL (TARGETS eMF DESTINATION ${CMAKE_INSTALL_PREFIX})
INSTALL (PROGRAMS ${eMF_SOURCE_DIR}/src/eMF-digest DESTINATION ${CMAKE_INSTALL_PREFIX})
