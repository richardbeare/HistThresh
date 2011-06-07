#Change PROJECT_NAME to the name of your project
PROJECT(HistThresh)

# required to avoid warnings while running cmake
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.2)

# cpack stuff
INCLUDE(CPack)

# set the name of the input image used to run the test
SET(INPUT_IMAGE ${CMAKE_CURRENT_SOURCE_DIR}/images/cthead1.png)

#include some macros from another file...
INCLUDE(IJMacros.txt)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)
FIND_PROGRAM(ITK_TEST_DRIVER itkTestDriver)
SET(TEST_COMMAND ${ITK_TEST_DRIVER} --add-before-env PATH ${CMAKE_CURRENT_BINARY_DIR})

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
SET(Required_Packages
  ITK
)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
 ${CMAKE_CURRENT_SOURCE_DIR}/Modules/maskAdaptor/
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
  ITKCommon
  ITKBasicFilters
  ITKStatistics
  ITKIO
)


# install devel files
OPTION(INSTALL_DEVEL_FILES "Install C++ headers" ON)
IF(INSTALL_DEVEL_FILES)
FILE(GLOB develFiles *.h *.txx) 
FOREACH(f ${develFiles})
  INSTALL_FILES(/include/InsightToolkit/BasicFilters FILES ${f})
ENDFOREACH(f)
ENDIF(INSTALL_DEVEL_FILES)



# option for wrapping
OPTION(BUILD_WRAPPERS "Wrap library" OFF)
IF(BUILD_WRAPPERS)
   SUBDIRS(Wrapping)
ENDIF(BUILD_WRAPPERS)
   
   

#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above. 
#You can build more than one executable per project

IF(BUILD_TESTING)

FOREACH(CurrentExe "testTriangle" "testIntermodes" "testKittlerIllingworth" "testHuang" "testIsoData" "testLi" "testMaxEntropy" "testMoments" "testRenyiEntropy" "testShanbhag" "testYen")
  ADD_EXECUTABLE(${CurrentExe} ${CurrentExe}.cxx)
  TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})
ENDFOREACH(CurrentExe)

ENDIF(BUILD_TESTING)

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.
# ADD_TEST(Testname ExecutableToRun arg1 arg2 arg3)

ADD_TEST(testTriangle ${TEST_COMMAND}
   testTriangle ${INPUT_IMAGE} outTriangle.png
)

ADD_TEST(testIntermodes ${TEST_COMMAND}
   testIntermodes ${INPUT_IMAGE} outIntermodes.png outMinimum.png
)
ADD_TEST(testKittlerIllingworth ${TEST_COMMAND}
   testKittlerIllingworth ${INPUT_IMAGE} outKittlerIllingworth.png
)
ADD_TEST(testHuang ${TEST_COMMAND}
   testHuang ${INPUT_IMAGE} outHuang.png
)
ADD_TEST(testLi ${TEST_COMMAND}
   testLi ${INPUT_IMAGE} outLi.png
)
ADD_TEST(testMaxEntropy ${TEST_COMMAND}
   testMaxEntropy ${INPUT_IMAGE} outMaxEntropy.png
)
ADD_TEST(testMoments ${TEST_COMMAND}
   testMoments ${INPUT_IMAGE} outMoments.png
)
ADD_TEST(testRenyiEntropy ${TEST_COMMAND}
   testRenyiEntropy ${INPUT_IMAGE} outRenyiEntropy.png
)
ADD_TEST(testShanbhag ${TEST_COMMAND}
   testShanbhag ${INPUT_IMAGE} outShanbhag.png
)

ADD_TEST(testYen ${TEST_COMMAND}
   testYen ${INPUT_IMAGE} outYen.png
)
