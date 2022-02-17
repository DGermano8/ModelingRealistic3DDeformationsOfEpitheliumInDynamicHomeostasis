/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

#include "CommandLineArguments.hpp"
int main( int argc, char *argv[] ) {
 CommandLineArguments::Instance()->p_argc = &argc;
 CommandLineArguments::Instance()->p_argv = &argv;
 return CxxTest::ErrorPrinter().run();
}
#include "projects/DeformableTissue_Cleaned/test/TestsWithAreaDep/TestPeriodicBendingForce_WntStrip.hpp"

static Test3dBoxModel suite_Test3dBoxModel;

static CxxTest::List Tests_Test3dBoxModel = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_Test3dBoxModel( "projects/DeformableTissue_Cleaned/test/TestsWithAreaDep/TestPeriodicBendingForce_WntStrip.hpp", 68, "Test3dBoxModel", suite_Test3dBoxModel, Tests_Test3dBoxModel );

static class TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts() : CxxTest::RealTestDescription( Tests_Test3dBoxModel, suiteDescription_Test3dBoxModel, 77, "TestPeriodicCubeWithGhosts" ) {}
 void runTest() { suite_Test3dBoxModel.TestPeriodicCubeWithGhosts(); }
} testDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts;

#include <cxxtest/Root.cpp>
