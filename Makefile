include ./makesrc
include ./makedeps

CXX = g++-4.2
COMPILE_ONLY_OPTIONS = $(INCLUDE) -O3
ALL_OPTIONS = $(INCLUDE) $(LIBS) -O3 

contourmatcher: contourmatcher_opencv \
				contourmatcher_opencv_tests

contourmatcher_opencv_tests: contourmatcher_opencv $(CM_OPENCV_TEST_DEPS)
	$(CXX) $(CM_OPENCV_TEST_SRC) $(ALL_OPTIONS)

contourmatcher_opencv: $(CM_OPENCV_DEPS)
	$(CXX) $(CM_OPENCV_SRC) -c $(COMPILE_ONLY_OPTIONS)