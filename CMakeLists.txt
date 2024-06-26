cmake_minimum_required (VERSION 3.9)

project(fsm)

find_package(GDAL)

option(USE_GDAL "GDAL is required to save FSM outputs" ON)

find_package(Boost COMPONENTS serialization system REQUIRED)
message(STATUS "Found Boost ${BOOST_VERSION}")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache) # Less useful to do it for linking, see edit2
endif(CCACHE_FOUND)

set(RICHDEM_NO_PROGRESS ON CACHE BOOL "Whether or not to show RichDEM progress bars")
add_subdirectory(submodules/dephier EXCLUDE_FROM_ALL)

add_library(fsm INTERFACE)
target_include_directories(fsm INTERFACE include/)

# include_directories(submodules)
option(USE_BOOST "Compile code using Boost to enable input-output of the full depression hierarchy." OFF)
if(USE_BOOST)
  target_link_libraries(fsm INTERFACE dephier Boost::serialization)

  add_executable(fsm_boost.exe src/main_boost.cpp)
  target_link_libraries(fsm_boost.exe PRIVATE fsm richdem)
  target_compile_features(fsm_boost.exe PUBLIC cxx_auto_type cxx_std_17)

else()
  target_link_libraries(fsm INTERFACE dephier)

  add_executable(fsm.exe src/main.cpp)
  target_link_libraries(fsm.exe PRIVATE fsm richdem)
  target_compile_features(fsm.exe PUBLIC cxx_auto_type cxx_std_17)

  message(WARNING "You are compiling without Boost. You will not be able to save or load the depression hierarchy if running FSM multiple times in the same region.")
endif()

add_executable(fsm_paper_tests.exe paper/fsm_paper_tests.cpp)
target_link_libraries(fsm_paper_tests.exe PRIVATE fsm richdem)
target_compile_features(fsm_paper_tests.exe PUBLIC cxx_auto_type cxx_std_17)

add_executable(fsm_paper_speed_comparison.exe paper/fsm_paper_speed_comparison.cpp)
target_link_libraries(fsm_paper_speed_comparison.exe PRIVATE fsm richdem)
target_compile_features(fsm_paper_speed_comparison.exe PUBLIC cxx_auto_type cxx_std_17)


add_executable(fsm_paper_scaling_tests.exe paper/fsm_paper_scaling_tests.cpp)
target_link_libraries(fsm_paper_scaling_tests.exe PRIVATE fsm richdem)
target_compile_features(fsm_paper_scaling_tests.exe PUBLIC cxx_auto_type cxx_std_17)

add_executable(fsm_unittests.exe
  unittests/test_main.cpp
  unittests/fsm_tests.cpp
)
target_link_libraries(fsm_unittests.exe PRIVATE fsm richdem)
target_compile_features(fsm_unittests.exe PUBLIC cxx_auto_type cxx_std_17)

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(fsm_unittests.exe PRIVATE -O0 -g --coverage)
  target_compile_definitions(fsm_unittests.exe PRIVATE -DCODE_COVERAGE)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(fsm_unittests.exe PRIVATE --coverage)
  else()
    target_link_libraries(fsm_unittests.exe PRIVATE --coverage)
  endif()

  add_custom_command(TARGET fsm_unittests.exe POST_BUILD
    COMMAND ./fsm_unittests.exe
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/gcoverage
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/lcoverage
    COMMAND gcovr -r ${CMAKE_CURRENT_SOURCE_DIR} --html --html-details -o ${CMAKE_BINARY_DIR}/gcoverage/coverage.html
    COMMAND lcov --capture --directory ${CMAKE_CURRENT_SOURCE_DIR} --output-file ${CMAKE_BINARY_DIR}/lcoverage/coverage.info
    COMMAND genhtml ${CMAKE_BINARY_DIR}/lcoverage/coverage.info --output-directory ${CMAKE_BINARY_DIR}/lcoverage/
  )
endif()
