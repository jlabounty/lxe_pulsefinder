cmake_minimum_required(VERSION 3.13)  # CMake version check
project(simple_example)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 17)            # Enable c++14 standard

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES  GlobalFitAlgorithm.cc ModelFunction.cc  test.cc)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(ROOT REQUIRED)
# include(${ROOT_USE_FILE})
# find_package(Minuit2 REQUIRED)

include_directories(${PROJECT_SOURCE_DIR} 
    ${ROOT_INCLUDE_DIRS} )

# target_include_directories(test PRIVATE standalone)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(test ${SOURCE_FILES} )
# target_link_libraries(test standalone )
target_link_libraries(test ${ROOT_LIBRARIES} )
target_link_libraries(test Eigen3::Eigen)