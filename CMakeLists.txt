project(eus_math)

cmake_minimum_required(VERSION 2.4.6)

find_package(Eigen REQUIRED)
find_package(OpenCV REQUIRED)

include_directories(${Eigen_INCLUDE_DIRS})
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${OpenCV_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/src)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_library(eus_eigen SHARED src/eus_eigen.cpp)

add_executable(eus_cv src/eus_cv.cpp)
target_link_libraries(eus_cv ${OpenCV_LIBRARIES})
