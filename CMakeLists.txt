set(VTK_DIR /usr/local/)

cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

#-----------------------------------------------------
project(LambdaPoly)

find_package(VTK CONFIG)
find_package(VTK REQUIRED)

if (NOT VTK_FOUND)
  message("Skipping LambdaPoly: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()

message (STATUS "VTK_VERSION: ${VTK_VERSION}")
message (STATUS "VTK_LIBRARIES: ${VTK_LIBRARIES}")

set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_C_FLAGS "-O0 -ggdb")
set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb")
set(CMAKE_C_FLAGS_RELEASE "-O0 -ggdb")
set(CMAKE_CXX_FLAGS "-O0 -ggdb")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "-O0 -ggdb")

include_directories(
		"/usr/local/include/vtk-9.0/"
)

link_directories(
		"/usr/local/lib/"
		"/usr/lib/x86_64-linux-gnu/"
)

# set(LD_LIBRARY_PATH "/usr/local/lib/")

add_executable(LambdaPoly LambdaPoly.cxx )
target_link_libraries(LambdaPoly ${VTK_LIBRARIES})
# target_link_libraries(LambdaPoly mpfr gmp)

# Turn this on for memory leak debugging
set(VTK_DEBUG_LEAKS ON)

vtk_module_autoinit(
    TARGETS LambdaPoly
    MODULES ${VTK_LIBRARIES} 
) 

