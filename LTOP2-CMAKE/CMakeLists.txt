cmake_minimum_required(VERSION 3.0)

# Define the name of your project
project(ltop2)

# Build some boost libs using an external project
include(boost.cmake)

# Add fftw3 library
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(LTOP2-CPP/fftw-3.3.8)

# Define your target
add_executable(${PROJECT_NAME} LTOP2-CPP/main.cpp LTOP2-CPP/gridmap.cpp LTOP2-CPP/helper.cpp LTOP2-CPP/interp.cpp LTOP2-CPP/farfield.cpp)

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)

add_custom_command(TARGET ltop2
        POST_BUILD
        COMMAND find . -name "*.sh" | xargs chmod a+x
        )

# Add necessary libs
target_link_libraries(${PROJECT_NAME} fftw3 boost::filesystem boost::iostreams)
target_include_directories(${PROJECT_NAME} PUBLIC LTOP2-CPP/fftw-3.3.8/api)
