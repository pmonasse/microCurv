cmake_minimum_required(VERSION 2.6)

# zlib is needed by libPNG and by libTIFF
find_package(ZLIB)
if(NOT ZLIB_FOUND)
  add_subdirectory(third_party/zlib-1.2.11)
endif(NOT ZLIB_FOUND)
include_directories(${ZLIB_INCLUDE_DIRS})

find_package(PNG)
find_package(TIFF)

if(NOT PNG_FOUND)
    add_subdirectory(third_party/libpng-1.6.29)
    include_directories(${PNG_INCLUDE_DIRS})
endif(NOT PNG_FOUND)
add_definitions(${PNG_DEFINITIONS})

if(NOT TIFF_FOUND)
    find_package(JPEG)
    if(NOT JPEG_FOUND)
        add_subdirectory(third_party/jpeg-9b)
    endif(NOT JPEG_FOUND)
    include_directories(${JPEG_INCLUDE_DIR})
    add_subdirectory(third_party/tiff-4.0.7)
    include_directories(${TIFF_INCLUDE_DIR})
endif(NOT TIFF_FOUND)

set(SRC
    cmdLine.h
    curv.cpp curv.h
    draw_curve.cpp draw_curve.h
    fill_curve.cpp fill_curve.h
    gass.cpp gass.h
    image.cpp image.h
    io_png.c io_png.h
    io_tiff.c io_tiff.h
    levelLine.cpp levelLine.h
    lltree.cpp lltree.h
    main.cpp
    microCurv.h microCurv.cpp
    xmtime.h)
add_executable(microCurv ${SRC})
target_link_libraries(microCurv ${PNG_LIBRARIES} ${TIFF_LIBRARIES})

find_package(OpenMP)
if(OPENMP_FOUND)
    set_target_properties(microCurv PROPERTIES COMPILE_FLAGS
                          ${OpenMP_CXX_FLAGS})
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set(CMAKE_EXE_LINKER_FLAGS ${OpenMP_CXX_FLAGS})
    endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
endif(OPENMP_FOUND)

if(UNIX)
    set_target_properties(microCurv PROPERTIES COMPILE_FLAGS
                          "-Wall -Wextra -Werror ${OpenMP_CXX_FLAGS}")
endif(UNIX)
