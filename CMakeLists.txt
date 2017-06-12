cmake_minimum_required(VERSION 2.8)

find_package(PNG)
find_package(TIFF)

if(NOT PNG_FOUND OR NOT TIFF_FOUND)
    find_package(ZLIB) # zlib is needed by libPNG and by libTIFF
    if(NOT ZLIB_FOUND)
        add_subdirectory(third_party/zlib-1.2.11)
        set(ZLIB_FOUND TRUE)
    endif()
    include_directories(${ZLIB_INCLUDE_DIRS})

    if(NOT PNG_FOUND)
        set(SKIP_INSTALL_ALL 1)
        add_subdirectory(third_party/libpng-1.6.29)
    endif()

    if(NOT TIFF_FOUND)
        find_package(JPEG)
        if(NOT JPEG_FOUND)
            add_subdirectory(third_party/jpeg-9b)
            set(JPEG_FOUND TRUE)
        endif()
        include_directories(${JPEG_INCLUDE_DIR})
        add_subdirectory(third_party/tiff-4.0.7)
    endif()
endif()

add_definitions(${PNG_DEFINITIONS})
include_directories(${PNG_INCLUDE_DIRS})
include_directories(${TIFF_INCLUDE_DIR})

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
endif()

if(UNIX)
    set_target_properties(microCurv PROPERTIES COMPILE_FLAGS
                          "-Wall -Wextra -Werror ${OpenMP_CXX_FLAGS}")
endif()
