cmake_minimum_required(VERSION 2.8)

find_package(PNG)
if(NOT PNG_FOUND OR NOT TIFF_FOUND)
    find_package(ZLIB) # zlib is needed by libPNG
    if(NOT ZLIB_FOUND)
        add_subdirectory(third_party/zlib-1.2.11)
        set(ZLIB_FOUND TRUE)
    endif()
    include_directories(${ZLIB_INCLUDE_DIRS})

    if(NOT PNG_FOUND)
        set(SKIP_INSTALL_ALL 1)
        add_subdirectory(third_party/libpng-1.6.37)
    endif()
endif()
add_definitions(${PNG_DEFINITIONS})
include_directories(${PNG_INCLUDE_DIRS})

set(SRC
    cmdLine.h
    draw_curve.cpp draw_curve.h
    fill_curve.cpp fill_curve.h
    levelLine.cpp levelLine.h
    lltree.cpp lltree.h)
add_executable(bilines ${SRC}
               io_png.c io_png.h
               bilines.cpp)
target_link_libraries(bilines ${PNG_LIBRARIES} ${TIFF_LIBRARIES})

if(UNIX)
    set_target_properties(bilines PROPERTIES COMPILE_FLAGS "-Wall -Wextra")
endif()

add_executable(test_extract ${SRC} test_extract.cpp)
