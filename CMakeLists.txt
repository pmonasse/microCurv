CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CPACK_PACKAGE_FILE_NAME "microCurv")
SET(CPACK_GENERATOR "TGZ")
INCLUDE(CPack)

FIND_PACKAGE(PNG REQUIRED)
FIND_PACKAGE(TIFF REQUIRED)
ADD_DEFINITIONS(${PNG_DEFINITIONS})

SET(SRC1
    cmdLine.h
    fill_curve.cpp fill_curve.h
    gass.cpp gass.h
    levelLine.cpp levelLine.h
    lltree.cpp lltree.h
    io_png.c io_png.h
    main.cpp)
ADD_EXECUTABLE(extractLines ${SRC1})
TARGET_LINK_LIBRARIES(extractLines ${PNG_LIBRARIES})

SET(SRC2
   gass.cpp gass.h
   levelLine.cpp levelLine.h
   main_gass.cpp)
ADD_EXECUTABLE(smoothLines ${SRC2})

SET(SRC3
    cmdLine.h
    curv.cpp curv.h
    draw_curve.cpp draw_curve.h
    fill_curve.cpp fill_curve.h
    gass.cpp gass.h
    io_png.c io_png.h
    io_tiff.c io_tiff.h
    levelLine.cpp levelLine.h
    lltree.cpp lltree.h
    mirror.cpp mirror.h
    microCurv.cpp)
ADD_EXECUTABLE(microCurv ${SRC3})
TARGET_LINK_LIBRARIES(microCurv ${PNG_LIBRARIES} ${TIFF_LIBRARIES})

LIST(APPEND SRC1 ${SRC2} ${SRC3})
LIST(REMOVE_DUPLICATES SRC1)
INSTALL(FILES ${SRC1} CMakeLists.txt DESTINATION .)
