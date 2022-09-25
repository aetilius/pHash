cmake_minimum_required (VERSION 2.8)

include(ExternalProject)
include(GNUInstallDirs)

project(pHash)
set(CMAKE_MACOSX_RPATH 1)
set(pHash_VERSION_MAJOR 1)
set(pHash_VERSION_MINOR 0)
set(pHash_VERSION_PATCH 0)

OPTION(PHASH_DYNAMIC   "Build pHash dynamic library"                     ON)
OPTION(PHASH_STATIC    "Build pHash static library"                     OFF)
OPTION(WITH_AUDIO_HASH "Audio hash support"                             OFF)
OPTION(WITH_VIDEO_HASH "Video hash support"                             OFF)
OPTION(PHASH_BINDINGS  "Compile foreign programming languages bindings" OFF)
OPTION(PHASH_EXAMPLES  "Compile examples"                               OFF)

set(HAVE_IMAGE_HASH 1)

if(WITH_AUDIO_HASH)
    set(HAVE_AUDIO_HASH 1)
endif()

if(WITH_VIDEO_HASH)
    set(HAVE_VIDEO_HASH 1)
endif()

find_path(CIMG_H_DIR NAMES CImg.h PATHS "${CMAKE_CURRENT_SOURCE_DIR}/third-party/CImg")
include_directories(${CIMG_H_DIR})

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/pHash.h.cmake ${CMAKE_CURRENT_SOURCE_DIR}/src/pHash.h)


if(NOT WIN32)
       set(MMAN_FILE "")
       set(DIRENT_FILE "")
    EXEC_PROGRAM(uname ARGS -m OUTPUT_VARIABLE BUILD_SYSTEM)
    EXEC_PROGRAM(uname ARGS -s OUTPUT_VARIABLE CMAKE_SYSTEM_NAME)
else()
    set(MMAN_FILE "win/mman.cpp")
    set(DIRENT_FILE "win/dirent.cpp")
    add_definitions("-D_EXPORTING")
    include_directories("win")
endif()

if(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

if(PHASH_BINDINGS)
  add_subdirectory(bindings)
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BUILD_TYPE}")

link_directories(/usr/local/lib)
include_directories(/usr/local/include)


set(LIBS_DEPS png z jpeg tiff)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
set(LIBS_DEPS ${LIBS_DEPS} Threads::Threads)

find_library(LIBMPG123 mpg123)
if(LIBMPG123)
    message("libmpg123 found at ${LIBMPG123}")
    set(HAVE_LIBMPG123 1)
    list(APPEND LIBS_DEPS ${LIBMPG123})
endif()

if(HAVE_AUDIO_HASH)
    set(EXTRA_SRC src/audiophash.cpp src/ph_fft.cpp)
    set(LIBS_DEPS ${LIBS_DEPS} sndfile samplerate vorbis vorbisenc ogg)
endif()

if(HAVE_VIDEO_HASH)
    list(APPEND EXTRA_SRC src/cimgffmpeg.cpp)
    set(LIBS_DEPS ${LIBS_DEPS} avcodec avformat swscale avutil)
endif()


if(PHASH_DYNAMIC)
	add_library(pHash SHARED ${MMAN_FILE} ${DIRENT_FILE} src/pHash.cpp ${EXTRA_SRC})
endif()

if(PHASH_STATIC)
	add_library(pHash STATIC ${MMAN_FILE} ${DIRENT_FILE} src/pHash.cpp ${EXTRA_SRC})
endif()

set_property(TARGET pHash PROPERTY VERSION "${pHash_VERSION_MAJOR}.${pHash_VERSION_MINOR}.${pHash_VERSION_PATCH}")
target_link_libraries(pHash ${LIBS_DEPS})

install(TARGETS pHash DESTINATION ${LIBDIR})
install(FILES src/pHash.h DESTINATION include)

if(PHASH_EXAMPLES)
set(EXAMPLEDIR examples)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

if(HAVE_IMAGE_HASH)
add_executable(TestDCT "${EXAMPLEDIR}/test_imagephash.cpp")
target_link_libraries (TestDCT pHash)
install(TARGETS TestDCT DESTINATION bin)

add_executable(TestRadish "${EXAMPLEDIR}/test_radish.cpp")
target_link_libraries (TestRadish pHash)
install(TARGETS TestRadish DESTINATION bin)

add_executable(TestMH "${EXAMPLEDIR}/test_mhimagehash.cpp")
target_link_libraries (TestMH pHash)
install(TARGETS TestMH DESTINATION bin)
endif()

if(HAVE_AUDIO_HASH)
add_executable(TestAudio "${EXAMPLEDIR}/test_audiophash.cpp")
target_link_libraries (TestAudio pHash)
install(TARGETS TestAudio DESTINATION bin)
endif()

if(HAVE_VIDEO_HASH)
add_executable(TestVideoHash "${EXAMPLEDIR}/test_dctvideohash.cpp")
target_link_libraries (TestVideoHash pHash)
install(TARGETS TestVideoHash DESTINATION bin)
endif()
endif() # PHASH_EXAMPLES

include (InstallRequiredSystemLibraries)

set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set (CPACK_PACKAGE_VERSION_MAJOR "${pHash_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${pHash_VERSION_MINOR}")

include (CPack)
