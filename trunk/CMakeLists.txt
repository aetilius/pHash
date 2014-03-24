cmake_minimum_required (VERSION 2.8)

include(ExternalProject)

project(pHash)

set(pHash_VERSION_MAJOR 0)
set(pHash_VERSION_MINOR 9.7)

OPTION(WITH_AUDIO_HASH "Audio hash support" 0)
OPTION(WITH_VIDEO_HASH "Video hash support" 0)

set(HAVE_IMAGE_HASH 1)

if(WITH_AUDIO_HASH)
	set(HAVE_AUDIO_HASH 1)
endif()

if(WITH_VIDEO_HASH)
	set(HAVE_VIDEO_HASH 1)
endif()

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/src/pHash.h.cmake ${CMAKE_SOURCE_DIR}/src/pHash.h)


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

add_subdirectory(bindings)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BUILD_TYPE}")

set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined ${CMAKE_SHARED_LINKER_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--no-undefined ${CMAKE_EXE_LINKER_FLAGS}")

add_library(pHash SHARED ${MMAN_FILE} ${DIRENT_FILE} src/pHash.cpp )


set(LIBS_DEPS 			png jpeg tiff pthread rt
				)


if(HAVE_AUDIO_HASH)
	add_library(pHash src/audiophash.cpp)
	set(LIBS_DEPS "${LIBS_DEPS} sndfile samplerate vorbis vorbisenc ogg")
endif()

target_link_libraries(pHash	${LIBS_DEPS})

install(TARGETS pHash DESTINATION lib)

set(EXAMPLEDIR examples)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
if(HAVE_IMAGE_HASH)
add_executable(TestDCT "${EXAMPLEDIR}/test_imagephash.cpp")
target_link_libraries (TestDCT pHash)
install(TARGETS TestDCT DESTINATION bin)

add_executable(TestRadish "${EXAMPLEDIR}/test_radish.cpp")
target_link_libraries (TestRadish pHash)
install(TARGETS TestRadish DESTINATION bin)

#add_executable(TestBMB "${EXAMPLEDIR}/test_bmbimagehash.cpp")
#target_link_libraries (TestBMB pHash)
#install(TARGETS TestBMB DESTINATION bin)

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

include (InstallRequiredSystemLibraries)

set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set (CPACK_PACKAGE_VERSION_MAJOR "${pHash_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${pHash_VERSION_MINOR}")

include (CPack)
