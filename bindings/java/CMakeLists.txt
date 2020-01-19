cmake_minimum_required(VERSION 3.7.0)
project(phash-jni-bindings VERSION 1.0.0)

set(USE_IMAGE_HASH ON CACHE BOOL "include image hash functions" FORCE)
set(USE_AUDIO_HASH ON CACHE BOOL "include audio hash functions" FORCE)
set(USE_VIDEO_HASH ON CACHE BOOL "include video hash functions" FORCE)
set(USE_TEXT_HASH ON CACHE BOOL  "include text hash functions" FORCE)

set(pHash-jni-SRCS phash-jni.cpp)
set(CMAKE_PROJECT_VERSION 1.0.0)

configure_file(pHash-jni.h.cmake pHash-jni.h)

find_library(pHashlib pHash PATHS ../../ NO_DEFAULT_PATH)
if (pHashlib-NOTFOUND)
  message(FATAL "pHash lib not found")
else()
  message(STATUS "found pHash lib!")
endif()

find_package(JNI 1.8 REQUIRED)
if (JNI_FOUND)
  message(STATUS "jni installation found!")
elseif (FindJNI-NOTFOUND)
  message(FATAL "no java jni installation found")
endif()

add_library(pHash-jni SHARED ${pHash-jni-SRCS})
set_property(TARGET pHash-jni PROPERTY CXX_STANDARD 11)
set_property(TARGET pHash-jni PROPERTY VERSION ${CMAKE_PROJECT_VERSION})

set_property(TARGET pHash-jni PROPERTY SKIP_BUILD_RPATH FALSE)
set_property(TARGET pHash-jni PROPERTY BUILD_WITH_INSTALL_RPATH FALSE)
set_property(TARGET pHash-jni PROPERTY INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set_property(TARGET pHash-jni PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

target_include_directories(pHash-jni PUBLIC ../../src ${JNI_INCLUDE_DIRS})
target_link_libraries(pHash-jni ${pHashlib} ${JNI_LIBRARIES})

install(TARGETS pHash-jni DESTINATION lib)

