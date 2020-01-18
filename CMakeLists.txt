cmake_minimum_required (VERSION 3.7.2)

include(ExternalProject)

project(pHash)

set(CMAKE_PROJECT_VERSION 1.0.0)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_RULE_MESSAGES OFF)


set(CMAKE_BUILD_TYPE  RelWithDebInfo CACHE STRING "build type" FORCE)
set(USE_IMAGE_HASH ON CACHE BOOLEAN "compile image hashes" FORCE)
set(USE_AUDIO_HASH ON CACHE BOOLEAN "compile audio hashes" FORCE)
set(USE_VIDEO_HASH ON CACHE BOOLEAN "compile video hash" FORCE)
set(USE_TEXT_HASH ON CACHE BOOLEAN "compile text hash" FORCE)

if (USE_IMAGE_HASH OR USE_VIDEO_HASH OR USE_TEXT_HASH)
  configure_file(src/pHash.h.cmake src/pHash.h)
endif()

if (USE_AUDIO_HASH)
  configure_file(src/audiophash.h.cmake src/audiophash.h)
endif()
if (USE_VIDEO_HASH)
  configure_file(src/cimgffmpeg.h.cmake src/cimgffmpeg.h)
endif()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS filesystem system program_options)
if (Boost_FOUND)
  if (NOT Boost_LIBRARIES)
	set(Boost_LIBRARIES Boost::filesystem Boost::system Boost::program_options)
  endif()
elseif()
  message(WARNING "No Boost libraries found.")
endif()

if (USE_IMAGE_HASH)
  message(STATUS "Using Image Hash")
  find_library(libpng png)
  if(libpng-NOTFOUND)
	message(WARNING "no libpng-dev found")
  else()
	set(USE_PNG 1)
  endif()
  find_library(libjpeg jpeg)
  if(libjpeg-NOTFOUND)
	message(WARNING "no libjpeg found")
  else()
	set(USE_JPEG 1)
  endif()
  find_library(libtiff tiff)
  if (libtiff-NOTFOUND)
	message(WARNING  "no libtiff found")
  else()
	set(USE_TIFF 1)
  endif()
  list(APPEND LIB_SRCS src/pHash.cpp src/bmbhash.cpp)
  list(APPEND LIB_DEPS ${libpng} ${libjpeg} ${libtiff})
endif()

if (USE_AUDIO_HASH)
  message(STATUS "Using Audio Hash")
  find_library(libsnd sndfile)
  if (libsnd-NOTFOUND)
	message(FATAL_ERROR "no libsndfile found")
  endif()
  find_library(libsamplerate samplerate)
  if (libsamplerate-NOTFOUND)
	message(FATAL_ERROR "no libsamplerate found")
  endif()
  find_library(libmpg123 mpg123)
  if (libmpg123-NOTFOUND)
	message(WARNING "no libmpg123 found")
  endif()
  find_library(libvorbis vorbis)
  if (libvorbis-NOTFOUND)
	message(FATAL_ERROR "no libvorbis found")
  endif()
  find_library(libogg ogg)
  if (libogg-NOTFOUND)
	message(FATAL_ERROR "no libogg found")
  endif()
  find_library(libm m)
  if (libm-NOTFOUND)
	message(FATAL_ERROR "no math library libm found")
  endif()
  
  list(APPEND LIB_SRCS src/audiophash.cpp  src/ph_fft.cpp)
  list(APPEND LIB_DEPS ${libm} ${libsnd} ${libsamplerate} ${libmpg123} ${libvorbis} ${libogg})
endif()

if (USE_VIDEO_HASH)
  message(STATUS "Using Video Hash")
  find_library(libavformat avformat)
  if (libavformat-NOTFOUND)
	message(FATAL_ERROR "no libavformat found")
  endif()
  find_library(libavcodec avcodec)
  if (libavcodec-NOTFOUND)
	message(FATAL_ERROR "no libavcodec found")
  endif()
  find_library(libswscale swscale)
  if (libswscale-NOTFOUND)
	message(FATAL_ERROR "no libswscale found")
  endif()
  find_library(libavutil avutil)
  if (libavutil-NOTFOUND)
	message(FATAL_ERROR "no libavutil found")
  endif()
  find_library(libswresample swresample)
  if (libswresample-NOTFOUND)
	message(FATAL_ERROR "no libswresample found")
  endif()
  list(APPEND LIB_SRCS src/pHash.cpp src/cimgffmpeg.cpp)
  list(APPEND LIB_DEPS ${libavformat} ${libavcodec} ${libswscale} ${libavutil} ${libswresample})
endif()

if (USE_TEXT_HASH)
  list(APPEND LIB_SRCS src/pHash.cpp)
  message(STATUS "Using Text Hash")
endif()

message(STATUS "phash library sources: " ${LIB_SRCS})
message(STATUS "phash library deps: "  ${LIB_DEPS})

add_library(pHash SHARED ${LIB_SRCS})
set_property(TARGET pHash PROPERTY VERSION ${CMAKE_PROJECT_VERSION})
target_include_directories(pHash PUBLIC src)
target_link_libraries(pHash	${LIBS_DEPS})

# example programs
if (USE_IMAGE_HASH)
add_executable(imghash examples/imghash.cpp)
target_include_directories(imghash PUBLIC src ${Boost_INCLUDE_DIRECTORIES})
target_link_libraries(imghash pHash ${LIB_DEPS} ${Boost_LIBRARIES})
endif()

install(TARGETS pHash DESTINATION lib)
if (USE_IMAGE_HASH OR USE_VIDEO_HASH OR USE_TEXT_HASH)
  install(FILES src/pHash.h DESTINATION include)
endif()
if (USE_AUDIO_HASH)
  install(FILES src/audiophash.h DESTINATION include)
endif()
if (USE_IMAGE_HASH)
  install(TARGETS imghash  DESTINATION bin)
endif()


# tests
enable_testing()

if(USE_AUDIO_HASH)
  add_executable(Testfft tests/fft-test.cpp)
  target_include_directories(Testfft PRIVATE src)
  target_link_libraries(Testfft pHash ${LIB_DEPS})
  add_test(NAME test_fft COMMAND Testfft)

  add_executable(TestAudioCount tests/audiohash-countsamples-test.cpp)
  target_include_directories(TestAudioCount PRIVATE src)
  target_link_libraries(TestAudioCount pHash ${LIB_DEPS})
  add_test(NAME test-audio-count COMMAND TestAudioCount resources/audio.wav resources/audioclip.wav 59531 23889)

  add_executable(TestReadAudio tests/audiohash-readaudio-test.cpp)
  target_include_directories(TestReadAudio PRIVATE src)
  target_link_libraries(TestReadAudio pHash ${LIB_DEPS})
  add_test(NAME test-read-audio COMMAND TestReadAudio resources/audio.wav resources/audioclip.wav 59531 23889)

  add_executable(TestAudioHash tests/audiohash-compare-test.cpp)
  target_include_directories(TestAudioHash PRIVATE src)
  target_link_libraries(TestAudioHash pHash ${LIB_DEPS})
  add_test(NAME test-audio-hash COMMAND TestAudioHash resources/audio.wav resources/audioclip.wav 434 155 0.985026 64 280)
endif(USE_AUDIO_HASH)

if(USE_IMAGE_HASH)
  add_executable(Testdct tests/dct-imagehash-test.cpp)
  target_include_directories(Testdct PRIVATE src)
  target_link_libraries(Testdct pHash ${LIB_DEPS})
  add_test(Test_dct Testdct resources/img1.jpg resources/img1.jpg 0)
  add_test(Test_dct2 Testdct resources/img1.jpg resources/img2.jpg  2)
  add_test(Test_dct3 Testdct resources/img1.jpg resources/img3.jpg 31)
  add_test(Test_dct4 Testdct resources/img2.jpg resources/img3.jpg 33)

  add_executable(Testbmb tests/bmb-imagehash-test.cpp)
  target_include_directories(Testbmb PRIVATE src)
  target_link_libraries(Testbmb pHash ${LIB_DEPS})
  add_test(Test_bmb1 Testbmb resources/img1.jpg resources/img1.jpg 0)
  add_test(Test_bmb2 Testbmb resources/img1.jpg resources/img2.jpg 0)
  add_test(Test_bmb3 Testbmb resources/img1.jpg resources/img3.jpg 150)
  add_test(Test_bmb4 Testbmb resources/img2.jpg resources/img3.jpg 150)

  add_executable(Testmh tests/mh-imagehash-test.cpp)
  target_include_directories(Testmh PRIVATE src)
  target_link_libraries(Testmh pHash ${LIB_DEPS})
  add_test(Test_mh1 Testmh resources/img1.jpg resources/img1.jpg 0)
  add_test(Test_mh2 Testmh resources/img1.jpg resources/img2.jpg 78)
  add_test(Test_mh3 Testmh resources/img1.jpg resources/img3.jpg 305)
  add_test(Test_mh4 Testmh resources/img2.jpg resources/img3.jpg 299)

  add_executable(Testrh tests/radial-imagehash-test.cpp)
  target_include_directories(Testrh PRIVATE src)
  target_link_libraries(Testrh pHash ${LIB_DEPS})
  add_test(Test_rh1 Testrh resources/img1.jpg resources/img1.jpg 1.0000)
  add_test(Test_rh2 Testrh resources/img1.jpg resources/img2.jpg 0.999954)
  add_test(Test_rh3 Testrh resources/img1.jpg resources/img3.jpg 0.843873)
  add_test(Test_rh4 Testrh resources/img2.jpg resources/img3.jpg 0.843654)
endif(USE_IMAGE_HASH)

if(USE_VIDEO_HASH)
  add_executable(Testvideohash tests/dct-videohash-test.cpp)
  target_include_directories(Testvideohash PRIVATE src)
  target_link_libraries(Testvideohash pHash ${LIB_DEPS})
  add_test(Test_vh1 Testvideohash resources/video.mp4 2)
  add_test(Test_vh2 Testvideohash resources/videoclip.mp4 1)
endif(USE_VIDEO_HASH)

if (USE_TEXT_HASH)
  add_executable(Testtxthash tests/txthash-compare-test.cpp)
  target_include_directories(Testtxthash PRIVATE src)
  target_link_libraries(Testtxthash pHash ${LIB_DEPS})
  add_test(Test_txt1 Testtxthash resources/alice.txt resources/aliceparag.txt 1273 5 4)
endif(USE_TEXT_HASH)

