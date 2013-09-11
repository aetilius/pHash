
########### next target ###############

SET(pHash_LIB_SRCS
   pHash.cpp
)

if(HAVE_AUDIO_HASH)
	set(pHash_LIB_SRCS "${pHash_LIB_SRCS} audiophash.cpp ph_fft.c")
endif()

if(HAVE_VIDEO_HASH)
	set(pHash_LIB_SRCS "${pHash_LIB_SRCS} cimgffmpeg.cpp")
endif()

add_library(pHash SHARED ${pHash_LIB_SRCS})

target_link_libraries(pHash)

#set_target_properties(pHash PROPERTIES VERSION 4.2.0 SOVERSION 4)
install(TARGETS pHash DESTINATION lib)


########### install files ###############

install(FILES  pHash.h DESTINATION include)

