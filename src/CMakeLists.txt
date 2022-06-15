FILE( GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cc )

list(REMOVE_ITEM SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/indbased_code.cc )
list(REMOVE_ITEM SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/runtests.cc )
list(REMOVE_ITEM SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/old_code.cc )

list(APPEND SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/ZIP/ZIP.cpp)

# add the executable
add_executable(${BEEPMBP} ${SRC_FILES})

# Set libraries compiled in debug mode to end in 'd'
SET_TARGET_PROPERTIES( ${BEEPMBP} PROPERTIES DEBUG_POSTFIX "d" )

SET( FDPAPI_URL "https://github.com/FAIRDataPipeline/cppDataPipeline/archive/refs/${CPPDATAPIPELINEREF}.zip" )

MESSAGE( STATUS "[FDPAPI]" )
MESSAGE( STATUS "\tFDPAPI Will be installed." )
MESSAGE( STATUS "\tURL: ${FDPAPI_URL}" )

include(FetchContent)
FetchContent_Declare(
    FDPApi
    URL ${FDPAPI_URL}
)
FetchContent_MakeAvailable(FDPApi)

find_package(ZLIB)

if(NOT ZLIB_FOUND)

    SET( zlib_URL "https://github.com/madler/zlib/archive/refs/tags/v1.2.12.zip" )

    MESSAGE( STATUS "[zlib]" )
    MESSAGE( STATUS "\tzlib Will be installed." )
    MESSAGE( STATUS "\tURL: ${zlib_URL}" )

    include(FetchContent)
    FetchContent_Declare(
        ZLib
        URL ${zlib_URL}
    )
    FetchContent_MakeAvailable(ZLib)

    TARGET_INCLUDE_DIRECTORIES( ${BEEPMBP} PUBLIC ${ZLIB_PUBLIC_HDRS} )
    TARGET_LINK_LIBRARIES( ${BEEPMBP} PUBLIC zlib )

else()
    TARGET_LINK_LIBRARIES( ${BEEPMBP} PUBLIC ZLIB::ZLIB )
endif()

TARGET_INCLUDE_DIRECTORIES( ${BEEPMBP} PUBLIC ${fdpapi_SOURCE_DIR}/include )
TARGET_INCLUDE_DIRECTORIES( ${BEEPMBP} PUBLIC ${easyzip_SOURCE_DIR}/include )
TARGET_INCLUDE_DIRECTORIES( ${BEEPMBP} PUBLIC ${CMAKE_SOURCE_DIR}/toml11 )

TARGET_LINK_LIBRARIES( ${BEEPMBP} PUBLIC ${MPI_C_LIBRARIES})
TARGET_LINK_LIBRARIES( ${BEEPMBP} PUBLIC fdpapi )


# Install the libraries
INSTALL( TARGETS ${BEEPMBP} 
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)