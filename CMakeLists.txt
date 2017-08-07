cmake_minimum_required(VERSION 3.1)
project (sm)

#REMEMBER: cmake -g "Unix Makefiles" .. -DCMAKE_BUILD_TYPE=Debug/Release

####### CMAKE SETTINGS ######

cmake_host_system_information(RESULT hname QUERY HOSTNAME)
#message(STATUS "Host name:" ${hname})

# Compiler options for the "debug" and "release" configurations
if(CMAKE_COMPILER_IS_GNUCXX)
    set (FLAGS_WARNINGS "-Wall -Wextra")
    set (FLAGS_ERRORS "-Wfatal-errors")
    set (FLAGS_DEBUG "-g3 -O0")
    
    if (hname STREQUAL "newturb")
        message(STATUS "Compiling on newturb")
        set (FLAGS_STATIC "")
    else()
        message(STATUS "NOT compiling on newturb")
        set (FLAGS_STATIC "-static")
    endif()

    set (FLAGS_OPTIMIZATION "-O3")
elseif(MSVC)
    set (FLAGS_WARNINGS "/W2") # More info: https://msdn.microsoft.com/en-us/library/thxezb7y.aspx
    set (FLAGS_ERRORS "")
    set (FLAGS_DEBUG "/Zi /Od")
    set (FLAGS_STATIC "")
    set (FLAGS_OPTIMIZATION "/O2")
else()
    MESSAGE( FATAL_ERROR "Compiler not supported. Only gcc and MSVC are currently supported (only for the makefile generation; the program should still compile fine)." )
endif()

# c++11
set (CMAKE_CXX_STANDARD 11)

# Debug flags:
set (CMAKE_CXX_FLAGS_DEBUG "${FLAGS_DEBUG} ${FLAGS_WARNINGS} ${FLAGS_ERRORS} ${FLAGS_STATIC}")
set (CMAKE_CXX_FLAGS_RELEASE "${FLAGS_OPTIMIZATION} ${FLAGS_WARNINGS} ${FLAGS_ERRORS} ${FLAGS_STATIC}")

# Generate a verbose makefile?
set (CMAKE_VERBOSE_MAKEFILE TRUE)

########## PROJECT ##########
# Look for all the source files
file(GLOB SOURCES "*.cpp")
file(GLOB HEADERS "*.hpp")

add_executable(sm ${SOURCES} ${HEADERS})

