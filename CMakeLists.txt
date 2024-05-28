cmake_minimum_required(VERSION 3.12...3.25)
project(
    stepbystep-secret
    LANGUAGES CXX C
)

# find_package(Boost)

## Declare tetgen library
add_library(tetgen OBJECT
    ExternalLibs/tetgen/tetgen1.4.3/src/tetgen.cxx
    ExternalLibs/tetgen/tetgen1.4.3/src/predicates.cxx
    )
target_include_directories(tetgen PUBLIC ExternalLibs/tetgen/tetgen1.4.3/inc)
target_compile_definitions(tetgen PRIVATE TETLIBRARY)

## Declare the executable
add_executable(TriMultPoly
    src/main.cpp
    src/Algorithm/DMWT.cpp
    src/DataStructure/TriangleInfo.cpp
    src/DataStructure/Tile.cpp
    src/DataStructure/Hole.cpp
    src/DataStructure/EdgeInfo.cpp
    src/DataStructure/CurveInfo.cpp
    src/DataStructure/Boundary.cpp
    src/Utility/Point3.cpp
    src/Utility/Vector3.cpp
    )
set_target_properties(TriMultPoly PROPERTIES CXX_STANDARD 11)
target_include_directories(TriMultPoly PRIVATE src/Algorithm src/DataStructure src/Utility)
target_compile_definitions(TriMultPoly PRIVATE __forceinline=inline)

target_link_libraries(TriMultPoly PRIVATE tetgen)
