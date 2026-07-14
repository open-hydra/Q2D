file(GLOB_RECURSE Q2D_SRC
    CONFIGURE_DEPENDS
    *.f90
)

add_library(Q2DL ${Q2D_SRC} ${FLINT_SRC})
add_library(Q2DL::Q2DL ALIAS Q2DL)
target_link_libraries(Q2DL PRIVATE FiNeR::FiNeR ORION)
link_FLINT_dependencies(Q2DL)

target_include_directories(Q2DL INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    # Use the imported target so compile *and* link flags propagate
    # transitively to every consumer of Q2DL (PUBLIC).
    # The old SET_TARGET_PROPERTIES(LINK_FLAGS) approach is deprecated and
    # its LINK_FLAGS are silently ignored for static libraries by CMake.
    target_link_libraries(Q2DL PUBLIC OpenMP::OpenMP_Fortran)
ENDIF(USE_OPENMP)

IF(USE_MPI)
    target_compile_definitions(Q2DL PUBLIC USE_MPI)
    target_link_libraries(Q2DL PUBLIC MPI::MPI_Fortran)
ENDIF(USE_MPI)
