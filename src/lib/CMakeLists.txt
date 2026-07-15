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
    # The imported target is used so that both the compile AND link flags
    # propagate transitively to every consumer of Q2DL (PUBLIC).
    # The SET_TARGET_PROPERTIES(LINK_FLAGS) alternative is deprecated and, for
    # static libraries, CMake silently ignores LINK_FLAGS.
    target_link_libraries(Q2DL PUBLIC OpenMP::OpenMP_Fortran)
ENDIF(USE_OPENMP)

IF(USE_MPI)
    target_compile_definitions(Q2DL PUBLIC USE_MPI)
    target_link_libraries(Q2DL PUBLIC MPI::MPI_Fortran)
ENDIF(USE_MPI)
