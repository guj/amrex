target_include_directories( amrex
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources( amrex
   PRIVATE
   AMReX_PlotFileUtilOPENPMD.H
   AMReX_PlotFileUtilOPENPMD.cpp
   AMReX_PlotFileOPENPMD.cpp
   AMReX_PlotFileOPENPMD_PTL.cpp	

   AMReX_ParticlesOPENPMD.H
 )
