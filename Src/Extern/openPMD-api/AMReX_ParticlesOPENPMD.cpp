#include <AMReX_VisMF.H>
#include <AMReX_AsyncOut.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FPC.H>
#include <AMReX_FabArrayUtility.H>
 
#include <fstream>
#include <iomanip>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#include <openPMD/openpmd.hpp>

namespace amrex {
  void WriteParticleOpenPMD (const amrex::Vector<ParticleDiag>& particle_diags,
		             const bool use_pinned_pc,
                             AMReX_BTDInfoParticle* btd_ptl_info = NULL)
  {

  }
}
