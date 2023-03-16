#include <AMReX_VisMF.H>
#include <AMReX_AsyncOut.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FPC.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_ParmParse.H>

#include <fstream>
#include <iomanip>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#include <openPMD/openpmd.hpp>

#include <regex>
namespace amrex {
  namespace openpmd_api {
    
    
    void AMReX_openPMDWriter::SetupPos(openPMD::ParticleSpecies& currSpecies,
				       const unsigned long long& np) const
    {
      //std::string options = "{}";
      //if (isBTD) options = "{ \"resizable\": true }";
      auto realType = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np}, m_openPMDDatasetOptions);
      auto idType = openPMD::Dataset(openPMD::determineDatatype< uint64_t >(), {np}, m_openPMDDatasetOptions);
      
      auto const positionComponents = helper::getParticlePositionComponentLabels();
      for( auto const& comp : positionComponents )
	{
	  currSpecies["position"][comp].resetDataset( realType );
	}
      
      auto const scalar = openPMD::RecordComponent::SCALAR;
      currSpecies["id"][scalar].resetDataset( idType );
    }

   
  }
}



