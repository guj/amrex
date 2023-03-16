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


namespace amrex
{  
  namespace openpmd_api
  {
    /*global*/std::unique_ptr< AMReX_openPMDHandler > m_OpenPMDHandler = nullptr;
    
    void InitHandler(const std::string& prefix, bool isBTD)
    {
      std::string filePath {""};
      if (prefix.size() == 0)
	{	  
	  ParmParse pp;
	  pp.query("openpmd_directory", filePath);	  
	}
      else {
	filePath = prefix;
      }
      
      if (m_OpenPMDHandler == nullptr)
	m_OpenPMDHandler.reset(new AMReX_openPMDHandler(filePath, isBTD));
      else if (m_OpenPMDHandler->m_Writer != nullptr)
	if (m_OpenPMDHandler->m_Writer->m_openPMDPrefix !=  filePath)	
	  m_OpenPMDHandler.reset(new AMReX_openPMDHandler(filePath, isBTD));
      // already using the directory, no action needed
    }

    void CloseHandler()
    {
      if (m_OpenPMDHandler == nullptr)
	return;

      m_OpenPMDHandler.reset(nullptr);
    }
    
    void SetStep(int ts)
    {
      if ((m_OpenPMDHandler == nullptr) || (m_OpenPMDHandler->m_Writer == nullptr))
	return;

      m_OpenPMDHandler->m_Writer->SetStep(ts);
    }

    void CloseStep(int ts)
    {
      if ((m_OpenPMDHandler == nullptr) || (m_OpenPMDHandler->m_Writer == nullptr))
	return;

      m_OpenPMDHandler->m_Writer->CloseStep(ts);      
    }
    
    void WriteSingleLevel (//const std::string &plotfilename,
                                          const MultiFab &mf,
                                          const Vector<std::string> &varnames,
                                          const Geometry &geom,
                                          Real t,
                                          int level_step,
                                          const std::string &versionName, 
                                          const std::string &levelPrefix,
                                          const std::string &mfPrefix,
                                          const Vector<std::string>& extra_dirs,
                                          AMReX_BTDInfoField* btd_info )
    {
      /*
      ParallelDescriptor::Barrier();
      
      amrex::Print()<<"===> To write single level in openpmd\n";
      if (m_OpenPMDHandler == nullptr)
	m_OpenPMDHandler.reset(new AMReX_openPMDHandler(plotfilename));
      else if (m_OpenPMDHandler->m_openPMDFileInput !=  plotfilename)	
	m_OpenPMDHandler.reset(new AMReX_openPMDHandler(plotfilename));

      // silly test
	//openPMD::Series series = openPMD::Series("1_structure.bp", openPMD::Access::CREATE);
      
        openPMD::ParticleSpecies electrons = m_OpenPMDHandler->m_Writer->GetIteration(1).particles["electrons"];

	openPMD::Record mass = electrons["mass"];
	openPMD::RecordComponent mass_scalar = mass[openPMD::RecordComponent::SCALAR];

	openPMD::Dataset dataset = openPMD::Dataset(openPMD::Datatype::DOUBLE, openPMD::Extent{1});
	mass_scalar.resetDataset(dataset);

	electrons["position"]["x"].resetDataset(dataset);
	electrons["position"]["x"].makeConstant(20.0);
	electrons["positionOffset"]["x"].resetDataset(dataset);
	electrons["positionOffset"]["x"].makeConstant(22.0);
      */
      
    }
    
    void WriteMultiLevel (//const std::string &plotfilename,
                                         int nlevels,
                                         const Vector<const MultiFab*> &mf,
                                         const Vector<std::string> &varnames,
                                         const Vector<Geometry> &geom,
                                         Real time,
                                         const Vector<int> &level_steps,
                                         const Vector<IntVect> &ref_ratio,
                                         const std::string &versionName,
                                         const std::string &levelPrefix,
                                         const std::string &mfPrefix,
                                         const Vector<std::string>& extra_dirs,
                                         AMReX_BTDInfoField* btd_info )
    {
      if ((m_OpenPMDHandler == nullptr) || (m_OpenPMDHandler->m_Writer == nullptr))
	return;

      m_OpenPMDHandler->m_Writer->WriteMesh(varnames,
					    mf, //amrex::GetVecOfConstPtrs(mf),
					    geom,
					    nlevels,
					    level_steps,
					    time);
    }
  } // namespace openpmd_api
} // namespace amrex

