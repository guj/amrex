#ifndef AMREX_PLOTFILE_OPENPMD_API_H
#define AMREX_PLOTFILE_OPENPMD_API_H

#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileDataImpl.H>

#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include <string>
#include <memory>

#ifdef AMREX_USE_OPENPMD_API
#include <openPMD/openPMD.hpp>
#endif

namespace amrex
{
  namespace openpmd_api
  {
    //static std::vector<double> getReversedVec( const Real* v );    
    
    struct AMReX_BTDInfoField
    {
      AMReX_BTDInfoField(int sid,  const amrex::Geometry& fbs):
         snapshotID(sid), fullBTDSnapshot(fbs)
      {}

      AMReX_BTDInfoField(int sid,  const amrex::Geometry&& fbs) = delete;
      int snapshotID;
      const amrex::Geometry& fullBTDSnapshot;
    };


    void InitHandler(const std::string& filePath, bool isBTD);
    
    void SetStep(int ts);
    void CloseStep(int ts);    
    void CloseHandler();
    
    void WriteSingleLevel (//const std::string &plotfilename,
                                          const MultiFab &mf,
                                          const Vector<std::string> &varnames,
                                          const Geometry &geom,
                                          Real t,
                                          int level_step,
                                          const std::string &versionName = "HyperCLaw-V1.1",
                                          const std::string &levelPrefix = "Level_",
                                          const std::string &mfPrefix = "Cell",
                                          const Vector<std::string>& extra_dirs = Vector<std::string>(),
                                          AMReX_BTDInfoField* btd_info = NULL);

     void WriteMultiLevel (//const std::string &plotfilename,
                                         int nlevels,
                                         const Vector<const MultiFab*> &mf,
                                         const Vector<std::string> &varnames,
                                         const Vector<Geometry> &geom,
                                         Real time,
                                         const Vector<int> &level_steps,
                                         const Vector<IntVect> &ref_ratio,
                                         const std::string &versionName = "HyperCLaw-V1.1",
                                         const std::string &levelPrefix = "Level_",
                                         const std::string &mfPrefix = "Cell",
                                         const Vector<std::string>& extra_dirs = Vector<std::string>(),
                                         AMReX_BTDInfoField* btd_info = NULL);
    
    template<class PC>
    void WriteParticles(PC const&  pc);


    ////////////////////////////    
    struct AMReX_VarNameParser {
      AMReX_VarNameParser(std::string varname);

      void GetFieldNameModeInt(const std::string& varname);
      void GetMeshCompNames (int meshLevel);

      inline std::vector< std::string > getFieldComponentLabels ()
      {
        using vs = std::vector< std::string >;
        if (m_ThetaMode)
	  {
	    // if we write individual modes
	    vs const fieldComponents{"r", "t", "z"};
	    return fieldComponents;
	  }
	else
	  {
	    // if we just write reconstructed fields at theta=0 or are Cartesian
	    // note: 1D3V and 2D3V simulations still have 3 components for the fields 
	    vs const fieldComponents{"x", "y", "z"};
	    return fieldComponents;
	  }
      }

      inline std::vector< std::string >
      getFieldAxisLabels () const
      {
        using vs = std::vector< std::string >;
	if (0) {
	  /* TODO warpx specific 
	  // Fortran order of the index labels for the AMReX FArrayBox
#if defined(WARPX_DIM_1D_Z)
	  vs const axisLabels{"z"};  // z varies fastest in memory
#elif defined(WARPX_DIM_XZ)
	  vs const axisLabels{"x", "z"};  // x varies fastest in memory
#elif defined(WARPX_DIM_RZ)
	  // when we write individual modes of a field (default)
	  vs const circAxisLabels{"r", "z"};  // r varies fastest in memory
	  // if we just write reconstructed 2D fields at theta=0
	  vs const cartAxisLabels{"x", "z"};  // x varies fastest in memory
	  vs const axisLabels = m_ThetaMode ? circAxisLabels : cartAxisLabels;
#elif defined(WARPX_DIM_3D)
	  vs const axisLabels{"x", "y", "z"};  // x varies fastest in memory	
#else
#   error Unknown WarpX dimensionality.
#endif
	  // revert to C order (fastest varying index last)
	  return {axisLabels.rbegin(), axisLabels.rend()};	  
	  */
	} else {
	  // temporary resolution
#if AMREX_SPACEDIM == 1
	  vs const axisLabels{"z"};  
#elif AMREX_SPACEDIM == 2
	  vs const axisLabels{"x", "z"};  
#elif AMREX_SPACEDIM == 3
	  vs const axisLabels{"x", "y", "z"};  // x varies fastest in memory
#else
#  error Unable to label more than 3d
	  std::cout<<"TODO  ... SPACEDIM = 0? "<<std::endl;
	  vs const axisLabels{}
#endif
	  // revert to C order (fastest varying index last)
	  return {axisLabels.rbegin(), axisLabels.rend()};
	  
	}
      } // getFieldAxisLabels

      
      /////  members      
      std::string m_FieldName;
      std::string m_CompName;
      bool m_ThetaMode = false;
      int m_ModeIndex=-1;
      
    };


    ////////////////////////////
    /*
    class AMReX_ParticleCounter
    {
    public:
      //using ParticleContainer = typename WarpXParticleContainer::ContainerLike<amrex::PinnedArenaAllocator>;
      using ParticleContainer  = ContainerLike<amrex::PinnedArenaAllocator>;
      using ParticleIter = typename amrex::ParIter<0, 0, PIdx::nattribs, 0, amrex::PinnedArenaAllocator>;
      
      AMReX_ParticleCounter (ParticleContainer* pc);
      unsigned long GetTotalNumParticles () {return m_Total;}
      
      std::vector<unsigned long long> m_ParticleOffsetAtRank;
      std::vector<unsigned long long> m_ParticleSizeAtRank;
    private:

      void GetParticleOffsetOfProcessor (const long& numParticles,
					 unsigned long long& offset,
					 unsigned long long& sum)  const ;
      
      
      int m_MPIRank = 0;
      int m_MPISize = 1;

      unsigned long long m_Total = 0;
      
      std::vector<unsigned long long> m_ParticleCounterByLevel;
};
    */
    

    ////////////////////////////    
    class AMReX_openPMDWriter {
    public:      
      AMReX_openPMDWriter(const std::string& prefix,
			  openPMD::IterationEncoding ie,
			  std::string filetype,
			  std::string openpmdOptions);
      
      virtual ~AMReX_openPMDWriter();
      //std::vector<bool> fieldPMLdirections);

      virtual openPMD::Iteration GetIteration (int const iteration) const
      {
	// double check reached here
	return m_Series->writeIterations()[iteration];	
      }
      virtual void SetStep(int ts);
      virtual void CloseStep(int ts);
      virtual void Init(openPMD::Access access);
      
      virtual void WriteMesh(const std::vector<std::string>& varnames,
			     const amrex::Vector<const amrex::MultiFab*>& mf,
			     const amrex::Vector<amrex::Geometry>& geom,
			     int output_levels,
			     const Vector<int> &iteration,
			     //const int iteration,
			     const double time ) const;

      template<class PC>
      void DumpParticles(PC const& pc,
			 const std::string& name,
			 int iteration,
			 const amrex::Vector<int>& write_real_comp,
			 const amrex::Vector<int>& write_int_comp,
			 const amrex::Vector<std::string>& real_comp_names,
			 const amrex::Vector<std::string>&  int_comp_names
			 ) const;

      void SetupPos(openPMD::ParticleSpecies& species, const unsigned long long& np) const;


      template<class PC>
      void
      SetupRealProperties(PC const& pc, 
			  openPMD::ParticleSpecies& currSpecies,
			  const amrex::Vector<int>& write_real_comp,
			  const amrex::Vector<std::string>& real_comp_names,
			  const amrex::Vector<int>& write_int_comp,
			  const amrex::Vector<std::string>& int_comp_names,
			  const unsigned long long np) const;
      /*
      */
      template<class PIt>
      void
      SaveRealProperty (PIt& pti, //int, int,                                                                        
			openPMD::ParticleSpecies& currSpecies,
			unsigned long long offset,
			const amrex::Vector<int>& write_real_comp,
			const amrex::Vector<std::string>& real_comp_names,
			const amrex::Vector<int>& write_int_comp,
			const amrex::Vector<std::string>& int_comp_names) const;
      

      std::unique_ptr<openPMD::Series> m_Series = nullptr;
      std::string m_openPMDPrefix = std::string();

      int m_openPMDMinDigits = 6;
      std::string m_openPMDDatasetOptions = "{}"; 
      
    private:
      void CompSetup(int lev,
		     openPMD::Container< openPMD::Mesh >& meshes,
		     amrex::Geometry& full_geom,
		     const std::vector<std::string>& varnames,
		     const amrex::MultiFab* mf) const;

      void CompStorage(int lev,
		       openPMD::Container< openPMD::Mesh >& meshes,
		       amrex::Geometry& full_geom,
		       const std::vector<std::string>& varnames,
		       const amrex::MultiFab* mf) const;


      void SetupFields (openPMD::Container< openPMD::Mesh >& meshes,			
			amrex::Geometry& full_geom) const;

      void SetupMeshComp (openPMD::Mesh& mesh,					     
			  amrex::Geometry& full_geom,
			  amrex::MultiFab const& mf,
			  const AMReX_VarNameParser& varName) const;

     
      void GetFileName(std::string& filepath);

      template<class PIt>
      void StoreAoS_Real(PIt& pti,
			 openPMD::ParticleSpecies& currSpecies,
			 const amrex::Vector<int>& write_real_comp,
			 const amrex::Vector<std::string>& real_comp_names,			 
			 unsigned long long offset) const;

      
      template<class PIt>
      void
      StoreAoS_Int(PIt& pti,
		   openPMD::ParticleSpecies& currSpecies,
		   const amrex::Vector<int>& write_int_comp,
		   const amrex::Vector<std::string>& int_comp_names,
		   unsigned long long offset) const;
      

      template<class PIt>
      void StoreSoAReal(PIt& pti,
			openPMD::ParticleSpecies& currSpecies,
			const amrex::Vector<int>& write_real_comp,
			const amrex::Vector<std::string>& real_comp_names,
			unsigned long long offset)  const;


      template<class PIt>
      void StoreSoAInt(PIt& pti,
		       openPMD::ParticleSpecies& currSpecies,
		       const amrex::Vector<int>& write_int_comp,
		       const amrex::Vector<std::string>& int_comp_names,
		       unsigned long long offset) const;
      
      int m_CurrentStep = -1;
      openPMD::IterationEncoding  m_openPMDEncoding = openPMD::IterationEncoding::fileBased;
      std::string m_openPMDFileType = "bp";
      std::string m_openPMDSeriesOptions = std::string();
    };

    ////////////////////////////    
    class AMReX_openPMDWriterBTD final : public AMReX_openPMDWriter
    {
    public:
      AMReX_openPMDWriterBTD(const std::string& prefix,
			     openPMD::IterationEncoding ie,
			     std::string filetype,
			     std::string openpmdOptions);

      ~AMReX_openPMDWriterBTD() override = default;
      openPMD::Iteration GetIteration (int const iteration) const override
      {
	// double check reached here
	return m_Series->iterations[iteration];
      }

      //void SetStep(int ts);
      void CloseStep(int ts) override;

      void Init(openPMD::Access access) override;

      void WriteMesh(const std::vector<std::string>& varnames,
		     const amrex::Vector<const amrex::MultiFab*>& mf,
		     const amrex::Vector<amrex::Geometry>& geom,
		     int output_levels,
		     const Vector<int> &iteration,
		     //const iteration,
		     const double time ) const override;

      //void WriteMesh() override;
    private:

    };


    ////////////////////////////
    class AMReX_openPMDHandler {
    public:
      AMReX_openPMDHandler(const std::string& prefix = std::string(),
			   bool isBTD=false);
      ~AMReX_openPMDHandler();      
      
      std::unique_ptr<AMReX_openPMDWriter> m_Writer = nullptr;
      
    private:
      void CreateWriter(const std::string& prefix = std::string());
      bool m_IsBTD;
      
    }; // class AMReX_openPMDHandler
    

    //static std::unique_ptr< AMReX_openPMDHandler > m_OpenPMDHandler = nullptr;
        
  }// name space openpmd_api
} // namespace amrex


#include <AMReX_PlotFileUtilOPENPMD_PTLImpl.H>

#endif // AMREX_PLOTFILE_OPENPMD_API_H
