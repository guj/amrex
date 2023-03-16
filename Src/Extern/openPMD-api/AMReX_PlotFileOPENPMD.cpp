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

    // from warpx Utils/RelativeCellPosition.cpp 
    std::vector< double >
    getRelativeCellPosition(amrex::MultiFab const& mf)
    {
      amrex::IndexType const idx_type = mf.ixType();
      
      std::vector< double > relative_position(AMREX_SPACEDIM, 0.0);
      // amrex::CellIndex::CELL means: 0.5 from lower corner for that index/direction
      // amrex::CellIndex::NODE means: at corner for that index/direction
      // WarpX::do_nodal means: all indices/directions on CellIndex::NODE
      for (int d = 0; d < AMREX_SPACEDIM; d++)
	{
	  if (idx_type.cellCentered(d))
            relative_position.at(d) = 0.5;
	}
      return relative_position;
    }

    
    static  std::vector<std::uint64_t>
    getReversedVec( const IntVect& v )
    {
      // Convert the IntVect v to and std::vector u
      std::vector<std::uint64_t> u = {
	AMREX_D_DECL(
		     static_cast<std::uint64_t>(v[0]),
		     static_cast<std::uint64_t>(v[1]),
		     static_cast<std::uint64_t>(v[2])
		     )
      };
      // Reverse the order of elements, if v corresponds to the indices of a
      // Fortran-order array (like an AMReX FArrayBox)
      // but u is intended to be used with a C-order API (like openPMD)
      std::reverse( u.begin(), u.end() );
      
      return u;
    }

    static std::vector<double> getReversedVec( const Real* v )
    {
      // Convert Real* v to and std::vector u
      std::vector<double> u = {
	AMREX_D_DECL(
		     static_cast<double>(v[0]),
		     static_cast<double>(v[1]),
		     static_cast<double>(v[2])
		     )
      };
      // Reverse the order of elements, if v corresponds to the indices of a
      // Fortran-order array (like an AMReX FArrayBox)
      // but u is intended to be used with a C-order API (like openPMD)
      
      std::reverse( u.begin(), u.end() );

      return u;
    }

    inline void
    setOpenPMDUnit ( openPMD::Mesh mesh, const std::string field_name )
    {
        if (field_name[0] == 'E'){
	  mesh.setUnitDimension({
	      {openPMD::UnitDimension::L,  1},
	      {openPMD::UnitDimension::M,  1},
	      {openPMD::UnitDimension::T, -3},
	      {openPMD::UnitDimension::I, -1},
	    });
	} else if (field_name[0] == 'B'){ 
            mesh.setUnitDimension({
		{openPMD::UnitDimension::M,  1},
		{openPMD::UnitDimension::I, -1},
		{openPMD::UnitDimension::T, -2}
	      });
	}else if (field_name[0] == 'j'){ // current
            mesh.setUnitDimension({
		{openPMD::UnitDimension::L, -2},
		{openPMD::UnitDimension::I,  1},
	      });
        } else if (field_name.substr(0,3) == "rho"){ // charge density
	  mesh.setUnitDimension({
	      {openPMD::UnitDimension::L, -3},
	      {openPMD::UnitDimension::I,  1},
	      {openPMD::UnitDimension::T,  1},
	    });
        }
    }
        
	    
    
    ////////////////////////////////////////
    //
    // Struct AMReX_VarNameParser
    //    parser var names to field and comp names
    //
    ////////////////////////////////////////

    AMReX_VarNameParser::AMReX_VarNameParser(std::string varname)
    :m_CompName(openPMD::MeshRecordComponent::SCALAR)
    {
	//auto [varname_no_mode, mode_index] = GetFieldNameModeInt(varname);
	GetFieldNameModeInt(varname);
	//bool var_in_theta_mode = mode_index != -1; // thetaMode or reconstructed Cartesian 2D slice    
	//std::string m_FieldName = varname_no_mode;
	
	// assume fields are scalar unless they match the following match of known vector fields  
    }

    void AMReX_VarNameParser::GetFieldNameModeInt (const std::string& varname)
    {
      // mode_index = -1 if varname isn't of form fieldName_mode_realOrImag
      // mode_index = 2 * mode - 1 + (realOrImag == 'imag')
      // in either case, there is a -1 in mode_index
      //int mode_index = -1;
      
      std::regex e_real_imag("(.*)_([0-9]*)_(real|imag)");
      std::smatch sm;
      std::regex_match(varname, sm, e_real_imag, std::regex_constants::match_default);
      
      if (sm.size() != 4 )
	{
	  m_ThetaMode = (m_ModeIndex != -1); // thetaMode or reconstructed Cartesian 2D slice
	  m_FieldName = varname;
	}
      else
	{
	  // sm = [varname, field_name, mode, real_imag]
	  int mode = std::stoi(sm[2]);
	  if (mode == 0) 
            m_ModeIndex = 0;
	  else
	    {
	      if (sm[3] == "imag") {
		m_ModeIndex += 1;
	      }
	      m_ModeIndex += 2 * mode;
	    }
	  m_ThetaMode = (m_ModeIndex != -1); // thetaMode or reconstructed Cartesian 2D slice
	  m_FieldName = std::string(sm[1]);
	}	
    }


    void AMReX_VarNameParser::GetMeshCompNames (int meshLevel) 
    {
      std::string varname = m_FieldName; 
      if (varname.size() >= 2u )
	{
	  std::string const varname_1st = varname.substr(0u, 1u); // 1st character
	  std::string const varname_2nd = varname.substr(1u, 1u); // 2nd character
	  
	  // Check if this field is a vector. If so, then extract the field name
	  
	  std::vector< std::string > const vector_fields = {"E", "B", "j"};
	  std::vector< std::string > const field_components = getFieldComponentLabels();
	  
	  for( std::string const& vector_field : vector_fields )
	    {
	      for( std::string const& component : field_components ) 
		{
		  if( vector_field.compare( varname_1st ) == 0 && component.compare( varname_2nd ) == 0 )		      
		    {
		      m_FieldName = varname_1st + varname.substr(2); // Strip component
		      m_CompName  = varname_2nd;
		    }
		}
	    }
	} 

      
      if ( 0 == meshLevel )
	return;
      
      m_FieldName += std::string("_lvl").append(std::to_string(meshLevel));
    }

     

    // create json options to pass to openpmd
    inline std::string
    getSeriesOptions (std::string const & operator_type,
                      std::map< std::string, std::string > const & operator_parameters,
                      std::string const & engine_type,
                      std::map< std::string, std::string > const & engine_parameters)
    {
      if (operator_type.empty() && engine_type.empty())
	return "{}";

      std::string options;
      std::string top_block;
      std::string end_block;
      std::string op_block;
      std::string en_block;
      
      std::string op_parameters;
      for (const auto& kv : operator_parameters) {
	if (!op_parameters.empty()) op_parameters.append(",\n");
	op_parameters.append(std::string(12, ' '))         /* just pretty alignment */
	  .append("\"").append(kv.first).append("\": ")    /* key */
	  .append("\"").append(kv.second).append("\""); /* value (as string) */
      }

      std::string en_parameters;
      for (const auto& kv : engine_parameters) {
	if (!en_parameters.empty()) en_parameters.append(",\n");
	en_parameters.append(std::string(12, ' '))         /* just pretty alignment */
	  .append("\"").append(kv.first).append("\": ")    /* key */
	  .append("\"").append(kv.second).append("\""); /* value (as string) */
      }
      
      top_block = R"END( 
{
  "adios2": {)END";
      end_block = R"END(
  }
})END";
    
      if (!operator_type.empty()) {
	op_block = R"END( 
    "dataset": { 
      "operators": [
        {
          "type": ")END";
	op_block += operator_type + "\"";
	if (!op_parameters.empty()) {
	  op_block += R"END(,
          "parameters": {
)END";		
	  op_block += op_parameters +
	    "\n          }";
	}
	op_block += R"END(
        }                                                                                                          
      ]    
    })END";
      
	if (!engine_type.empty() || !en_parameters.empty())
	  op_block += ",";
      }
      
      if (!engine_type.empty() || !en_parameters.empty())
	{
	  en_block = R"END( 
    "engine": {)END";	     
	  if (!engine_type.empty()) {
	    en_block += R"END(
      "type": ")END";
	    en_block += engine_type + "\"";	       
	    if(!en_parameters.empty())
	      en_block += ",";
	  }
	  if (!en_parameters.empty()) {
	    en_block += R"END( 
      "parameters": { 
)END";
	    en_block += en_parameters +
	      "\n      }";
	  }
	  en_block += R"END( 
    })END";
	}
      
      options = top_block + op_block + en_block + end_block;
      return options;    
    }


    ////////////////////////////////////////
    //
    // Class AMReX_openPMDHandler
    //    files are saved as prefix/openpmd.bp
    //
    ////////////////////////////////////////

    AMReX_openPMDHandler::AMReX_openPMDHandler(const std::string& prefix, // match to diag_name in warpx
					       bool isBTD)
      :m_Writer(nullptr),       
       m_IsBTD(isBTD)
    {
      CreateWriter(prefix);
    }

    AMReX_openPMDHandler::~AMReX_openPMDHandler ()
    {
      std::cout<<"... ::~handler"<<std::endl;
      if( m_Writer)
      {
          m_Writer.reset( nullptr );
      }
    }

    void AMReX_openPMDHandler::CreateWriter(const std::string& prefix)
    {
      ParmParse pp_prefix(prefix);
	
      // choose backend (e.g. ADIOS, ADIOS2 or HDF5). Default depends on openPMD-api configuration     
      std::string openpmd_backend {"default"};
      pp_prefix.query("openpmd_backend", openpmd_backend);
      
      std::string  openpmd_encoding {"f"};      
      bool encodingDefined = pp_prefix.query("openpmd_encoding", openpmd_encoding);
      openPMD::IterationEncoding encoding = openPMD::IterationEncoding::groupBased;
     
      if ( 0 == openpmd_encoding.compare("v") )
        encoding = openPMD::IterationEncoding::variableBased;
      else if ( 0 == openpmd_encoding.compare("g") )
        encoding = openPMD::IterationEncoding::groupBased;
      else if ( 0 == openpmd_encoding.compare("f") )
        encoding = openPMD::IterationEncoding::fileBased;


      if (m_IsBTD)
      {
        if ( ( openPMD::IterationEncoding::fileBased != encoding ) &&
             ( openPMD::IterationEncoding::groupBased != encoding ) )
        {
       	   std::string warnMsg = prefix+" Unable to support BTD with streaming. Using GroupBased ";
           encoding = openPMD::IterationEncoding::groupBased;
        }
      }

      auto lf_collect = [&](const char* key,
			    const std::string& parameter_tag,
			    std::string&  key_type,
			    std::map< std::string, std::string>& result)->void
      {
	//std::string key_type;
	pp_prefix.query(key, key_type);
	std::string const key_prefix = prefix + parameter_tag;
	ParmParse pp;
	auto entr = pp.getEntries(key_prefix);

	auto const prefix_len = key_prefix.size() + 1;
	for (std::string k : entr) {
	  std::string v;
	  pp.get(k.c_str(), v);
	  k.erase(0, prefix_len);
	  result.insert({k, v});
	}
      };

      std::string  operator_type;
      std::map< std::string, std::string > operator_parameters;
      lf_collect("adios2_operator.type", ".adios2_operator.parameters", operator_type, operator_parameters);

      std::string  engine_type;
      std::map< std::string, std::string > engine_parameters;      
      lf_collect("adios2_engine.type", ".adios2_engine.parameters", engine_type,  engine_parameters);      

      std::string options=getSeriesOptions(operator_type, operator_parameters,
					   engine_type, engine_parameters);
      
      if (m_IsBTD)
	m_Writer = std::make_unique<AMReX_openPMDWriterBTD>(prefix, encoding, openpmd_backend, options);
      else
	m_Writer = std::make_unique<AMReX_openPMDWriter>(prefix, encoding, openpmd_backend, options);

      pp_prefix.query("file_min_digits", m_Writer->m_openPMDMinDigits);
      
    } // CreateWriter()

    
    ////////////////////////////////////////
    //
    // Classs AMReX_openPMDWriter
    //
    ////////////////////////////////////////
    AMReX_openPMDWriter::AMReX_openPMDWriter (const std::string& prefix,
					      openPMD::IterationEncoding ie,
					      std::string filetype,
					      std::string options)
      :m_openPMDPrefix(prefix),
       m_openPMDEncoding(ie),
       m_openPMDFileType(filetype),
       m_openPMDSeriesOptions(options)
                                              //std::vector<bool> fieldPMLdirections // warpx specific
    {
if( m_openPMDFileType == "default" )
#if openPMD_HAVE_ADIOS2==1
    m_openPMDFileType = "bp";
#elif openPMD_HAVE_ADIOS1==1
    m_openPMDFileType = "bp";
#elif openPMD_HAVE_HDF5==1
    m_openPMDFileType = "h5";
#else
    m_openPMDFileType = "json";
#endif
    }

    AMReX_openPMDWriter::~AMReX_openPMDWriter ()
    {
      std::cout<<" .. ::~writer "<<std::endl;
      if( m_Series )
      {
          m_Series->flush();
          m_Series.reset( nullptr );
      }
    }

    void AMReX_openPMDWriter::SetStep(int ts)
    {
      m_CurrentStep = ts;
      Init(openPMD::Access::CREATE);
    }

    void AMReX_openPMDWriter::CloseStep(int ts)
    {
      bool callClose = true;

      if (m_Series) {
	GetIteration(m_CurrentStep).close();
      }      
    }

    void AMReX_openPMDWriter::Init(openPMD::Access access)
    {
      std::string filepath = m_openPMDPrefix;
      GetFileName(filepath);

      if ( m_openPMDEncoding == openPMD::IterationEncoding::fileBased )
        m_Series = nullptr;
      else if ( m_Series != nullptr )
        return;

      if (amrex::ParallelDescriptor::NProcs() > 1)
	{
#if defined(AMREX_USE_MPI)
        m_Series = std::make_unique<openPMD::Series>(
						     filepath, access,
						     amrex::ParallelDescriptor::Communicator(),
						     m_openPMDSeriesOptions
						     );
#else
        amrex::Abort(Utils::TextMsg::Err("AMReX did not build with MPI support!"));
#endif
	}
      else
	{
	  m_Series = std::make_unique<openPMD::Series>(filepath, access, m_openPMDSeriesOptions);
	}
      
      m_Series->setIterationEncoding( m_openPMDEncoding );

      m_Series->setMeshesPath( "fields" );
      // conform to ED-PIC extension of openPMD
      
      uint32_t const openPMD_ED_PIC = 1u;
      m_Series->setOpenPMDextension( openPMD_ED_PIC );
      // meta info
      
      m_Series->setSoftware( "AMReX", amrex::Version() );            
    }
   
    void AMReX_openPMDWriter::CompSetup(int lev,
					openPMD::Container< openPMD::Mesh >& meshes,
					amrex::Geometry& full_geom,
					const std::vector<std::string>& varnames,					
				        const amrex::MultiFab* curr_mf) const
    {
      int const ncomp = curr_mf->nComp();
      for ( int icomp=0; icomp<ncomp; icomp++ )
	{
	  std::string const & varname = varnames[icomp];
	  amrex::openpmd_api::AMReX_VarNameParser curr(varname);
	  curr.GetMeshCompNames( lev ); 
	  {
	    std::cout<<"Level:  "<<lev<<"  Domain:  "<<full_geom.Domain()<<std::endl;;
	    if (curr.m_CompName == openPMD::MeshRecordComponent::SCALAR)
	      {
		if ( ! meshes.contains(curr.m_FieldName) )
		  {
		    auto mesh = meshes[curr.m_FieldName];
		    SetupMeshComp(  mesh, full_geom, *curr_mf, curr);
		  }
	      }
	    else
	      {
		auto mesh = meshes[curr.m_FieldName];
		if ( ! mesh.contains(curr.m_CompName) )
		  SetupMeshComp(  mesh, full_geom, *curr_mf, curr );		
	      }
	  }
	} // icomp setup loop 
    }

    void AMReX_openPMDWriter::CompStorage(int lev,
					  openPMD::Container< openPMD::Mesh >& meshes,
					  amrex::Geometry& full_geom,
					  const std::vector<std::string>& varnames,
					  const amrex::MultiFab* curr_mf) const
    {
      int const ncomp = curr_mf->nComp();
      amrex::Box const & global_box = full_geom.Domain();
      
      for ( int icomp=0; icomp<ncomp; icomp++ )
	{
	  std::string const & varname = varnames[icomp];
	  amrex::openpmd_api::AMReX_VarNameParser curr(varname);
	  curr.GetMeshCompNames( lev );
	  
	  auto mesh = meshes[curr.m_FieldName];
	  auto mesh_comp = mesh[curr.m_CompName];
	  
	  for( amrex::MFIter mfi(*curr_mf); mfi.isValid(); ++mfi )
            {
	      amrex::FArrayBox const& fab = (*curr_mf)[mfi];
	      amrex::Box const& local_box = fab.box();
	      
	      // Determine the offset and size of this chunk
	      amrex::IntVect const box_offset = local_box.smallEnd() - global_box.smallEnd();
	      auto chunk_offset = getReversedVec( box_offset );
	      auto chunk_size = getReversedVec( local_box.size() );
	      
	      if (curr.m_ThetaMode)
		{
		  chunk_offset.emplace(chunk_offset.begin(), curr.m_ModeIndex);
		  chunk_size.emplace(chunk_size.begin(), 1);
		}
#ifdef AMREX_USE_GPU
	      if (fab.arena()->isManaged() || fab.arena()->isDevice())
		{
		  amrex::BaseFab<amrex::Real> foo(local_box, 1, amrex::The_Pinned_Arena());
		  std::shared_ptr<amrex::Real> data_pinned(foo.release());
		  amrex::Gpu::dtoh_memcpy_async(data_pinned.get(), fab.dataPtr(icomp), local_box.numPts()*sizeof(amrex::Real));
		  // intentionally delayed until before we .flush(): amrex::Gpu::streamSynchronize();
		  mesh_comp.storeChunk(data_pinned, chunk_offset, chunk_size);
		}
	      else
#endif
		{
		  amrex::Real const *local_data = fab.dataPtr(icomp);
		  mesh_comp.storeChunkRaw(local_data, chunk_offset, chunk_size);				       
		}
            }
        } // icomp store loop     

    }
    
    void AMReX_openPMDWriter::WriteMesh(const std::vector<std::string>& varnames,
					const amrex::Vector<const amrex::MultiFab*>& mf,
					const amrex::Vector<amrex::Geometry>& geom,
					int output_levels,
					const Vector<int> &iteration,
					//const int iteration,
					const double time ) const

    {
      openPMD::Iteration series_iteration = GetIteration(m_CurrentStep);
      series_iteration.open();

      auto meshes = series_iteration.meshes;
      series_iteration.setTime( time );
            
      if ( varnames.size()==0 ) return;

      for (int lev=0; lev < output_levels; lev++)
	{
	  amrex::Geometry full_geom = geom[lev];
	  
	  if ( 0 == lev )
	    SetupFields(meshes, full_geom);
	  
	  CompSetup(lev, meshes, full_geom, varnames, mf[lev]);
	  CompStorage(lev, meshes, full_geom, varnames, mf[lev]);
#ifdef AMREX_USE_GPU
	  amrex::Gpu::streamSynchronize();
#endif
	  
	  m_Series->flush();	  
      } // for lev loop 
    }

    void AMReX_openPMDWriter::GetFileName(std::string& filepath)
    {
      if (filepath.size() == 0)
	filepath.append(".");
	  
      filepath.append("/");
      // transform paths for Windows
#ifdef _WIN32
      filepath = detail::replace_all(filepath, "/", "\\");
#endif

      std::string filename = "openpmd";
      //
      // OpenPMD supports timestepped names
      //                                                                                                                                           
      if (m_openPMDEncoding == openPMD::IterationEncoding::fileBased)
	{
	  std::string fileSuffix = std::string("_%0") + std::to_string(m_openPMDMinDigits) + std::string("T");
	  filename = filename.append(fileSuffix);
	}
      filename.append(".").append(m_openPMDFileType);
      filepath.append(filename);
      //return filename;      
    }

    void AMReX_openPMDWriter::SetupFields (openPMD::Container< openPMD::Mesh >& meshes,	
					   amrex::Geometry& full_geom) const
    {
      // meta data for ED-PIC extension
      auto const period = full_geom.periodicity(); // TODO double-check: is this the proper global bound or of some level?
      
      std::vector<std::string> fieldBoundary(6, "reflecting");
      std::vector<std::string> particleBoundary(6, "absorbing");
      fieldBoundary.resize(AMREX_SPACEDIM * 2);
      particleBoundary.resize(AMREX_SPACEDIM * 2);

      /*
	  TODO m_fieldPMLdirection is warpx specific  
      for (auto i = 0u; i < fieldBoundary.size() / 2u; ++i)
	if (m_fieldPMLdirections.at(i))
	  fieldBoundary.at(i) = "open";
      */
      
      for (auto i = 0u; i < fieldBoundary.size() / 2u; ++i)
	if (period.isPeriodic(i)) {
	  fieldBoundary.at(2u * i) = "periodic";
	  fieldBoundary.at(2u * i + 1u) = "periodic";
	  particleBoundary.at(2u * i) = "periodic";
	  particleBoundary.at(2u * i + 1u) = "periodic";
	}
      
      /* TODO  warpx specific
      meshes.setAttribute("fieldSolver", []() {
	switch (WarpX::electromagnetic_solver_id) {
	case ElectromagneticSolverAlgo::Yee :
	  return "Yee";
	case ElectromagneticSolverAlgo::CKC :
	  return "CK";
	case ElectromagneticSolverAlgo::PSATD :
	  return "PSATD";
	default:
	  return "other";
	}
      }())
	;
      */
      meshes.setAttribute("fieldBoundary", fieldBoundary);
      meshes.setAttribute("particleBoundary", particleBoundary);
      /* TODO warpx specific
      meshes.setAttribute("currentSmoothing", []() {
          if (WarpX::use_filter) return "Binomial";
          else return "none";
      }());
      if (WarpX::use_filter)
          meshes.setAttribute("currentSmoothingParameters", []() {
              std::stringstream ss;
              ss << "period=1;compensator=false";
#if (AMREX_SPACEDIM >= 2)
              ss << ";numPasses_x=" << WarpX::filter_npass_each_dir[0];
#endif
#if defined(WARPX_DIM_3D)
              ss << ";numPasses_y=" << WarpX::filter_npass_each_dir[1];
              ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[2];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
              ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[1];
#elif defined(WARPX_DIM_1D_Z)
              ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[0];
#endif
              std::string currentSmoothingParameters = ss.str();
              return currentSmoothingParameters;
          }());
      meshes.setAttribute("chargeCorrection", []() {
          if (WarpX::do_dive_cleaning) return "hyperbolic"; // TODO or "spectral" or something? double-check                                   
          else return "none";
      }());
       if (WarpX::do_dive_cleaning)
          meshes.setAttribute("chargeCorrectionParameters", "period=1");
      */
    }

    void AMReX_openPMDWriter::SetupMeshComp (openPMD::Mesh& mesh,					     
					     amrex::Geometry& full_geom,
					     amrex::MultiFab const& mf,
					     const AMReX_VarNameParser& varName
					     ) const
    {
      auto mesh_comp = mesh[varName.m_CompName];
      amrex::Box const & global_box = full_geom.Domain();
      auto global_size = getReversedVec(global_box.size() );
      amrex::Print()<<" TODO double check mesh comps (scalar vs not scalar) name="<<varName.m_CompName<<" global[0]="<<global_size[0]<<std::endl;
      // - Grid spacing
      std::vector<double> const grid_spacing = getReversedVec(full_geom.CellSize());
      
      // - Global offset
      std::vector<double> const global_offset = getReversedVec(full_geom.ProbLo());
      /*
	TODO warpx specific
#if defined(WARPX_DIM_RZ)
      auto & warpx = WarpX::GetInstance();
      if (var_in_theta_mode) {
	global_size.emplace(global_size.begin(), warpx.ncomps);
      }
#endif
      */
      // - AxisLabels
      
      //std::vector<std::string> axis_labels = getFieldAxisLabels(varName.m_ThetaMode);
      std::vector<std::string> axis_labels = varName.getFieldAxisLabels();
      
      // Prepare the type of dataset that will be written
      openPMD::Datatype const datatype = openPMD::determineDatatype<amrex::Real>();
      auto const dataset = openPMD::Dataset(datatype, global_size);
      mesh.setDataOrder(openPMD::Mesh::DataOrder::C);
      if (varName.m_ThetaMode) {
        mesh.setGeometry("thetaMode");
	/*
	  TODO  warpx specific
        mesh.setGeometryParameters("m=" + std::to_string(WarpX::n_rz_azimuthal_modes) + ";imag=+");
	*/
      }
      mesh.setAxisLabels(axis_labels);
      mesh.setGridSpacing(grid_spacing);
      mesh.setGridGlobalOffset(global_offset);
      mesh.setAttribute("fieldSmoothing", "none");
      mesh_comp.resetDataset(dataset);
      
      setOpenPMDUnit( mesh, varName.m_FieldName );
      
      // TODO: getRelativeCellPosition(mf) is from WarpX
      auto relative_cell_pos = getRelativeCellPosition(mf);     // AMReX Fortran index order    
      std::reverse( relative_cell_pos.begin(), relative_cell_pos.end() ); // now in C order
      mesh_comp.setPosition( relative_cell_pos );      
    }
    
   
					    
    ////////////////////////////////////////
    //
    // Classs AMReX_openPMDWriter<BTD>
    //
    ////////////////////////////////////////
    AMReX_openPMDWriterBTD::AMReX_openPMDWriterBTD(const std::string& prefix,
						   openPMD::IterationEncoding ie,
						   std::string filetype,
						   std::string options)
    //std::vector<bool> fieldPMLdirections // warpx specific
      :AMReX_openPMDWriter(prefix, ie, filetype, options)
    {
      m_openPMDDatasetOptions="{ \"resizable\": true }";
      amrex::Print()<<" TODO: make sure dataset option resizable IS used for BTD\n";
    }

    void AMReX_openPMDWriterBTD::CloseStep(int ts)
    {
      std::cout<<" TODO.. BTD ..  close step "<<std::endl;
      /*
      // default close is true
      bool callClose = true;
      
      // close BTD file only when isLastBTDFlush is true
      if (isBTD and !isLastBTDFlush) callClose = false;
      if (callClose) {
        if (m_Series) {
	  GetIteration(m_CurrentStep, isBTD).close();
        }	
      }
      */
    }

    void AMReX_openPMDWriterBTD::Init(openPMD::Access access)
    {
      std::cout<<"TODO.. BTD ..  Init "<<std::endl;
      if (m_Series != nullptr)
	return;
      AMReX_openPMDWriter::Init(access);
    }

    void AMReX_openPMDWriterBTD::WriteMesh(const std::vector<std::string>& varnames,
					   const amrex::Vector<const amrex::MultiFab*>& mf,
					   const amrex::Vector<amrex::Geometry>& geom,
					   int output_levels,
					   const Vector<int> &iteration,
					   //const int iteration,
					   const double time ) const
   
    {
      std::cout<<"TODO ... BTD.. fields .."<<std::endl;
    }

  } // namespace openpmd_api
} // namespace amrex

