#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include <AMReX_AmrParticles.H>

#include "mypc.H"
#include "trilinear_deposition_K.H"

using namespace amrex;

struct TestParams {
    int nx;
    int ny;
    int nz;
    int nlevs;
    int max_grid_size;
    int nppc;
    bool verbose;
  
    int nplotfile=1;
};

void checkMFBox(const TestParams& parms,
		Vector<const MultiFab*> outputMF)		
{
      for (int lev=0; lev < parms.nlevs; lev++)
	{
	  auto curr_mf = outputMF[lev];
	  int const ncomp = curr_mf->nComp();
	  std::cout<<" checking boxes, ncomp="<<ncomp<<std::endl;

	  for ( int icomp=0; icomp<ncomp; icomp++ )
	    {	      
	      for( amrex::MFIter mfi(*curr_mf); mfi.isValid(); ++mfi )
		{
		  amrex::FArrayBox const& fab = (*curr_mf)[mfi];
		  amrex::Box const& local_box = fab.box();
		  
		  std::cout<<"  .. checking:   icomp="<<icomp<< " local box="<<local_box<<std::endl;
		  std::cout<<"   "<<*(fab.dataPtr())<<std::endl;
		}
	    }
	}
}



void testParticleMesh (TestParams& parms, int nghost)
{
    Vector<IntVect> rr(parms.nlevs-1);
    for (int lev = 1; lev < parms.nlevs; lev++)
        rr[lev-1] = IntVect(AMREX_D_DECL(2,2,2));

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box base_domain(domain_lo, domain_hi);

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;

    Vector<Geometry> geom(parms.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < parms.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(parms.nlevs);
    Vector<DistributionMapping> dm(parms.nlevs);

    Box domain = base_domain;
    IntVect size = IntVect(AMREX_D_DECL(parms.nx, parms.ny, parms.nz));
    for (int lev = 0; lev < parms.nlevs; ++lev)
    {
        ba[lev].define(domain);
        ba[lev].maxSize(parms.max_grid_size);
        dm[lev].define(ba[lev]);
        domain.grow(-size/4);   // fine level cover the middle of the coarse domain
        domain.refine(2);
    }

    Vector<MultiFab> density1(parms.nlevs);
    Vector<MultiFab> density2(parms.nlevs);
    for (int lev = 0; lev < parms.nlevs; lev++) {
        density1[lev].define(ba[lev], dm[lev], 1, nghost);
        density1[lev].setVal(0.0);
        density2[lev].define(ba[lev], dm[lev], 1, nghost);
        density2[lev].setVal(0.0);
    }

    MyParticleContainer myPC(geom, dm, ba, rr);
    myPC.SetVerbose(false);

    int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
    amrex::Print() << "Total number of particles    : " << num_particles << '\n' << '\n';

    bool serialize = true;
    int iseed = 451;
    double mass = 10.0;

    MyParticleContainer::ParticleInitData pdata = {{mass}, {}, {}, {}};
    myPC.InitRandom(num_particles, iseed, pdata, serialize);

    //
    // Here we provide an example of one way to call ParticleToMesh
    //
    amrex::ParticleToMesh(myPC, GetVecOfPtrs(density1), 0, parms.nlevs-1,
        [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& p,
                              amrex::Array4<amrex::Real> const& rho,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) noexcept
        {
            ParticleInterpolator::Linear interp(p, plo, dxi);

            interp.ParticleToMesh(p, rho, 0, 0, 1,
                [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& part, int comp)
                {
                    return part.rdata(comp);  // no weighting
                });
        });

    //
    // Here we provide an example of another way to call ParticleToMesh
    //
    int start_part_comp = 0;
    int start_mesh_comp = 0;
    int        num_comp = 1;

    amrex::ParticleToMesh(myPC,GetVecOfPtrs(density2),0,parms.nlevs-1,
                          TrilinearDeposition{start_part_comp,start_mesh_comp,num_comp});

    //
    // Now write the output from each into separate plotfiles for comparison
    //

    Vector<std::string> varnames;
    varnames.push_back("density");

    Vector<std::string> particle_varnames;
    particle_varnames.push_back("mass");

    Vector<int> level_steps;
    level_steps.push_back(0);
    level_steps.push_back(0);

    int output_levs = parms.nlevs;

    Vector<const MultiFab*> outputMF(output_levs);
    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputMF[lev] = &density1[lev];
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }
    /*
    WriteMultiLevelPlotfile("plt00000_v1", output_levs, outputMF,
                            varnames, geom, 0.0, level_steps, outputRR);
    myPC.WritePlotFile("plt00000_v1", "particle0", particle_varnames);
    
    for (int lev = 0; lev < output_levs; ++lev) {
        outputMF[lev] = &density2[lev];
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }
    WriteMultiLevelPlotfile("plt00000_v2", output_levs, outputMF,
                            varnames, geom, 0.0, level_steps, outputRR);
    myPC.WritePlotFile("plt00000_v2", "particle0", particle_varnames);
    */

    checkMFBox(parms, outputMF);

    // call count ptls to prepare ahead of time
    myPC.CountParticles();

    std::string fname = "";
    bool isBTD = false;
    openpmd_api::InitHandler(fname, isBTD);
    
    for (int ts = 0; ts < parms.nplotfile; ts++)
      {
	openpmd_api::SetStep(ts);
	openpmd_api::WriteParticles(myPC);
	
	char name[512]; 
	std::snprintf(name, sizeof name, "plotfile_%05d",  ts);

#ifdef PLOT_FILES_With_1_VAR  // if myPC.H is init with (1) 
	// write mesh with plotfile 
	WriteMultiLevelPlotfile(name, output_levs, outputMF,
				varnames, geom, 0.0, level_steps, outputRR);
	// write ptl with plotfile, partticle_varnames has size 1, so myPC.H has to define with (1)
	myPC.WritePlotFile("plt00000_v1", "particle0", particle_varnames);
#endif
	openpmd_api::WriteMultiLevel(//fname,
				     parms.nlevs,
				     outputMF, // amrex::GetVecOfConstPtrs(testField.m_mf),
				     varnames, // testField.m_Varnames,
				     geom, //  testField.m_Geom,
				     0.0, // testField.m_Time,
				     level_steps, //testField.m_Level_steps,
				     outputRR //testField.m_Ref_ratio
				     );
	
	openpmd_api::CloseStep(ts);
	
	amrex::Print()<<"Timestep: "<<ts<<" done \n";
      }

    openpmd_api::CloseHandler();
}




void testMeshOnly (TestParams& parms, int nghost)
{
  std::cout<<" ======= \n Testing mesh only!!!\n ======= "<<std::endl;
    Vector<IntVect> rr(parms.nlevs-1);
    for (int lev = 1; lev < parms.nlevs; lev++)
        rr[lev-1] = IntVect(AMREX_D_DECL(2,2,2));

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));

    const Box base_domain(domain_lo, domain_hi);

    Vector<BoxArray> grids(parms.nlevs);
    Vector<DistributionMapping> dm(parms.nlevs);

    Box domain = base_domain;
    IntVect size = IntVect(AMREX_D_DECL(parms.nx, parms.ny, parms.nz));
    for (int lev = 0; lev < parms.nlevs; ++lev)
    {
        grids[lev].define(domain);
        grids[lev].maxSize(parms.max_grid_size);
        dm[lev].define(grids[lev]);

	domain.grow(-size/4);   // fine level cover the middle of the coarse domain
	domain.refine(2);	
    }

    // Geom

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;
    
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    Vector<Geometry> geom(parms.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < parms.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<MultiFab> density1(parms.nlevs);
    for (int lev = 0; lev < parms.nlevs; lev++) {
        density1[lev].define(grids[lev], dm[lev], 1, nghost);
        density1[lev].setVal(0.0);
    }

    //
    // Now write the output from each into separate plotfiles for comparison
    //

    Vector<std::string> varnames;
    varnames.push_back("density");

    Vector<int> level_steps;
    level_steps.push_back(0);
    level_steps.push_back(0);

    int output_levs = parms.nlevs;

    Vector<const MultiFab*> outputMF(output_levs);
    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputMF[lev] = &density1[lev];
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }

    
    checkMFBox(parms, outputMF);
    
    std::string fname = "";
    // call count ptls to prepare ahead of time

    bool isBTD = false;
    openpmd_api::InitHandler(fname, isBTD);
    //openpmd_api::ok();
    
    for (int ts = 0; ts < parms.nplotfile; ts++)
      {
	openpmd_api::SetStep(ts);
	
	char name[512]; 
	std::snprintf(name, sizeof name, "plotfile_%05d",  ts);

	openpmd_api::WriteMultiLevel(//fname,
				     parms.nlevs,
				     outputMF, // amrex::GetVecOfConstPtrs(testField.m_mf),
				     varnames, // testField.m_Varnames,
				     geom, //  testField.m_Geom,
				     0.0, // testField.m_Time,
				     level_steps, //testField.m_Level_steps,
				     outputRR //testField.m_Ref_ratio
				     );
	
	openpmd_api::CloseStep(ts);
	
	amrex::Print()<<"Timestep: "<<ts<<" done \n";
      }

    openpmd_api::CloseHandler();

}




int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  ParmParse pp;

  TestParams parms;

  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);
  pp.get("nppc", parms.nppc);
  pp.get("nlevs", parms.nlevs);
  //pp.get("nplotfile", parms.nplotfile);
    
  if (parms.nppc < 1 && ParallelDescriptor::IOProcessor())
    amrex::Abort("Must specify at least one particle per cell");

  parms.verbose = false;
  pp.query("verbose", parms.verbose);

  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles per cell : ";
    std::cout << parms.nppc  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }

  
  int nghost = 0 ;
  std::cout<<"  TODO: RESOLVE!!!  if nghost=1  there is tile offset be  at -1  "<<std::endl;
#ifdef TEST_MESH_ONLY  
  testMeshOnly(parms, nghost);
#else
  testParticleMesh(parms, nghost);
#endif

  amrex::Finalize();
}
