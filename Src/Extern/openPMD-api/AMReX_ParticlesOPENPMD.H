#ifndef AMREX_PTL_OPENPMD_API_H
#define AMREX_PTL_OPENPMD_API_H

//#include <openPMD/openpmd.hpp>
  
struct AMReX_BTDInfoParticle
{
  AMReX_BTDInfoParticle(bool isLast,  const amrex::Vector<int>&  past):
    isLastBTDFlush(isLast), totalParticlesFlushedAlready(past)
  {}
  
  AMReX_BTDInfoParticle(bool isLast,  const amrex::Vector<int>&&  past) = delete;
  
  bool isLastBTDFlush;
  const amrex::Vector<int>& totalParticlesFlushedAlready;
};


struct AMReX_PtlCounter
{
  int m_MPIRank = 0;
  int m_MPISize = 1;
  
  unsigned long long m_Total = 0;
  
  std::vector<unsigned long long> m_ParticleCounterByLevel;

  unsigned long GetTotalNumParticles () const { return m_Total;}
  
  std::vector<unsigned long long> m_ParticleOffsetAtRank;
  std::vector<unsigned long long> m_ParticleSizeAtRank;
};


void CountParticles();

/*
// TODO: warpx style signature and may not be the right one for AMReX
template <class F>
void WritePlotFileOPENPMD (const std::string& dir,
			   const std::string& name,
			   const Vector<int>& write_real_comp,
			   const Vector<int>& write_int_comp,
			   const Vector<std::string>& real_comp_names,
			   const Vector<std::string>&  int_comp_names,
			   F&& f) const;

template <class F>
void WriteBinaryParticleDataOPENPMD (const std::string& dir,
				     const std::string& name,
				     const Vector<int>& write_real_comp,
				     const Vector<int>& write_int_comp,
				     const Vector<std::string>& real_comp_names,
				     const Vector<std::string>&  int_comp_names,
				     F&& f, bool is_checkpoint=false) const;

*/
AMReX_PtlCounter m_PtlCounter;
#endif // AMREX_PTL_OPENPMD_API_H


