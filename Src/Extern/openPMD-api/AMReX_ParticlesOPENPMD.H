#ifndef AMREX_PTL_OPENPMD_API_H
#define AMREX_PTL_OPENPMD_API_H

//#include <openPMD/openpmd.hpp>
  

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

AMReX_PtlCounter m_PtlCounter;
#endif // AMREX_PTL_OPENPMD_API_H


