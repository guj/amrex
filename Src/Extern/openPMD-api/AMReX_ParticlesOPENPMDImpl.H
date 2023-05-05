#ifndef AMREX_PTL_OPENPMD_API_IMPL_H
#define AMREX_PTL_OPENPMD_API_IMPL_H



template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::CountParticles()
{
  m_PtlCounter.m_MPISize = amrex::ParallelDescriptor::NProcs();
  m_PtlCounter.m_MPIRank = amrex::ParallelDescriptor::MyProc();

  m_PtlCounter.m_ParticleCounterByLevel.resize(this->finestLevel()+1);
  m_PtlCounter.m_ParticleOffsetAtRank.resize(this->finestLevel()+1);
  m_PtlCounter.m_ParticleSizeAtRank.resize(this->finestLevel()+1);

  auto lf_GetParticleOffsetOfProcessor = [&](const long& numParticles,
					     unsigned long long& offset,
					     unsigned long long& sum)  -> void
  {    
    std::vector<long> result(m_PtlCounter.m_MPISize,  0);
    amrex::ParallelGather::Gather (numParticles, result.data(), -1, amrex::ParallelDescriptor::Communicator());
    
    sum = 0;
    offset = 0;
    for (int i=0;  i<result.size();  i++) {
      sum +=  result[i];
      if ( i < m_PtlCounter.m_MPIRank)
	offset +=  result[i];
    }
  };
  
  for (auto currentLevel = 0; currentLevel <= this->finestLevel(); currentLevel++)
    {
      long numParticles = 0; // numParticles in this processor
      for (ParIter pti(*this, currentLevel); pti.isValid(); ++pti) {
	auto numParticleOnTile = pti.numParticles();
	numParticles += numParticleOnTile;
      }

      unsigned long long offset=0; // offset of this level
      unsigned long long sum=0; // numParticles in this level (sum from all processors)
      lf_GetParticleOffsetOfProcessor(numParticles, offset,  sum);

      m_PtlCounter.m_ParticleCounterByLevel[currentLevel] = sum;
      m_PtlCounter.m_ParticleOffsetAtRank[currentLevel] = offset;
      m_PtlCounter.m_ParticleSizeAtRank[currentLevel] = numParticles;
      
      // adjust offset, it should be numbered after particles from previous levels
      for (auto lv=0; lv<currentLevel; lv++)
	{
	  m_PtlCounter.m_ParticleOffsetAtRank[currentLevel] += m_PtlCounter.m_ParticleCounterByLevel[lv];
	}
      
      m_PtlCounter.m_Total += sum;
    }
}


#endif // AMREX_PTL_OPENPMD_API_H

