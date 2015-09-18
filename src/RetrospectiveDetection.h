#ifndef RETROSPECTIVEDETECTION_H
#define RETROSPECTIVEDETECTION_H

#include "Stencils.h"
#include "History.h"
#include <set>

class Mesh;

class BroadPhase
{
public:
	virtual void findCollisionCandidates(
      const History &h, 
      const Mesh &m, 
      double outerEta, 
      const std::set<int> &fixedVerts,
      std::set<VertexFaceStencil> &vfs, 
      std::set<EdgeEdgeStencil> &ees) = 0;
	virtual ~BroadPhase() {}
};

class NarrowPhase
{
public:
  // Whether to record earliest time of intereference (or conversely whether to
  // return early on any interference)
  bool record_min_t;
  NarrowPhase():record_min_t(false){};
	virtual void findCollisions(const History &h, const std::set<std::pair<VertexFaceStencil, double> > &candidateVFS, const std::set<std::pair<EdgeEdgeStencil, double> > &candidateEES,
		std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees) = 0;
	virtual ~NarrowPhase() {}
};

#endif
