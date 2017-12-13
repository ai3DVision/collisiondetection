#ifndef CTCDNARROWPHASE_H
#define CTCDNARROWPHASE_H

#include "RetrospectiveDetection.h"

// Alec: *NarrowPhase classes are really just namespaces (these functions are
// effectively static and there are no members). record_min_t should just be a
// parameter...
class CTCDNarrowPhase : public NarrowPhase
{
public:
	virtual void findCollisions(const History &h, const std::set<std::pair<VertexFaceStencil, double> > &candidateVFS, const std::set<std::pair<EdgeEdgeStencil, double> > &candidateEES,
		std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees);

//private:
	bool checkVFS(const History &h, VertexFaceStencil & vfs, double eta);
	bool checkEES(const History &h, EdgeEdgeStencil & ees, double eta);
};

#endif
