#include "CTCDNarrowPhase.h"
#include "History.h"
#include "CTCD.h"
#include <set>
#include <iostream>

using namespace std;

void CTCDNarrowPhase::findCollisions(const History &h, const set<pair<VertexFaceStencil, double> > &candidateVFS, const set<pair<EdgeEdgeStencil, double> > &candidateEES,
		set<VertexFaceStencil> &vfs, set<EdgeEdgeStencil> &ees)
{
	for(set<pair<VertexFaceStencil, double> >::const_iterator it = candidateVFS.begin(); it != candidateVFS.end(); ++it)
	{
    VertexFaceStencil vf =  it->first;
		if(checkVFS(h, vf, it->second))
			vfs.insert(vf);
	}
	for(set<pair<EdgeEdgeStencil, double> >::const_iterator it = candidateEES.begin(); it != candidateEES.end(); ++it)
	{
    EdgeEdgeStencil ee = it->first;
		if(checkEES(h, ee, it->second))
			ees.insert(ee);
	}
}

bool CTCDNarrowPhase::checkVFS(const History &h, VertexFaceStencil & vfs, double eta)
{
	vector<int> verts;
	verts.push_back(vfs.p);
	verts.push_back(vfs.q0);
	verts.push_back(vfs.q1);
	verts.push_back(vfs.q2);
	vector<StitchedEntry> sh;
	h.stitchCommonHistory(verts, sh);

	for(vector<StitchedEntry>::iterator it = sh.begin(); it != sh.end(); ++it)
	{
		vector<StitchedEntry>::iterator next = it+1;
		if(next == sh.end())
			break;

		double tinterval = next->time - it->time;

		double t;
		if(CTCD::vertexFaceCTCD(it->pos[0], it->pos[1], it->pos[2], it->pos[3],
						next->pos[0], next->pos[1], next->pos[2], next->pos[3],
						eta, t))
		{
      if(record_min_t)
      {
        vfs.min_t = std::min(vfs.min_t,t);
      }else
      {
        return true;
      }
		}
			
		// Vertex-face edges
		for(int edge=0; edge<3; edge++)
		{
			if(CTCD::vertexEdgeCTCD(it->pos[0], it->pos[1+(edge%3)], it->pos[1+ ((edge+1)%3)],
						next->pos[0], next->pos[1+(edge%3)], next->pos[1+ ((edge+1)%3)],
						eta, t))
			{
				if(record_min_t)
        {
          vfs.min_t = std::min(vfs.min_t,t);
        }else
        {
         return true;
        }
			}
		}
		// Vertex-face vertices
		for(int vert=0; vert<3; vert++)
		{
			if(CTCD::vertexVertexCTCD(it->pos[0], it->pos[1+vert],
						  next->pos[0], next->pos[1+vert],
						  eta, t))
			{
				if(record_min_t)
        {
          vfs.min_t = std::min(vfs.min_t,t);
        }else
        {
         return true;
        }
			}
		}
	}
  if(record_min_t)
  {
    return std::isfinite(vfs.min_t);
  }else
  {
    return false;
  }
}

bool CTCDNarrowPhase::checkEES(const History &h, EdgeEdgeStencil & ees, double eta)
{
	vector<int> verts;
	verts.push_back(ees.p0);
	verts.push_back(ees.p1);
	verts.push_back(ees.q0);
	verts.push_back(ees.q1);
	vector<StitchedEntry> sh;
	h.stitchCommonHistory(verts, sh);

	for(vector<StitchedEntry>::iterator it = sh.begin(); it != sh.end(); ++it)
	{
		vector<StitchedEntry>::iterator next = it+1;
		if(next == sh.end())
			break;
	
		double t;
		if(CTCD::edgeEdgeCTCD(it->pos[0], it->pos[1], it->pos[2], it->pos[3],
					next->pos[0], next->pos[1], next->pos[2], next->pos[3],
					      eta, t))
		{			
			if(record_min_t)
      {
        ees.min_t = std::min(ees.min_t,t);
      }else
      {
       return true;
      }
		}

		// Edge-edge vertices
		if(CTCD::vertexEdgeCTCD(it->pos[0], it->pos[2], it->pos[3], next->pos[0], next->pos[2], next->pos[3], eta, t))
		{
			if(record_min_t)
      {
        ees.min_t = std::min(ees.min_t,t);
      }else
      {
       return true;
      }
		}
		if(CTCD::vertexEdgeCTCD(it->pos[1], it->pos[2], it->pos[3], next->pos[1], next->pos[2], next->pos[3], eta, t))
		{
			if(record_min_t)
      {
        ees.min_t = std::min(ees.min_t,t);
      }else
      {
       return true;
      }
		}
		if(CTCD::vertexEdgeCTCD(it->pos[2], it->pos[0], it->pos[1], next->pos[2], next->pos[0], next->pos[1], eta, t))
		{
			if(record_min_t)
      {
        ees.min_t = std::min(ees.min_t,t);
      }else
      {
       return true;
      }
		}
		if(CTCD::vertexEdgeCTCD(it->pos[3], it->pos[0], it->pos[1], next->pos[3], next->pos[0], next->pos[1], eta, t))
		{
			if(record_min_t)
      {
        ees.min_t = std::min(ees.min_t,t);
      }else
      {
       return true;
      }
		}

		// edge vertex-edge vertex
		if(CTCD::vertexVertexCTCD(it->pos[0], it->pos[2], next->pos[0], next->pos[2], eta, t))
		{
			if(record_min_t)
      {
        ees.min_t = std::min(ees.min_t,t);
      }else
      {
       return true;
      }
		}
		if(CTCD::vertexVertexCTCD(it->pos[0], it->pos[3], next->pos[0], next->pos[3], eta, t))
		{
			if(record_min_t)
      {
        ees.min_t = std::min(ees.min_t,t);
      }else
      {
       return true;
      }
		}
		if(CTCD::vertexVertexCTCD(it->pos[1], it->pos[2], next->pos[1], next->pos[2], eta, t))
		{
			if(record_min_t)
      {
        ees.min_t = std::min(ees.min_t,t);
      }else
      {
       return true;
      }
		}
		if(CTCD::vertexVertexCTCD(it->pos[1], it->pos[3], next->pos[1], next->pos[3], eta, t))
		{
			if(record_min_t)
      {
        ees.min_t = std::min(ees.min_t,t);
      }else
      {
       return true;
      }
		}
	}
  if(record_min_t)
  {
    return std::isfinite(ees.min_t);
  }else
  {
    return false;
  }
}
