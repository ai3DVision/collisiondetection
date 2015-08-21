#include "KDOPBroadPhase.h"
#include <set>
#include "History.h"
#include "Mesh.h"
#include <iostream>

using namespace std;
using namespace Eigen;

KDOPBroadPhase::KDOPBroadPhase()
{
  DOPaxis.push_back(Vector3d(1.0, 0, 0));
  DOPaxis.push_back(Vector3d(0, 1.0, 0));
  DOPaxis.push_back(Vector3d(0, 0, 1.0));
  DOPaxis.push_back(Vector3d(1.0, 1.0, 0));
  DOPaxis.push_back(Vector3d(1.0, -1.0, 0));
  DOPaxis.push_back(Vector3d(0, 1.0, 1.0));
  DOPaxis.push_back(Vector3d(0, 1.0, -1.0));
  DOPaxis.push_back(Vector3d(1.0, 0, 1.0));
  DOPaxis.push_back(Vector3d(1.0, 0, -1.0));
  DOPaxis.push_back(Vector3d(1.0, 1.0, 1.0));
  DOPaxis.push_back(Vector3d(1.0, 1.0, -1.0));
  DOPaxis.push_back(Vector3d(1.0, -1.0, 1.0));
  DOPaxis.push_back(Vector3d(1.0, -1.0, -1.0));
  if(K>DOPaxis.size())
    exit(0);

  for(int i=0; i<(int)DOPaxis.size(); i++)
    {
      DOPaxis[i] /= DOPaxis[i].norm();
    }
}

void KDOPBroadPhase::findCollisionCandidates(const History &h, const Mesh &m, double outerEta, set<VertexFaceStencil> &vfs, set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts)
{
	vfs.clear();
	ees.clear();
	KDOPNode *tree = buildKDOPTree(h, m, outerEta);
	intersect(tree, tree, m, vfs, ees, fixedVerts);
	delete tree;
}

KDOPNode *KDOPBroadPhase::buildKDOPTree(const History &h, const Mesh &m, double outerEta)
{
  vector<KDOPNode *> leaves;
  vector<bool> referenced(m.vertices.size()/3,false);
  for(int i=0; i<(int)m.faces.cols(); i++)
  {
    KDOPLeafNode *node = new KDOPLeafNode;
    node->type = KDOPLeafNode::PRIMITIVE_TYPE_FACE;
    node->id = i;
    int verts[3];
    for(int j=0; j<3; j++)
    {
      verts[j] = m.faces.coeff(j, i);
      referenced[m.faces.coeff(j, i)] = true;
    }
    for(int j=0; j<K; j++)
    {
      node->mins[j] = std::numeric_limits<double>::infinity();
      node->maxs[j] = -std::numeric_limits<double>::infinity();
    }
    for(int j=0; j<3; j++)
    {
      for(vector<HistoryEntry>::const_iterator it = h.getVertexHistory(verts[j]).begin(); it != h.getVertexHistory(verts[j]).end(); ++it)
      {
        for(int k=0; k<K; k++)
        {
          node->mins[k] = min(it->pos.dot(DOPaxis[k]) - outerEta, node->mins[k]);
          node->maxs[k] = max(it->pos.dot(DOPaxis[k]) + outerEta, node->maxs[k]);
        }
      }
    }
    leaves.push_back(node);		
  }
  // also add leaf nodes for unreferenced vertices
  for(int i = 0;i<(int)m.vertices.size()/3;i++)
  {
    if(!referenced[i])
    {
      KDOPLeafNode *node = new KDOPLeafNode;
      node->type = KDOPLeafNode::PRIMITIVE_TYPE_VERTEX;
      node->id = i;
      // initialize extents
      for(int j=0; j<K; j++)
      {
        node->mins[j] = std::numeric_limits<double>::infinity();
        node->maxs[j] = -std::numeric_limits<double>::infinity();
      }
      // loop through history
      for(vector<HistoryEntry>::const_iterator it = h.getVertexHistory(i).begin(); it != h.getVertexHistory(i).end(); ++it)
      {
        for(int k=0; k<K; k++)
        {
          node->mins[k] = min(it->pos.dot(DOPaxis[k]) - outerEta, node->mins[k]);
          node->maxs[k] = max(it->pos.dot(DOPaxis[k]) + outerEta, node->maxs[k]);
        }
      }
      leaves.push_back(node);		
    }
  }
  return buildKDOPInterior(leaves);
}

KDOPNode *KDOPBroadPhase::buildKDOPInterior(vector<KDOPNode *> &children)
{
	int nchildren = children.size();
	assert(nchildren > 0);
	if(nchildren == 1)
		return children[0];

	KDOPInteriorNode *node = new KDOPInteriorNode;
	
	for(int i=0; i<K; i++)	
	{
		node->mins[i] = std::numeric_limits<double>::infinity();
		node->maxs[i] = -std::numeric_limits<double>::infinity();
	}

	for(vector<KDOPNode *>::iterator it = children.begin(); it != children.end(); ++it)
	{
	  for(int j=0; j<K; j++)
	    {
	      node->mins[j] = min((*it)->mins[j], node->mins[j]);
	      node->maxs[j] = max((*it)->maxs[j], node->maxs[j]);
	    }
	}
	double lengths[K];
	for(int i=0; i<K; i++)
	{
	  lengths[i] = node->maxs[i] - node->mins[i];
	}
	int greatest = -1;
	double greatestlen = 0;
	for(int i=0; i<K; i++)
	  {
	    if(lengths[i] > greatestlen)
	      {
		greatestlen = lengths[i];
		greatest = i;
	      }
	  }
	node->splitaxis = greatest;

	sort(children.begin(), children.end(), NodeComparator(node->splitaxis));

	vector<KDOPNode *> left;
	int child=0;
	for(; child<nchildren/2; child++)
		left.push_back(children[child]);
	node->left = buildKDOPInterior(left);
	vector<KDOPNode *> right;
	for(; child<nchildren; child++)
		right.push_back(children[child]);
	node->right = buildKDOPInterior(right);
	return node;
}

void KDOPBroadPhase::intersect(KDOPNode *left, KDOPNode *right, const Mesh &m, std::set<VertexFaceStencil> &vfs, std::set<EdgeEdgeStencil> &ees, const std::set<int> &fixedVerts)
{
	for(int axis=0; axis<K; axis++)
	{
	  if(left->maxs[axis] < right->mins[axis] ||
	     right->maxs[axis] < left->mins[axis])
	    return;		
	}
	if(!left->isLeaf())
	{
		KDOPInteriorNode *ileft = (KDOPInteriorNode *)left;
		intersect(ileft->left, right, m, vfs, ees, fixedVerts);
		intersect(ileft->right, right, m, vfs, ees, fixedVerts);	
	}
	else if(!right->isLeaf())
	{
		KDOPInteriorNode *iright = (KDOPInteriorNode *)right;
		intersect(left, iright->left, m, vfs, ees, fixedVerts);
		intersect(left, iright->right, m, vfs, ees, fixedVerts);
	}
	else
	{
		KDOPLeafNode *lleft = (KDOPLeafNode *)left;
		KDOPLeafNode *lright = (KDOPLeafNode *)right;
    if(
      lleft->type == KDOPLeafNode::PRIMITIVE_TYPE_FACE &&
      lright->type == KDOPLeafNode::PRIMITIVE_TYPE_FACE)
    {
		  if(
          m.neighboringFaces(lleft->id, lright->id))
		    return;

		  // 6 vertex-face and 9 edge-edge
		  for(int i=0; i<3; i++)
		  {
		  	bool alllfixed = true;
		  	bool allrfixed = true;
		  	alllfixed = alllfixed && fixedVerts.count(m.faces.coeff(i, lleft->id));
		  	allrfixed = allrfixed && fixedVerts.count(m.faces.coeff(i, lright->id));
		  	for(int j=0; j<3; j++)
		  	{
		  		alllfixed = alllfixed && fixedVerts.count(m.faces.coeff(j, lright->id));
		  		allrfixed = allrfixed && fixedVerts.count(m.faces.coeff(j, lleft->id));
		  	}
		  	if(!alllfixed)
		  		vfs.insert(VertexFaceStencil(m.faces.coeff(i, lleft->id), m.faces.coeff(0, lright->id), m.faces.coeff(1, lright->id), m.faces.coeff(2, lright->id),lright->id));
		  	if(!allrfixed)
		  		vfs.insert(VertexFaceStencil(m.faces.coeff(i, lright->id), m.faces.coeff(0, lleft->id), m.faces.coeff(1, lleft->id), m.faces.coeff(2, lleft->id),lleft->id));
		  	for(int j=0; j<3; j++)
		  	{
		  		bool allefixed = true;
		  		allefixed = allefixed && fixedVerts.count(m.faces.coeff(i, lleft->id));
		  		allefixed = allefixed && fixedVerts.count(m.faces.coeff((i+1)%3, lleft->id));
		  		allefixed = allefixed && fixedVerts.count(m.faces.coeff(j, lright->id));
		  		allefixed = allefixed && fixedVerts.count(m.faces.coeff((j+1)%3, lright->id));
		  		if(!allefixed)
		  			ees.insert(
              EdgeEdgeStencil(
                m.faces.coeff(i, lleft->id), 
                m.faces.coeff((i+1)%3, lleft->id), 
                m.faces.coeff(j, lright->id), 
                m.faces.coeff((j+1)%3, lright->id),
                lleft->id,
                lright->id,
                (i+2)%3,
                (j+2)%3));
		  	}
		  }		
    }else if(
      lleft->type == KDOPLeafNode::PRIMITIVE_TYPE_VERTEX &&
      lright->type == KDOPLeafNode::PRIMITIVE_TYPE_VERTEX)
    {
      // both vertices... I guess ignore these...
    }else
    {
      const auto & addvf = [&](const int v, const int f)
      {
		  	bool allffixed = true;
		  	allffixed = allffixed && fixedVerts.count(v);
		  	for(int j=0; j<3; j++)
		  	{
		  		allffixed = allffixed && fixedVerts.count(m.faces.coeff(j,f));
		  	}
		  	if(!allffixed)
		  		vfs.insert(VertexFaceStencil(v, m.faces.coeff(0,f), m.faces.coeff(1,f), m.faces.coeff(2,f),f));
      };
      if(lleft->type == KDOPLeafNode::PRIMITIVE_TYPE_FACE)
      {
        assert(lright->type == KDOPLeafNode::PRIMITIVE_TYPE_VERTEX);
        addvf(lright->id,lleft->id);
      }else
      {
        assert(lleft->type == KDOPLeafNode::PRIMITIVE_TYPE_VERTEX);
        assert(lright->type == KDOPLeafNode::PRIMITIVE_TYPE_FACE);
        addvf(lleft->id,lright->id);
      }
    }

	}
}
