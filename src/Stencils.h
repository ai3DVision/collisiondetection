#ifndef STENCILS_H
#define STENCILS_H

#include <cassert>
#include <utility>
#include <limits>

struct VertexFaceStencil
{
  int p;
  int q0,q1,q2;
  // Id of face (q0,q1,q2)
  int fq;
  bool isnew;
  // time of first recorded interference
  double min_t;

  VertexFaceStencil(int p_, int q0_, int q1_, int q2_, int fq_) : 
    isnew(true),
    min_t(std::numeric_limits<double>::infinity())
  {
    p = p_;
    fq=fq_;
    if(q0_ <= q1_)
    {
      if(q1_ <= q2_)
      {
        q0 = q0_;
        q1 = q1_;
        q2 = q2_;
      }
      else if(q0_ <= q2_)
      {
        q0 = q0_;
        q1 = q2_;
        q2 = q1_;
      }
      else
      {
        q0 = q2_;
        q1 = q0_;
        q2 = q1_;
      }
    }
    else if(q0_ <= q2_)
    {
      q0 = q1_;
      q1 = q0_;
      q2 = q2_;
    }
    else if(q1_ <= q2_)
    {
      q0 = q1_;
      q1 = q2_;
      q2 = q0_;
    }
    else
    {
      q0 = q2_;
      q1 = q1_;
      q2 = q0_;
    }

    assert(q0 <= q1 && q1 <= q2);
  }

  bool operator<(const VertexFaceStencil &other) const
  {
    if(p < other.p)
	return true;
    else if(p > other.p)
	return false;

    if(q0 < other.q0)
	return true;
    else if(q0 > other.q0)
	return false;

    if(q1 < other.q1)
	return true;
    else if(q1 > other.q1)
	return false;

    return q2 < other.q2;
  }
};

struct EdgeEdgeStencil
{
  int p0,p1;
  int q0,q1;
  // faces of (p0,p1) and (q0,q1)
  int fp,fq;
  // Corner across from (p0,p1) in fp and from (q0,q1) in fq
  int cp,cq;
  bool isnew;
  // time of first recorded interference
  double min_t;

  EdgeEdgeStencil(
    int p0_, int p1_, int q0_, int q1_,int fp_,int fq_,int cp_, int cq_) : 
    isnew(true),
    min_t(std::numeric_limits<double>::infinity())
  {
    fp = fp_;
    fq = fq_;
    cp = cp_;
    cq = cq_;
    if(p0_ <= p1_)
    {
      p0 = p0_;
      p1 = p1_;
    }
    else
    {
      p0 = p1_;
      p1 = p0_;
    }

    if(q0_ <= q1_)
    {
      q0 = q0_;
      q1 = q1_;
    }
    else
    {
      q0 = q1_;
      q1 = q0_;
    }

    if(p0 > q0)
    {
      std::swap(q0, p0);
      std::swap(q1, p1);
      std::swap(fq, fp);
      std::swap(cq, cp);
    }

    assert(p0 <= p1 && q0 <= q1 && p0 <= q0);
  }

  bool operator<(const EdgeEdgeStencil &other) const
  {
    if(p0 < other.p0)
	return true;
    else if(p0 > other.p0)
	return false;
   
    if(p1 < other.p1)
	return true;
    else if(p1 > other.p1)
	return false;

    if(q0 < other.q0)
	return true;
    else if(q0 > other.q0)
	return false;

    return q1 < other.q1;
  }
};

#endif
