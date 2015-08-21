// mex( ...
//   '-v', ...
//   '-largeArrayDims', ...
//   'ctcd.cpp', ...
//   '-I/usr/local/igl/libigl/include', ...
//   '-I../../', ...
//   '-I/usr/local/include/eigen3/', ...
//   '-DMEX', ...
//   '-L../bin/','-lcollisions', ...
//   'CXXFLAGS=$CXXFLAGS -std=c++11');
//

#include <mex.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/C_STR.h>

#include <collisiondetection/src/Mesh.h>
#include <collisiondetection/src/History.h>
#include <collisiondetection/src/KDOPBroadPhase.h>
#include <collisiondetection/src/CTCDNarrowPhase.h>
#include <collisiondetection/src/SeparatingPlaneNarrowPhase.h>
#include <collisiondetection/src/stencils.h>
#include <set>
#include <iostream>

#ifdef MEX
#  include "mex.h"
#endif

#include <string>

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  // This is useful for debugging whether Matlab is caching the mex binary
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);

  using namespace std;
  using namespace Eigen;
  using namespace igl;
  using namespace igl::matlab;

  double outer_eta = 1e-7;
  double inner_eta = 1e-8;
  Matrix<double,Dynamic,3,RowMajor> V0,V1;
  Matrix<int,Dynamic,3,RowMajor> F;
  if(nrhs < 3)
  {
    mexErrMsgTxt("nrhs < 3");
  }
  parse_rhs_double(prhs,V0);
  parse_rhs_index(prhs+1,F);
  parse_rhs_double(prhs+2,V1);
  {
    int i = 3;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),"Parameter names should be strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      const auto requires_arg = 
        [](const int i, const int nrhs, const char * name)
      {
        mexErrMsgTxt((i+1)<nrhs,
          C_STR("Parameter '"<<name<<"' requires argument"));
      };
      const auto validate_char = 
        [](const int i, const mxArray * prhs[], const char * name)
      {
        mexErrMsgTxt(mxIsChar(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires char argument"));
      };
      const auto validate_double= 
        [](const int i, const mxArray * prhs[], const char * name)
      {
        mexErrMsgTxt(mxIsDouble(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires Double argument"));
      };
      const auto validate_logical= 
        [](const int i, const mxArray * prhs[], const char * name)
      {
        mexErrMsgTxt(mxIsLogical(prhs[i]),
          C_STR("Parameter '"<<name<<"' requires Logical argument"));
      };
      const auto validate_scalar = 
        [](const int i, const mxArray * prhs[], const char * name)
      {
        mexErrMsgTxt(mxGetN(prhs[i])==1 && mxGetM(prhs[i])==1,
          C_STR("Parameter '"<<name<<"' requires scalar argument"));
      };
      if(strcmp("OuterEta",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_scalar(i,prhs,name);
        validate_double(i,prhs,name);
        outer_eta = mxGetScalar(prhs[i]);
      }else if(strcmp("InnerEta",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_scalar(i,prhs,name);
        validate_double(i,prhs,name);
        inner_eta = mxGetScalar(prhs[i]);
      }else
      {
        mexErrMsgTxt(false,"Unknown parameter");
      }
      i++;
    }
  }

  mexErrMsgTxt(V0.cols()==3,"V0 must be #V0 by 3");
  mexErrMsgTxt(F.cols()==3,"F must be #F by 3");
  mexErrMsgTxt(V0.rows() == V1.rows() && V1.cols()==3,"V1 must be #V0 by 3");

  // Build data structures: mesh and history (trajectories)
  Mesh m;
  m.faces = Eigen::Map<Eigen::Matrix3Xi>(F.data(),3,F.rows());
  // These are needed but ignored
  m.vertices = Eigen::Map<Eigen::VectorXd>(V0.data(),V0.size(),1);
  History h( Eigen::Map<Eigen::VectorXd>(V0.data(),V0.size(),1));
  h.finishHistory( Eigen::Map<Eigen::VectorXd>(V1.data(),V1.size(),1));
  // Broadphase to collect candidate vertex-face and edge-edge collisions
  std::set<VertexFaceStencil> can_vfs;
  std::set<EdgeEdgeStencil> can_ees;
  std::set<std::pair<VertexFaceStencil,double> > can_vfs_eta;
  std::set<std::pair<EdgeEdgeStencil,double> > can_ees_eta;
  // No fixed vertices
  const std::set<int> fixedVerts;
  KDOPBroadPhase().findCollisionCandidates(h,m,outer_eta,can_vfs,can_ees,fixedVerts);
  // Narrowphase 
  for(const auto & vf : can_vfs)
  {
    can_vfs_eta.emplace(vf,inner_eta);
  }
  for(const auto & ee : can_ees)
  {
    can_ees_eta.emplace(ee,inner_eta);
  }
  std::set<VertexFaceStencil> vfs;
  std::set<EdgeEdgeStencil> ees;
  CTCDNarrowPhase().findCollisions(h,can_vfs_eta,can_ees_eta,vfs,ees);
  //std::cout<<"CTCDNarrowPhase:"<<std::endl;
  //std::cout<<"  "<<vfs.size()<<" vf collisions"<<std::endl;
  //std::cout<<"  "<<ees.size()<<" ee collisions"<<std::endl;
  Eigen::MatrixX2i VF,EE;
  VF.resize(vfs.size(),2);
  {
    int i = 0;
    for(const auto vf : vfs)
    {
      VF(i,0) = vf.p;
      VF(i,1) = vf.fq;
      i++;
    }
  }
  EE.resize(ees.size(),2);
  {
    int i = 0;
    for(const auto ee : ees)
    {
      EE(i,0) = ee.fp+ee.cp*F.rows();
      EE(i,1) = ee.fq+ee.cq*F.rows();
      i++;
    }
  }

  // Prepare left-hand side
  switch(nlhs)
  {
    case 2:
    {
      // Treat indices as reals
      plhs[1] = mxCreateDoubleMatrix(EE.rows(),EE.cols(), mxREAL);
      double * EEp = mxGetPr(plhs[1]);
      MatrixX2d EEd = (EE.cast<double>().array()+1).matrix();
      copy(&EEd.data()[0],&EEd.data()[0]+EEd.size(),EEp);
      // Fallthrough
    }
    case 1:
    {
      // Treat indices as reals
      plhs[0] = mxCreateDoubleMatrix(VF.rows(),VF.cols(), mxREAL);
      double * VFp = mxGetPr(plhs[0]);
      MatrixX2d VFd = (VF.cast<double>().array()+1).matrix();
      copy(&VFd.data()[0],&VFd.data()[0]+VFd.size(),VFp);
      // Fallthrough
    }
    default:break;
  }

  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
}

