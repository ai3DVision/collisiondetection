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
// clang++ -o ctcd -g -std=c++11 ctcd.cpp -I/usr/local/include/eigen3/ -I/usr/local/igl/libigl/include -I../../ -L../bin/ -lcollisions

#ifdef MEX
#  include <mex.h>
#  undef assert
#  define assert( isOK ) ( (isOK) ? (void)0 : (void) mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )
#  include <igl/matlab/MexStream.h>
#  include <igl/matlab/mexErrMsgTxt.h>
#  include <igl/matlab/parse_rhs.h>
#else
#  include <igl/read_triangle_mesh.h>
#  include <igl/readDMAT.h>
#  include <igl/get_seconds.h>
#endif
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

#ifdef MEX
void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
  // This is useful for debugging whether Matlab is caching the mex binary
  //mexPrintf("%s %s\n",__TIME__,__DATE__);
  igl::matlab::MexStream mout;
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  using namespace igl::matlab;
#else
int main(int argc, char * argv[])
{
#endif

  using namespace std;
  using namespace Eigen;
  using namespace igl;

  double outer_eta = 1e-7;
  double inner_eta = 1e-8;
  Matrix<double,Dynamic,3,RowMajor> V0,V1;
  VectorXi fixedVerts,C;
  Matrix<int,Dynamic,3,RowMajor> F;
#ifdef MEX
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
      }else if(strcmp("Fixed",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_double(i,prhs,name);
        parse_rhs_index(prhs+i,fixedVerts);
      }else if(strcmp("Components",name) == 0)
      {
        requires_arg(i,nrhs,name);
        i++;
        validate_double(i,prhs,name);
        parse_rhs_index(prhs+i,C);
        mexErrMsgTxt(C.size() == 0 || C.size() == V0.rows(),
          "C should be empty or #C should equal #V");
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
#else
  switch(argc)
  {
    default:
      cerr<<R"(Usage:
./ctcd mesh-at-0.obj mesh-at-1.obj components.dmat
)";
      exit(EXIT_FAILURE);
    case 4:
      readDMAT(argv[3],C);
      // fall through
    case 3:
      read_triangle_mesh(argv[1],V1,F);
      read_triangle_mesh(argv[2],V0,F);
      break;
  }
  assert(C.size() == 0 || C.rows() == V0.rows());
  assert(C.size() == 0 || C.minCoeff() == 0);
#endif

#ifndef MEX
#define MAX_RUNS 4
  const double t_before = get_seconds();
  for(int r = 0;r<MAX_RUNS;r++)
  {
#endif
  //cout<<"Building data structures..."<<endl;
  //cout<<"  "<<V0.rows()<<" vertices..."<<endl;
  //cout<<"  "<<F.rows()<<" faces..."<<endl;
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
  const std::set<int> sfixedVerts(fixedVerts.data(),fixedVerts.data()+fixedVerts.size());
  //std::cout<<"KDOPBroadPhase:"<<std::endl;
  KDOPBroadPhase().findCollisionCandidates(
    h,m,outer_eta,sfixedVerts,C,can_vfs,can_ees);
  //std::cout<<"  "<<can_vfs.size()<<" candidate vf collisions"<<std::endl;
  //std::cout<<"  "<<can_ees.size()<<" candidate ee collisions"<<std::endl;
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
  //std::cout<<"CTCDNarrowPhase:"<<std::endl;

  const bool at_most_one = false;
  Matrix<bool,Dynamic,1> V_flag = Matrix<bool,Dynamic,1>::Zero(V0.rows(),1);
  Matrix<bool,Dynamic,1> F_flag = Matrix<bool,Dynamic,1>::Zero(F.rows(),1);
  //Matrix<bool,Dynamic,1> HE_flag = Matrix<bool,Dynamic,1>::Zero(3*F.rows(),1);
  if(at_most_one)
  {
    CTCDNarrowPhase narrow;
    for(set<pair<VertexFaceStencil, double> >::const_iterator it = can_vfs_eta.begin(); it != can_vfs_eta.end(); ++it)
    {
      VertexFaceStencil vf =  it->first;
      if(!V_flag(vf.p) && !F_flag(vf.fq) && narrow.checkVFS(h, vf, it->second))
      {
        V_flag(vf.p) = true;
        F_flag(vf.fq) = true;
        vfs.insert(vf);
      }
    }
    for(set<pair<EdgeEdgeStencil, double> >::const_iterator it = can_ees_eta.begin(); it != can_ees_eta.end(); ++it)
    {
      EdgeEdgeStencil ee = it->first;
      const int hep = ee.fp+ee.cp*F.rows();
      const int heq = ee.fq+ee.cq*F.rows();
      //if(!HE_flag(hep) && !HE_flag(heq) && narrow.checkEES(h, ee, it->second))
      //{
      //  HE_flag(hep) = true;
      //  HE_flag(heq) = true;
      //  ees.insert(ee);
      //}
      if(!F_flag(ee.fp) && !F_flag(ee.fq) && narrow.checkEES(h, ee, it->second))
      {
        F_flag(ee.fp) = true;
        F_flag(ee.fq) = true;
        ees.insert(ee);
      }
    }
  }else
  {
    CTCDNarrowPhase().findCollisions(h,can_vfs_eta,can_ees_eta,vfs,ees);
  }

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
  //cout<<"VF "<<VF.rows()<<endl;
  //cout<<"EE "<<EE.rows()<<endl;
#ifndef MEX
  }
  cout<<(get_seconds()-t_before)/(double)MAX_RUNS<<endl;
#endif

#ifdef MEX
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
#else
}
#endif

