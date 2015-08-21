% CTCD Detect continious time collisions occurring for a mesh with faces F and
% vertices moving from V0 to V1.
%
% Inputs:
%   V0  #V0 by 3 list of vertex positions at t=0
%   F  #F by 3 list of triangle indices into V0 (and V1)
%   V1  #V0 by 3 list of vertex positions at t=1
% Outputs:
%   VF  #VF by 2 list of pairs of vertex indices and face indices which undergo
%     interference.
%   EE  #EE by 2 list of pairs of half-edges (indices into 1:#F*3) which
%     undergo interference.
%  
