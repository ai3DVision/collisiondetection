% CTCD Detect continious time collisions occurring for a mesh with faces F and
% vertices moving from V0 to V1.
%
% Inputs:
%   V0  #V0 by 3 list of vertex positions at t=0
%   F  #F by 3 list of triangle indices into V0 (and V1)
%   V1  #V0 by 3 list of vertex positions at t=1
%   Optional:
%     'InnerEta' followed by inner eta value {1e-7}
%     'OuterEta' followed by outer eta value {1e-8}
%     'Fixed' followed by a list of "fixed" vertices (collisions with these
%        vertices are ignored.
%     'Components' followed by a #V list of component ids intersections within
%       the same component are ignored, [] means no components (all
%       intersections found).
% Outputs:
%   VF  #VF by 2 list of pairs of vertex indices and face indices which undergo
%     interference.
%   EE  #EE by 2 list of pairs of half-edges (indices into 1:#F*3) which
%     undergo interference.
%  
