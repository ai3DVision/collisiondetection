function [SV,SF,G,hitG,side] = sweep_surface(VV,F)

  function in = in_bb(G,V,eta)
    in = all(bsxfun(@ge,G,min(V)-eta),2) & all(bsxfun(@le,G,max(V)+eta),2);
  end

  sx = 100;
  fprintf('building grid...\n');
  [G,side,r] = voxel_grid(reshape(permute(VV,[1 3 2]),[],3),sx,'Pad',2);
  eta = min(r);
  hitG = false(size(G,1),1);
  fprintf('dtcd...\n');

  for t = 1:size(VV,3)
    not_hit = find(~hitG);
    not_hit = not_hit(in_bb(G(not_hit,:),VV(:,:,t),eta));
    fprintf('  %d\n',numel(not_hit));
    hitG(not_hit) = signed_distance(G(not_hit,:),VV(:,:,t),F,'SignedDistanceType','pseudonormal')<eta;
  end

  %clf;
  %hold on;
  %s = scatter3(G(:,1),G(:,2),G(:,3),'.');
  %tsh = tsurf(F,VV(:,:,1),'FaceColor',[0.3 0.4 0.8]);
  %hold off;
  %camproj('persp');
  %colormap([0.8 0.8 0.8;0 0 0]);
  %axis equal;

  fprintf('ctcd...\n');
  for t = 2:size(VV,3)
    not_hit = find(~hitG);
    not_hit = not_hit(in_bb(G(not_hit,:),[VV(:,:,t);VV(:,:,t-1)],eta));
    fprintf('  %d\n',size(not_hit,1));
    V0 = [VV(:,:,t-1);G(not_hit,:)];
    V1 = [VV(:,:,t  );G(not_hit,:)];
    % There can't be edge-edge collisions
    VF = ctcd(V0,F,V1,'OuterEta',eta,'InnerEta',eta);
    hitG(not_hit(VF(VF(:,1)>size(VV,1),1)-size(VV,1))) = true;
    %s.CData = 1*hitG;
    %tsh.Vertices = V1;
    %drawnow;
  end
  resh = @(X) reshape(X,side([2 1 3]));
  fprintf('isosurface...\n');
  [SF,SV] = isosurface( resh(G(:,1)), resh(G(:,2)), resh(G(:,3)), resh(hitG));
end

