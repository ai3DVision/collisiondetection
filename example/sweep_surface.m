function [SV,SF,G,hitG,side] = sweep_surface(VV,F)

  function in = in_bb(G,V,eta)
    in = all(bsxfun(@ge,G,min(V)-eta),2) & all(bsxfun(@le,G,max(V)+eta),2);
  end

  sx = 40;
  fprintf('building grid...\n');
  [G,side,r] = voxel_grid(reshape(permute(VV,[1 3 2]),[],3),sx,'Pad',2);
  eta = min(r);
  hitG = false(size(G,1),1);

  fprintf('dtcd...\n');
  for t = 1:size(VV,3)-1
    nss = 0;
    substeps = linspace(0,1,2+nss);
    if t ~= 1
      substeps = substeps(2:end);
    end
    % linear sub-stepping
    for s = substeps
      not_hit = find(~hitG);
      VVs = VV(:,:,t) + s*(VV(:,:,t+1) - VV(:,:,t));
      not_hit = not_hit(in_bb(G(not_hit,:),  VVs,eta));
      new_hit = signed_distance(G(not_hit,:),VVs,F,'SignedDistanceType','pseudonormal')<eta;
      fprintf('  %d -> %d\n',numel(not_hit),sum(new_hit));
      hitG(not_hit) = new_hit;
    end
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
    V0 = [VV(:,:,t-1);G(not_hit,:)];
    V1 = [VV(:,:,t  );G(not_hit,:)];
    % There can't be edge-edge collisions
    C = [ones(size(VV,1),1);2*ones(size(not_hit,1),1)];
    VF = ctcd(V0,F,V1,'OuterEta',eta,'InnerEta',eta,'Components',C);
    new_hit = VF(VF(:,1)>size(VV,1),1)-size(VV,1);
    hitG(not_hit(new_hit)) = true;
    fprintf('  %d -> %d\n',size(not_hit,1),size(new_hit,1));
    %s.CData = 1*hitG;
    %tsh.Vertices = V1;
    %drawnow;
  end
  resh = @(X) reshape(X,side([2 1 3]));

  fprintf('filter...\n');
  hitG = resh(hitG);
  dilate_side = [3 3 3];
  hitG = imdilate(hitG,ones(dilate_side));
  smooth_k = 3;
  hitG = smooth3(hitG,'box',smooth_k);

  fprintf('isosurface...\n');
  level = 0.5;
  [SF,SV] = isosurface(resh(G(:,1)), resh(G(:,2)), resh(G(:,3)),hitG,level);
  % isosurface uses wrong orientation
  SF = fliplr(SF);
  decimate_factor = 0.1;
  if decimate_factor < 1
    fprintf('decimate...\n');
    [SV,SF] = decimate_cgal(SV,SF,decimate_factor,'Adaptive',false);
  end
end

