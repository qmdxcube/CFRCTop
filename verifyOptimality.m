function verifyOptimality(nelx,nely,he,U,D,xPhys,pxPhys,pyPhys)
%% Finite element analysis prepration
nelm = nelx*nely; % Number of elements
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelm,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelm,1);
Ud = U(edofMat');
a = he/2; b = he/2;
B = zeros(3,8);
dNx = [-1 1 1 -1]/4/a;
dNy = [-1 -1 1 1]/4/b;
B(1,1:2:8) = dNx; B(2,2:2:8) = dNy;
B(3,1:2:8) = dNy; B(3,2:2:8) = dNx;
strain = B*Ud;
strain_3d = reshape(strain, 3, 1, []);
stress_3d = pagemtimes(D, strain_3d);
stress = squeeze(stress_3d);
correlation = zeros(nelm,1);
l1 = pxPhys(:); m1 = pyPhys(:);
for e=1:nelm
    [pVectors,pStress] = eig([stress(1,e),stress(3,e);stress(3,e),stress(2,e)]);
    [~,idx] = sort(abs(diag(pStress)), 'descend');
    pVectors_sorted = pVectors(:,idx);
    correlation(e,1) = abs(dot(pVectors_sorted(:,1),[l1(e),m1(e)]));
end
correlation = reshape(correlation,nely,nelx);
correlation(xPhys<0.5)=0;
figure;
colormap(jet); imagesc(correlation); caxis([0 1]); axis equal; axis off; drawnow;
cb = colorbar('north');
axPos = get(gca, 'Position');
cbPos = cb.Position;
cbPos(2) = axPos(2) + axPos(4) - 0.16;
set(cb, 'Position', cbPos);
exportgraphics(gcf, 'MBB_correlation.png', 'Resolution', 3000);
end