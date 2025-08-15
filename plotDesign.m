function plotDesign(xPhys,pxPhys,pyPhys)
figure;
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
figure;
s = 0.7;
[nely,nelx]=size(xPhys);
[X,Y] = meshgrid(1:nelx,1:nely);
pxPhys(xPhys<0.5)=0;
pyPhys(xPhys<0.5)=0;
q = quiver(X-s*pxPhys/2,-Y-s*pyPhys/2,pxPhys,pyPhys,s);
q.ShowArrowHead = 'off';axis equal; axis off;grid on;box on
end

