function plotDesign(xPhys,pxPhys,pyPhys)
figure
colormap(gray);imagesc(1-xPhys);axis equal;axis tight;axis off;pause(1e-6);
hold on
xMax =2*max(xPhys(:));
quiver(-pxPhys.*xPhys/xMax,pyPhys.*xPhys/xMax,0,'y.','LineWidth',1);
quiver(pxPhys.*xPhys/xMax,-pyPhys.*xPhys/xMax,0,'y.','LineWidth',1);
end

