function CFRCTopU(nelx,nely,he,thickness,volfrac,rmin,rmin_p)
%% Material properties
E1 = 134.3e9; E2 = 8.5e9; v21 = 0.34; G12 = 5.8e9;
D0 = inv([1/E1,-v21/E1,0;-v21/E1,1/E2,0;0,0,1/G12]);
%% Define support and load
nelm = nelx*nely; % Number of elements
ndof = 2*(nelx+1)*(nely+1); % Number of dofs
alldofs = 1:ndof; % All dofs
fixeddofs = union(1:2:2*(nely+1),2*(nelx+1)*(nely+1)); % fixed dofs
freedofs = setdiff(alldofs,fixeddofs);% free dofs
F = sparse(2,1,-1000,ndof,1); % Spatial distribution of the load
U = zeros(ndof,1);
%% Finite element analysis prepration
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelm,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelm,1);
iK = kron(edofMat,ones(8,1))';
jK = kron(edofMat,ones(1,8))';
KET = zeros(8,8,6);% Template Stiffness Matrices (TSMs) in matrix form
KETc = zeros(64,6);% Template Stiffness Matrices (TSMs) in vector form
for i=1:3
    for j=1:i
        n = i*(i-1)/2+j;
        DT = zeros(3,3);
        DT(i,j) = 1;
        DT(j,i) = 1;
        KEt = BasicKe(he/2,he/2,thickness,DT);
        KET(:,:,n) = KEt;
        KETc(:,n) = KEt(:);
    end
end
%% Initial designs for topology and fiber orientation
x = volfrac*ones(nely,nelx);
px = ones(nely,nelx);
py = zeros(nely,nelx);
xPhys = zeros(nely,nelx);
pxPhys = zeros(nely,nelx);
pyPhys = zeros(nely,nelx);
%% FILTER Preparation
[dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
h = max(0,rmin-sqrt(dx.^2+dy.^2));
Hs = conv2(ones(nely,nelx),h,'same');
[dpy,dpx] = meshgrid(-ceil(rmin_p)+1:ceil(rmin_p)-1,-ceil(rmin_p)+1:ceil(rmin_p)-1);
hp = max(0,rmin-sqrt(dpx.^2+dpy.^2));
Hps = conv2(ones(nely,nelx),hp,'same');
%% Optimization parameters
iter = 0;
change = 1;
loopbeta = 0;
maxiter = 500;
beta = 0;
p = [px(:);py(:)];
pold1 = p;
pold2 = p;
pL = -ones(2*nelm,1);
pU = ones(2*nelm,1);
compliance = zeros(maxiter,1);
volumefrac = zeros(maxiter,1);
%%START ITERATION
while (change>0.01 && iter<maxiter)
    iter = iter+1;
    loopbeta = loopbeta+1;
    %% FILTERING AND PROJECTION/NORMALIZATION
    xTilde = conv2(x,h,'same')./Hs;
    pxTilde = conv2(px,hp,'same')./Hps;
    pyTilde = conv2(py,hp,'same')./Hps;
    xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
    pxPhys = pxTilde./sqrt(pxTilde.^2+pyTilde.^2);
    pyPhys = pyTilde./sqrt(pxTilde.^2+pyTilde.^2);
    %% FE-ANALYSIS
    l1 = pxPhys(:); m1 = pyPhys(:);
    [D,Dl,Dm] = ComputeElasticityMatrices(D0,l1,m1); % Compute elasticity matrices
    [xmp,dxmp] = MaterialInterpolation(xPhys); % Material interpolation
    sK = zeros(64,nelm);
    for i = 1:3
        for j = 1:i
            n = i*(i-1)/2+j;
            Dij = D(i,j,:);
            sK = sK+KETc(:,n)*(Dij(:).*xmp(:))';
        end
    end
    K = sparse(iK,jK,sK);
    K = (K+K')/2;
    Kf = K(freedofs,freedofs);
    U(freedofs,1) = Kf\F(freedofs,1);
    compliance(iter,1) = F'*U;
    volumefrac(iter,1) = mean(xPhys(:));
    %% SENSITIVITY ANALYSIS
    ce = zeros(nelm,1);
    cex = zeros(nelm,1);
    cey = zeros(nelm,1);
    Ud = U(edofMat);
    for i = 1:3
        for j = 1:i
            n = i*(i-1)/2+j;
            KE = KET(:,:,n);
            uku = sum((Ud*KE).*Ud,2);
            Dij = D(i,j,:);
            Dlij = Dl(i,j,:);
            Dmij = Dm(i,j,:);
            ce = ce+Dij(:).*uku;
            cex = cex+Dlij(:).*uku;
            cey = cey+Dmij(:).*uku;
        end
    end
    dcdxPhys = -dxmp.*reshape(ce,nely,nelx);
    dcdpxPhys = -xmp.*reshape(cex,nely,nelx);
    dcdpyPhys = -xmp.*reshape(cey,nely,nelx);
    %% FILTERING OF SENSITIVITIES
    dx = beta*exp(-beta*xTilde)+exp(-beta);
    pPhys = zeros(2,1,nelm);
    pPhys(1,1,:) = l1;
    pPhys(2,1,:) = m1;
    dp = pagemtimes(pPhys,'none',pPhys,'transpose');
    I = zeros(size(dp));
    I(1,1,:) = 1.0;
    I(2,2,:) = 1.0;
    dcdpPhys = zeros(2,1,nelm);
    dcdpPhys(1,1,:) = dcdpxPhys(:);
    dcdpPhys(2,1,:) = dcdpyPhys(:);
    dcdpTiles = pagemtimes(I-dp,dcdpPhys);
    dcdpxTiles_t = dcdpTiles(1,1,:);
    dcdpyTiles_t = dcdpTiles(2,1,:);
    Lp = sqrt(pxTilde(:).^2+pyTilde(:).^2);
    dcdpxTiles = dcdpxTiles_t(:)./Lp;
    dcdpyTiles = dcdpyTiles_t(:)./Lp;
    dvdxPhys = ones(nely,nelx);
    dcdx = conv2(dcdxPhys.*dx./Hs,h,'same');
    dvdx = conv2(dvdxPhys.*dx./Hs,h,'same');  
    %% DESIGN VARIABLE UPDATE
    [x,xPhys,change] = updateDensityVariablesU(x,dcdx,dvdx,h,Hs,beta,volfrac);
    dfdpx = conv2(reshape(dcdpxTiles,nely,nelx)./Hps,hp,'same');
    dfdpy = conv2(reshape(dcdpyTiles,nely,nelx)./Hps,hp,'same');
    dfdp = [dfdpx(:);dfdpy(:)];
    [p,pold1,pold2,pL,pU] = updateOrientationVariables(dfdp,p,pold1,pold2,pL,pU,iter);
    px = reshape(p(1:nelm,1),nely,nelx);
    py = reshape(p(nelm+(1:nelm),1),nely,nelx);
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',iter,compliance(iter,1),volumefrac(iter,1),change);
    %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if  beta < 256 && (loopbeta >= 50 || change < 0.01)
        beta = max(2*beta,beta+1);
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
end
compliance = compliance(1:iter,1);
volumefrac = volumefrac(1:iter,1);
%% VALIDATION OF THE OPTIMALITY OF FIBER ORIENTATION
verifyOptimality(nelx,nely,he,U,D,xPhys,pxPhys,pyPhys)
%% VISUALIZATION OF THE OPTIMIZED DESIGN
plotHistory(compliance,volumefrac)
plotDesign(xPhys,pxPhys,pyPhys)
save('CFRCTop_results.mat','compliance','volumefrac','xPhys','pxPhys','pyPhys');
end