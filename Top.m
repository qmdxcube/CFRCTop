%%Top(nelx,nely,he,thickness,volfrac,rmin)
function Top(nelx,nely,he,thickness,volfrac,rmin)
%% Material properties
E = 70e9; v = 0.3; G = E/2/(1+v);
D0 = inv([1/E,-v/E,0;-v/E,1/E,0;0,0,1/G]);
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
KE = BasicKe(he/2,he/2,thickness,D0);
%% Initial designs for topology and fiber orientation
x = volfrac*ones(nely,nelx);
xTilde = zeros(nely,nelx);
xPhys = zeros(nely,nelx);
%% FILTER Preparation
[H,Hs] = PreFILTER(x,rmin);
%% Optimization parameters
iter = 0;
change = 1;
loopbeta = 0;
maxiter = 500;
beta = 0;
compliance = zeros(maxiter,1);
volumefrac = zeros(maxiter,1);
%%START ITERATION
while (change>0.01 && iter<maxiter)
    iter = iter+1;
    loopbeta = loopbeta+1;
    %% FILTERING AND PROJECTION/NORMALIZATION
    xTilde(:) = (H*x(:))./Hs;
    xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
    %% FE-ANALYSIS
    [xmp,dxmp] = MaterialInterpolation(xPhys(:)); % Material interpolation
    sK = KE(:)*xmp';
    K = sparse(iK,jK,sK);
    K = (K+K')/2;
    Kf = K(freedofs,freedofs);
    U(freedofs,1) = Kf\F(freedofs,1);
    compliance(iter,1) = F'*U;
    volumefrac(iter,1) = mean(xPhys(:));
    %% SENSITIVITY ANALYSIS
    Ud = U(edofMat);
    ce = sum((Ud*KE).*Ud,2);
    dcdxPhys = -dxmp.*ce;
    %% FILTERING OF SENSITIVITIES
    dx = beta*exp(-beta*xTilde)+exp(-beta);
    dvdxPhys = ones(nelm,1);
    dcdx = reshape(H*(dcdxPhys.*dx(:)./Hs),nely,nelx);
    dvdx = reshape(H*(dvdxPhys.*dx(:)./Hs),nely,nelx);
    %% DESIGN VARIABLE UPDATE
    [x,xPhys,change] = updateDensityVariables(x,dcdx,dvdx,H,Hs,beta,volfrac);
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
%% VISUALIZATION OF THE OPTIMIZED DESIGN
plotHistory(compliance,volumefrac)
figure;
colormap(jet); imagesc(xPhys); caxis([0 1]); axis equal; axis off; drawnow;
save('Top_results.mat','compliance','volumefrac','xPhys');
end