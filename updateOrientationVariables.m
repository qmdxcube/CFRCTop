%% -------------------------------------- MMA UPDATE SCHEME (UNCONSTRAINED)
function [zNew,zold1,zold2,L,U] = updateOrientationVariables(dfdz,z,zold1,zold2,L,U,Iter) 
zMin = -1.0;zMax = 1.0;move = 0.2;
AsymInit = 0.2; AsymInc = 1.2; AsymDecr = 0.7;
xmin = max(zMin,z-move); xmax = min(zMax,z+move);
% Compute asymptotes L and U:
if Iter<=2
  L = z - AsymInit*(xmax-xmin);  U = z + AsymInit*(xmax-xmin);    
else
  sgn = (z-zold1).*(zold1-zold2);
  s = ones(size(z)); s(sgn>0) = AsymInc; s(sgn<0) = AsymDecr;
  L = z - s.*(zold1 - L); U = z + s.*(U - zold1);
end
% Compute bounds alpha and beta
alpha = 0.9*L + 0.1*z;   beta = 0.9*U + 0.1*z;
alpha = max(xmin,alpha); beta = min(xmax,beta);
% Solve unconstrained subproblem
feps = 0.000001; 
p = (U-z).^2.*(max(dfdz,0)+0.001*abs(dfdz)+feps./(U-L)); 
q = (z-L).^2.*(-min(dfdz,0)+0.001*abs(dfdz)+feps./(U-L));
zCnd = (L.*p - U.*q + (U - L).*sqrt(p.*q))./(p - q);
zNew = max(alpha,min(beta,zCnd));
zold2 = zold1; zold1=z;