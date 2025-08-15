function [x, xPhys, change] = updateDensityVariables(x, dcdx, dvdx, h, Hs, beta, vfmax)
l1 = 0; l2 = 1e9; move = 0.2;
while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dcdx./dvdx/lmid)))));
    xTilde = conv2(xnew,h,'same')./Hs;
    xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
%     xPhys = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta)));
    if mean(xPhys(:)) >= vfmax
        l1 = lmid;
    else
        l2 = lmid;
    end
end
change = max(abs(xnew(:)-x(:)));
x = xnew;