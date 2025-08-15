function [D,dDl,dDm]=ComputeElasticityMatrices(D0,l1,m1)
Tc=zeros(3,3,numel(l1));
dTcl=zeros(3,3,numel(l1));
dTcm=zeros(3,3,numel(l1));
Tc(1,1,:)=l1.^2;Tc(1,2,:)=m1.^2;Tc(1,3,:)=-2*l1.*m1;
Tc(2,1,:)=m1.^2;Tc(2,2,:)=l1.^2;Tc(2,3,:)=2*m1.*l1;
Tc(3,1,:)=l1.*m1;Tc(3,2,:)=-l1.*m1;Tc(3,3,:)=l1.^2-m1.^2;
dTcl(1,1,:)=2*l1;dTcl(1,2,:)=0;dTcl(1,3,:)=-2*m1;
dTcl(2,1,:)=0;dTcl(2,2,:)=2*l1;dTcl(2,3,:)=2*m1;
dTcl(3,1,:)=m1;dTcl(3,2,:)=-m1;dTcl(3,3,:)=2*l1;
dTcm(1,1,:)=0;dTcm(1,2,:)=2*m1;dTcm(1,3,:)=-2*l1;
dTcm(2,1,:)=2*m1;dTcm(2,2,:)=0;dTcm(2,3,:)=2*l1;
dTcm(3,1,:)=l1;dTcm(3,2,:)=-l1;dTcm(3,3,:)=-2*m1;
TcD0 = pagemtimes(Tc, D0);
D = pagemtimes(TcD0, pagetranspose(Tc));
dDl = pagemtimes(pagemtimes(dTcl, D0), pagetranspose(Tc)) + ...
      pagemtimes(TcD0, pagetranspose(dTcl));
dDm = pagemtimes(pagemtimes(dTcm, D0), pagetranspose(Tc)) + ...
      pagemtimes(TcD0, pagetranspose(dTcm));