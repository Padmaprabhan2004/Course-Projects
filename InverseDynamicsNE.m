function tau= InverseDynamicsNE(theta,dtheta,ddtheta, ...
    g,ftip,Mlist,Glist,Slist)

%Mlist is list of link frames {i} relative to {i-1} at home position.
%Slist is the screw axes Ai in base frame. So its screw wrt ground. But we
%need screw axes wrt {i}.
n=size(theta,1);
Mi=eye(4);
Ai=zeros(6,n);
AdTi=zeros(6,6,n+1);
Vi=zeros(6,n+1);
Vdi=zeros(6,n+1);
%define V0 as (0,-g)
Vdi(4:6,1)=-g;
Fi=ftip;
tau=zeros(n,1);
for i=1:n
    Mi_zero=Mi*Mlist(:,:,i);
    Ai(:,i)=Adjoint(TransInv(Mi_zero))*Slist(:,i);
    Ti_iprev=MatrixExp6(-theta(i)*VecTose3(Ai(:,i)))*TransInv(Mlist(:,:,i));
    AdTi(:,:,i)=Adjoint(Ti_iprev);
    Vi(:,i+1)=AdTi(:,:,i)*Vi(:,i)+Ai(:,i)*dtheta(i);
    Vdi(:,i+1)=AdTi(:,:,i)*Vi(:,i)+Ai(:,i)*ddtheta(i)+ad(Vi(:,i+1))*Ai(:,i)*dtheta(i);
end
for i=n:-1:1
    Fi=AdTi(:,:,i+1)'*Fi+Glist(:,:,i)*Vdi(:,i+1)-ad(Vi(:,i+1))'*(Glist(:,:,i)*Vi(:,i+1));
    tau(i)=Fi'*Ai(:,i);
end
end