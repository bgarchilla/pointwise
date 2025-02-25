function [A,M,Dx,Dy,Mp,Sp,GD11,GD12,GD22,J,Se,Me,Xi,Eta]=makronvq(T,z,epsilon,Ic,phi,ccphi)
% FOR QUADRATICS WITH SKEW-SYMMETRIC NONLINEAR TERM
% Matrices coef. constantes para -epsilon*\Delta +b\cdot\nabla
% sin ciclo;
%  INPUT
%    T0  nt1 x nt  array with the indexes of the vertex of every
%              triangle. T(i,j), i=1,2,3, j=1,...,nt is the index of
%              the i-th  vertex of triangle number j. Here nt=2*N*N
%              nt1 >= ne, ne=3 for linears, 6 for quadratics, etc.
%    z  nn x 2 array with the xy coordinatex of every vertex in the
%               triangulation. Entries z(n,1) and z(n,2), are, respectively,
%               the x and y coordinates of node number n.
%               Here nn=(N+1)*(N+1)
%    epsilon  real variable with the diffusion parameter
%    Ic, (optinal) indexes of triangles whos map is not the standar affine
%               map, but that given by phi
%    phi  function of the form [x,yxix,xiy,etx,ety,J]=phi(xi,eta,z1,z2,c);
%               returning the transformed of points of coordinates xi
%               and eta in the standard triangle to x and y, toguether
%               with the corresponding derivatives and the Jacobian
%               of the transformation.
%    ccphi (optional) extra argument for phi
% OUTPUT
%    A  nn x nn matrix of the bilinear form associated with
%               with Neumann bc.
%    M  nn x nn mass matrix
%    Dx,Dy,     n1 x nn matrices with the bilinear forms
%                (\psi_i, partial_x\phi_j), (\psi_i, partial_y\phi_j)
%               i=1,2,...,n1, j=1,2,...,nn,
%               (divergence of velocity)
%    GD11,GD12,GD22 nn x nn matrices for the grad div term, i.e., the
%               matrices of bilinear forms (phi_x,phi_x), (phi_x,phi_y)
%               and (phi_y,phy_y)
%    J  nt x 1  vector with the Jacobians of every element.
%    Se ne*nt x ne element matrices of a(u,v)= (grad(u), grad(v))
%    Me ne*nt x ne element mass matrices
%    Mp n1xn1   with the mass matrix for linears.
%    Xi nt x 2*nq element matrix with the values on the nq quadrature nodes
%                 on the reference triangle of the partial derivatives with
%                 respect to x (first nq columns) and y (last nq columns)
%                 of the (horizontal) coordinate xi of the reference triangle
%    Eta          idem but with ehe vertical corrdinate eta of the reference
%                 triangle
%
%
%   Writen by Bosco Garcia-Archilla (last modified: January 2025).
%
%   This code comes with no guarantee or warranty of any kind.
%
%   If you use this code, please cite
%
%   B. Garcia-Archilla and J. Novo,
%     Pointwise error bounds in POD methods without difference quotients,
%     Journal of Scientific Computing (to appear)
%     Paper accepted for publication on February 17, 2025.
%     Please check volume, pages and year of publication with Journal for proper citation.


%nodos de cuadratura
c1=[1 1]/3;
c2=(6+sqrt(15))*[1 1]/21;
c3=[(9-2*sqrt(15)),(6+sqrt(15))]/21;
c4=[c3(2) c3(1)];
c5=(6-sqrt(15))*[1 1]/21;
c6=[(9+2*sqrt(15)),(6-sqrt(15))]/21;
c7=[ c6(2) c6(1)];
c=[c1;c2;c3;c4;c5;c6;c7]';
cx=c(1,:);
cy=(c(2,:));

% pesos en los nodos de cuadratura
w1=0.1125;
w2=(155+sqrt(15))/2400;
w3=w2;
w4=w2;
w5=(155-sqrt(15))/2400;
w6=w5;
w7=w5;

w=[w1 w2 w3 w4 w5 w6 w7]';
% the basis functions on the quadrature nodes
Nc=[(1-cx-cy).*(1-2*cx-2*cy); cx.*(2*cx-1); cy.*(2*cy-1); ...
         4*cx.*(1-cx-cy); 4*cx.*cy; 4*cy.*(1-cx-cy)];
Ncxi=[-(1-2*cx-2*cy)-2*(1-cx-cy); 4*cx-1; zeros(size(cx));...
          4*(1-cx-cy)-4*cx; 4*cy; -4*cy];
Ncet=[-(1-2*cx-2*cy)-2*(1-cx-cy); zeros(size(cx)); 4*cy-1;...
          -4*cx; 4*cx; 4*(1-cx-cy)-4*cy];
      ne=6; ne0=3; % for cuadratrics in velocity and linears in pressure
% NOTA: Ncxi y Ncet comprobadas el 30 julio de 2012.

%Ic=find(T(ne+1,:)); % Ic triangles with side 2-3 curved
%Is=find(~T(ne+1,:)); % Is triangles with all sides straight
   nt=size(T,2); nn=size(z,1);
   x1=z(:,1);y1=z(:,2);x=x1(T(1:ne,:)');y=y1(T(1:ne,:)');
   clear x1 y1;% b=reshape(b,2,1);
   

if nargin>3;
    Is=ones(size(T,2),1); Is(Ic)=0; Is=find(Is);
    Iec=reshape(T(1:ne,Ic),length(Ic)*ne,1);
    Kec=kron((Ic-1)*ne,ones(ne,1))+kron(ones(size(Ic)),[1:ne]');
    [xac,yac,xixc,xiyc,etxc,etyc,Jc]=phi(cx,cy,x(Ic,1:ne),y(Ic,1:ne),ccphi);
    Jcsigned=Jc; Jc=abs(Jc);
else
    Is=[1:size(T,2)]'; Ic=[];
    Iec=[];
end
Ies=reshape(T(1:ne,Is),length(Is)*ne,1);
Kes=kron((Is-1)*ne,ones(ne,1))+kron(ones(size(Is)),[1:ne]');


    % Partial derivatives of x and y wrt xi eta
   xxi=x(:,2)-x(:,1); xet=x(:,3)-x(:,1);
   yxi=y(:,2)-y(:,1); yet=y(:,3)-y(:,1);

% The Jacobians
   J=(x(:,2)-x(:,1)).*(y(:,3)-y(:,1))-(x(:,3)-x(:,1)).*(y(:,2)-y(:,1));
   
   % Parciales de xi y eta con respecto a x e y
   xix=(y(:,3)-y(:,1))./J; xiy=(x(:,1)-x(:,3))./J;
   etx=(y(:,1)-y(:,2))./J; ety=(x(:,2)-x(:,1))./J;
%   J=J(Is); xix=xis(Is); xiy=xiy(Is); etx=etx(Is); ety=ety(Is);

% THE ELEMENTARY MATRICES
   % parciales con respecto a xi y eta de funciones base
   u=[-1 1 0]'; v=[-1 0 1]'; e=ones(ne,1);
   S1=[
3 1 0 -4 0 0
1 3 0 -4 0 0
0 0 0 0 0 0
-4 -4 0 8 0 0
0 0 0 0 8 -8
0 0 0 0 -8 8]/6;
S2 =[
6 1 1 -4 0 -4
1 0 -1 -4 4 0
1 -1 0 0 4 -4
-4 -4 0 8 -8 8
0 4 4 -8 8 -8
-4 0 -4 8 -8 8]/6;
S3 =[
3 0 1 0 0 -4
0 0 0 0 0 0
1 0 3 0 0 -4
0 0 0 8 -8 0
0 0 0 -8 8 0
-4 0 -4 0 0 8]/6;
S12 =[
3 0 1 0 0 -4
1 0 -1 -4 4 0
0 0 0 0 0 0
-4 0 0 4 -4 4
0 0 4 -4 4 -4
0 0 -4 4 -4 4]/6;

S4=[ ...
 6 -1 -1  0 -4  0; ...
-1  6 -1  0  0 -4; ...
-1 -1  6 -4  0  0; ...
 0  0 -4 32 16 16; ...
-4  0  0 16 32 16; ...
 0 -4  0 16 16 32]/360;

% la difusi'on
   Se=zeros(ne*nt,ne);
   Se=(kron(abs(J).*(xix.^2+xiy.^2),S1)+...
       kron(abs(J).*(etx.^2+ety.^2),S3)+...
       kron(abs(J).*(xix.*etx+xiy.*ety),S2));
   if nargin>3
       Se(Kec,:)=(kron(ones(size(Ic)),Ncxi).*kron((xixc.^2+xiyc.^2).*Jc,e))*(diag(w)*(Ncxi'))...
           + (kron(ones(size(Ic)),Ncet).*kron((etxc.^2+etyc.^2).*Jc,e))*(diag(w)*(Ncet'))...
           + (kron(ones(size(Ic)),Ncxi).*kron((xixc.*etxc+xiyc.*etyc).*Jc,e))*(diag(w)*(Ncet'))...
           + (kron(ones(size(Ic)),Ncet).*kron((xixc.*etxc+xiyc.*etyc).*Jc,e))*(diag(w)*(Ncxi'));
   end
   % Transformed of quadrature nodes at 
   % xc=x(:,1)*ones(size(cx))+xxi*cx+xet*cy; yc=y(:,1)*ones(size(cx))+yxi*cx+yet*cy;
   
    Ae=epsilon*Se;
    bC=zeros( nt , 2*length(w));

   % la masa
   Me=kron(abs(J),S4);
   if nargin>3
       Me(Kec,:)=(kron(ones(size(Ic)),Nc).*kron(Jc,e))*(diag(w)*Nc');
   end
   Me=reshape(Me,(ne^2)*nt,1);
   
   % The indexes
     I1=reshape(T(1:ne,:),ne*nt,1);
     I=I1*e'; J1=kron(T(1:ne,:)',e);
     I=reshape(I,ne*ne*nt,1); J1=reshape(J1,ne*ne*nt,1);
   
%  THE FEM MATRICES
A=sparse(I,J1,Ae,nn,nn); M=sparse(I,J1,Me,nn,nn);

Me=reshape(Me,ne*nt,ne);
Xi=[xix,xiy];
Eta=[etx,ety];
% PRUEBA
%D11=kron(J.*bxi,(e*u')/6)+kron(J.*bet,(e*v')/6); D11=reshape(D11,9*nt,1);
%        kron(J.*bxi,(e*u')/6)+kron(J.*bet,(e*v')/6);
%S=reshape(Se,9*nt,1);
%D11=sparse(I,J1,D11,nn,nn); S=sparse(I,J1,S,nn,nn);
%save d11 D11 S
% end PRUEBA

% El grad-div
   
   GD11e=(kron(abs(J).*(xix.^2),S1)+...
       kron(abs(J).*(etx.^2),S3)+...
       kron(abs(J).*(xix.*etx),S2));
%        Se(Kec,:)=(kron(ones(size(Ic)),Ncxi).*kron((xixc.^2+xiyc.^2).*Jc,e))*(diag(w)*(Ncxi'))...
%            + (kron(ones(size(Ic)),Ncet).*kron((etxc.^2+etyc.^2).*Jc,e))*(diag(w)*(Ncet'))...
%            + (kron(ones(size(Ic)),Ncxi).*kron((xixc.*etxc+xiyc.*etyc).*Jc,e))*(diag(w)*(Ncet'))...
%            + (kron(ones(size(Ic)),Ncet).*kron((xixc.*etxc+xiyc.*etyc).*Jc,e))*(diag(w)*(Ncxi'));
   if nargin>3
       GD11e(Kec,:)=(kron(ones(size(Ic)),Ncxi).*kron((xixc.^2).*Jc,e))*(diag(w)*(Ncxi'))...
           + (kron(ones(size(Ic)),Ncet).*kron((etxc.^2).*Jc,e))*(diag(w)*(Ncet'))...
           + (kron(ones(size(Ic)),Ncxi).*kron((xixc.*etxc).*Jc,e))*(diag(w)*(Ncet'))...
           + (kron(ones(size(Ic)),Ncet).*kron((xixc.*etxc).*Jc,e))*(diag(w)*(Ncxi'));
   end
   GD11=sparse(I,J1,reshape(GD11e,numel(GD11e),1),nn,nn);

   GD22e=(kron(abs(J).*(xiy.^2),S1)+...
       kron(abs(J).*(ety.^2),S3)+...
       kron(abs(J).*(xiy.*ety),S2));
   if nargin>3
       GD22e(Kec,:)=(kron(ones(size(Ic)),Ncxi).*kron((xiyc.^2).*Jc,e))*(diag(w)*(Ncxi'))...
           + (kron(ones(size(Ic)),Ncet).*kron((etyc.^2).*Jc,e))*(diag(w)*(Ncet'))...
           + (kron(ones(size(Ic)),Ncxi).*kron((xiyc.*etyc).*Jc,e))*(diag(w)*(Ncet'))...
           + (kron(ones(size(Ic)),Ncet).*kron((xiyc.*etyc).*Jc,e))*(diag(w)*(Ncxi'));
   end
   GD22=sparse(I,J1,reshape(GD22e,numel(GD11e),1),nn,nn);
   
   GD12e=(kron(abs(J).*(xix.*xiy),S1)+...
       kron(abs(J).*(etx.*ety),S3)+...
       kron(abs(J).*(xix.*ety),S12)+...
       kron(abs(J).*(xiy.*etx),S12'));
   if nargin>3
       GD12e(Kec,:)=(kron(ones(size(Ic)),Ncxi).*kron((xixc.*xiyc).*Jc,e))*(diag(w)*(Ncxi'))...
           + (kron(ones(size(Ic)),Ncet).*kron((etxc.*etyc).*Jc,e))*(diag(w)*(Ncet'))...
           + (kron(ones(size(Ic)),Ncxi).*kron((xixc.*etyc).*Jc,e))*(diag(w)*(Ncet'))...
           + (kron(ones(size(Ic)),Ncet).*kron((xiyc.*etxc).*Jc,e))*(diag(w)*(Ncxi'));
   end
       
   GD12=sparse(I,J1,reshape(GD12e,numel(GD12e),1),nn,nn);
   
% The divergence;
ne0=3;
N0c=[1-cx-cy; cx; cy]; e0=ones(ne0,1); nn0=max(max(T(1:ne0,:)));
II01=reshape(T(1:ne0,:),nt*ne0,1);
I01=reshape(kron(II01,e'),nt*ne0,ne,1);
J01=reshape(kron(T(1:ne,:)',e0),ne0*ne*nt,1);
Ke0c=kron((Ic-1)*ne0,e0)+kron(ones(size(Ic)),[1:ne0]');
Dxe=kron(abs(J).*xix,N0c)*(diag(w)*Ncxi')+kron(abs(J).*etx,N0c)*(diag(w)*Ncet');
Dye=kron(abs(J).*xiy,N0c)*(diag(w)*Ncxi')+kron(abs(J).*ety,N0c)*(diag(w)*Ncet');
if nargin>3
    Dxe(Ke0c,:)=(kron(abs(Jc).*xixc,e0).*kron(ones(size(Ic)),N0c))*(diag(w)*Ncxi')...
        + (kron(abs(Jc).*etxc,e0).*kron(ones(size(Ic)),N0c))*(diag(w)*Ncet');
    Dye(Ke0c,:)=(kron(abs(Jc).*xiyc,e0).*kron(ones(size(Ic)),N0c))*(diag(w)*Ncxi')...
        + (kron(abs(Jc).*etyc,e0).*kron(ones(size(Ic)),N0c))*(diag(w)*Ncet');
end

Dx=sparse(I01,J01,reshape(Dxe,numel(Dxe),1),nn0,nn);
Dy=sparse(I01,J01,reshape(Dye,numel(Dye),1),nn0,nn);




% Pressure matrices. pxi pet derivatis wrt xi and eta of linear basis
% functions
% PRESSURE MASS MATRIX
S04=(1/24)*[2 1 1;1 2 1; 1 1 2];
M0e=kron(abs(J),S04);
if nargin>5
    M0e(Ke0c,:)=(kron(ones(size(Ic)),N0c).*kron(abs(Jc),e0))*(diag(w)*N0c');
end
I10=reshape(T(1:ne0,1:nt),nt*ne0,1);
J0=kron(T(1:ne0,1:nt)',e0);
IIq=reshape(I10*e0',numel(I10)*numel(e0),1); 
JJq=reshape(J0,numel(J0),1);
Mp=sparse(IIq,JJq,reshape(M0e,numel(M0e),1),nn0,nn0);

% PRESSURE - LAPLACIAN
ex0=[-1 1 0]'; ey0=[-1 0 1]';
S01=0.5*ex0*ex0'; S03=0.5*ey0*ey0'; S02=0.5*ex0*ey0' + 0.5*ey0*ex0';

   S0e=(kron(abs(J).*(xix.^2+xiy.^2),S01)+...
       kron(abs(J).*(etx.^2+ety.^2),S03)+...
       kron(abs(J).*(xix.*etx+xiy.*ety),S02));
   if nargin>5
       S0e(Ke0c,:)=kron(Jc.*(xixc.^2+xiyc.^2),ex0)*(w(:)*ex0')...
           + kron(Jc.*(etxc.^2+etyc.^2),ey0)*(w(:)*ey0')...
           + kron(Jc.*(xixc.*etxc+xiyc.*etyc),ex0)*(w(:)*ey0')...
            + kron(Jc.*(xixc.*etxc+xiyc.*etyc),ey0)*(w(:)*ex0');
   end
Sp=sparse(IIq,JJq,reshape(S0e,numel(S0e),1),nn0,nn0);

