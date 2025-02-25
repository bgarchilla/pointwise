function [pn,pun,pvn]=uwv(u,v,T,z,J,r)
% J determinants
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
nn=length(z);
if r<=2;
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
else
    %Gauss de 6 nodos.

    %nodos de cuadratura
    c1=0.238619186083197;
    c2=0.661209386466265;
    c3=0.932469514203152;
    c=[-c3 -c2 -c1 c1 c2 c3]';
    c=0.5+c/2;
    % Pesos en los nodos de cuadratura
    w1=0.467913934572691;
    w2=0.360761573048139;
    w3=0.171324492379170;
    w=[w3 w2 w1 w1 w2 w3]/2;w=w';
    w1=w;
    %formula exacta hasta el grado 11
    ry=(1-c)*c';
    rx=c*ones(size(c'));
    rx=reshape(rx,36,1);
    ry=reshape(ry,36,1);

    wx=w*ones(size(w'));
    wy=(1-c)*w';
    w=wx.*wy;
    w=reshape(w,36,1);
    cx=rx';cy=ry';
end
xi=cx; eta=cy;
if r==1
    % Matrices en tri\'angulo de referencia
    % r=1;
    ne=3; % ne = nodos por elemento
    S4=(1/24)*[2 1 1; 1 2 1; 1 1 2];
    Nc=[1-xi-eta; xi; eta];
%     Ncxi=[-ones(size(xi));ones(size(xi)); zeros(size(eta))];
%     Ncet=[-ones(size(eta));zeros(size(xi)); ones(size(eta))];
elseif r==2
    ne=6;  % ne = nodos por elemento
    S4=[ 6 -1 -1  0 -4  0; -1  6 -1  0  0 -4; -1 -1  6 -4  0  0; ...
        0  0 -4 32 16 16; -4  0  0 16 32 16; 0 -4  0 16 16 32]/360;
    Nc=[(1-xi-eta).*(1-2*xi-2*eta); xi.*(2*xi-1); eta.*(2*eta-1); ...
        4*xi.*(1-xi-eta); 4*xi.*eta; 4*eta.*(1-xi-eta)];
%     Ncxi=[-(1-2*xi-2*eta)-2*(1-xi-eta);(2*xi-1)+2*xi;zeros(size(eta));...
%         -4*xi+4*(1-xi-eta); 4*eta; -4*eta];
%     Ncet=[-(1-2*xi-2*eta)-2*(1-xi-eta);zeros(size(xi)); (2*eta-1)+2*eta; ...
%         -4*xi; 4*xi; -4*eta+4*(1-xi-eta)];
elseif r==3
    ne=10;  % ne = nodos por elemento

    S4=[76    11    11    18     0    27    27     0    18    36
        11    76    11     0    18    18     0    27    27    36
        11    11    76    27    27     0    18    18     0    36
        18    0     27   540  -189  -135   -54  -135   270   162
        0   18     27  -189   540   270  -135   -54  -135   162
        27   18      0  -135   270   540  -189  -135   -54   162
        27    0     18   -54  -135  -189   540   270  -135   162
        0   27     18  -135   -54  -135   270   540  -189   162
        18   27      0   270  -135   -54  -135  -189   540   162
        36   36     36   162   162   162   162   162   162  1944]/13440;
    Nc=[nc1(cx,cy); nc2(cx,cy);nc3(cx,cy);nc4(cx,cy);nc5(cx,cy);...
        nc6(cx,cy);nc7(cx,cy);nc8(cx,cy);nc9(cx,cy);nc10(cx,cy)];
%     Ncxi=[nc1x(cx,cy); nc2x(cx,cy);nc3x(cx,cy);nc4x(cx,cy);nc5x(cx,cy);...
%         nc6x(cx,cy);nc7x(cx,cy);nc8x(cx,cy);nc9x(cx,cy);nc10x(cx,cy)];
%     Ncet=[nc1y(cx,cy); nc2y(cx,cy);nc3y(cx,cy);nc4y(cx,cy);nc5y(cx,cy);...
%         nc6y(cx,cy);nc7y(cx,cy);nc8y(cx,cy);nc9y(cx,cy);nc10y(cx,cy)];;
else
    error('r debe ser 1, 2, o 3')
end
e=ones(ne,1);
U=u(T(1:ne,:)'); V=v(T(1:ne,:)');
Uc=U*Nc; Vc=V*Nc;
Pc= (kron(abs(J),Nc).*kron(Vc.*(Uc.^2),e))*w;
II=T(1:ne,:); II=II(:);
pn=full(sparse(II,ones(size(II)),Pc,nn,1));

if nargout>1
    Pcu=kron(abs(J),Nc).*kron(Vc.*(2*Uc),e)*(diag(w)*Nc');
    Pcv=kron(abs(J),Nc).*kron(Uc.^2,e)*(diag(w)*Nc');
    JJ=kron(T(1:ne,:)',e); I2=kron(II,e');
    pun=sparse(I2,JJ,Pcu,nn,nn); pvn=sparse(I2,JJ,Pcv,nn,nn);
end
end

function f=nc1(x,y)
L=(1-x-y);
f=0.5*(3*L-1).*(3*L-2).*L;
end

function f=nc10(x,y)
L1=(1-x-y);L2=x;L3=y;
f=27*L1.*L2.*L3;
end

function f=nc10x(x,y)
L1=(1-x-y);L2=x;L3=y;
%f=27*L1.*L2.*L3;
f=27*(-L2.*L3+L1.*L3);
end

function f=nc10y(x,y)
L1=(1-x-y);L2=x;L3=y;
%f=27*L1.*L2.*L3;
f=27*(-L2.*L3+L1.*L2);
end

function f=nc1x(x,y)
L=(1-x-y);
Lx=-1;
%f=0.5*(3*L-1).*(3*L-2).*L;

f=(0.5*Lx)*(3*(3*L-2).*L+3*(3*L-1).*L+(3*L-1).*(3*L-2));
end

function f=nc1y(x,y)
L=(1-x-y);
Lx=-1;
%f=0.5*(3*L-1).*(3*L-2).*L;

f=(0.5*Lx)*(3*(3*L-2).*L+3*(3*L-1).*L+(3*L-1).*(3*L-2));
end

function f=nc2(x,y)
L=x;
f=0.5*(3*L-1).*(3*L-2).*L;
end

function f=nc2x(x,y)
L=x;Ld=1;
%f=0.5*(3*L-1).*(3*L-2).*L;
f=(0.5*Ld)*(3*(3*L-2).*L+3*(3*L-1).*L+(3*L-1).*(3*L-2));
end

function f=nc2y(x,y)
L=x;Ld=0;
%f=0.5*(3*L-1).*(3*L-2).*L;
f=zeros(size(x));
end

function f=nc3(x,y)
L=y;
f=0.5*(3*L-1).*(3*L-2).*L;
end

function f=nc3x(x,y)
L=y;Lx=0;
%f=0.5*(3*L-1).*(3*L-2).*L;
f=zeros(size(x));
end

function f=nc3y(x,y)
L=y;Ly=1;
%f=0.5*(3*L-1).*(3*L-2).*L;
f=(0.5*Ly)*(3*(3*L-2).*L+3*(3*L-1).*L+(3*L-1).*(3*L-2));
end

function f=nc4(x,y)
L1=(1-x-y);L2=x;
f=4.5*L1.*L2.*(3*L1-1);
end

function f=nc4x(x,y)
L1=(1-x-y);L2=x;
L1d=-ones(size(x));L2d=ones(size(x));
%f=4.5*L1.*L2.*(3*L1-1);
f=4.5*(L1d.*L2.*(3*L1-1)+L1.*L2d.*(3*L1-1)+L1.*L2.*(3*L1d));
end

function f=nc4y(x,y)
L1=(1-x-y);L2=x;
L1d=-ones(size(x));L2d=zeros(size(x));
%f=4.5*L1.*L2.*(3*L1-1);
f=4.5*(L1d.*L2.*(3*L1-1)+L1.*L2.*(3*L1d));
end

function f=nc5(x,y)
L1=(1-x-y);L2=x;
f=4.5*L1.*L2.*(3*L2-1);
end

function f=nc5x(x,y)
L1=(1-x-y);L2=x;
L1d=-ones(size(x));L2d=ones(size(x));
%f=4.5*L1.*L2.*(3*L2-1);
f=4.5*(L1d.*L2.*(3*L2-1)+L1.*L2d.*(3*L2-1)+L1.*L2.*(3*L2d));
end

function f=nc5y(x,y)
L1=(1-x-y);L2=x;
L1d=-ones(size(x));L2d=zeros(size(x));
%f=4.5*L1.*L2.*(3*L2-1);
f=4.5*(L1d.*L2.*(3*L2-1)+L1.*L2.*(3*L2d));
end

function f=nc6(x,y)
L1=x;L2=y;
f=4.5*L1.*L2.*(3*L1-1);
end

function f=nc6x(x,y)
L1=x;L2=y;
L1d=ones(size(x));L2d=zeros(size(x));
%f=4.5*L1.*L2.*(3*L1-1);
f=4.5*(L1d.*L2.*(3*L1-1)+L1.*L2.*(3*L1d));
end

function f=nc6y(x,y)
L1=x;L2=y;
L1d=zeros(size(x));L2d=ones(size(x));
%f=4.5*L1.*L2.*(3*L1-1);
f=4.5*(L1.*L2d.*(3*L1-1));
end

function f=nc7(x,y)
L1=x;L2=y;
f=4.5*L1.*L2.*(3*L2-1);
end

function f=nc7x(x,y)
L1=x;L2=y;
L1d=ones(size(x));
%f=4.5*L1.*L2.*(3*L2-1);
f=4.5*L1d.*L2.*(3*L2-1);
end

function f=nc7y(x,y)
L1=x;L2=y;
L1d=ones(size(x));
%f=4.5*L1.*L2.*(3*L2-1);
f=4.5*(L1.*(3*L2-1)+3*L1.*L2);
end

function f=nc8(x,y)
L1=y;L2=1-x-y;
f=4.5*L1.*L2.*(3*L1-1);
end

function f=nc8x(x,y)
L1=y;L2=1-x-y;
%f=4.5*L1.*L2.*(3*L1-1);
f=-4.5*(L1.*(3*L1-1));
end

function f=nc8y(x,y)
L1=y;L2=1-x-y;
%f=4.5*L1.*L2.*(3*L1-1);
f=4.5*(L2.*(3*L1-1)-L1.*(3*L1-1)+3*L1.*L2);
end

function f=nc9(x,y)
L1=y;L2=1-x-y;
f=4.5*L1.*L2.*(3*L2-1);
end

function f=nc9x(x,y)
L1=y;L2=1-x-y;
%f=4.5*L1.*L2.*(3*L2-1);
f=-4.5*(L1.*(3*L2-1)+3*L1.*L2);
end

function f=nc9y(x,y)
L1=y;L2=1-x-y;
%f=4.5*L1.*L2.*(3*L2-1);
f=4.5*(L2.*(3*L2-1)-L1.*(3*L2-1)-3*L1.*L2);
end


