function [Ahn,Mhn,J]=matrices(T,z,r)
% Laplacian and mass matrix.
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

   nt=size(T,2); nn=size(z,1);

if r==1
    % Matrices en tri\'angulo de referencia
    % r=1;
    ne=3; % ne = nodos por elemento
    S1=0.5*[1 -1 0; -1 1 0; 0 0 0];
    S12=0.5*[1 0 -1;-1 0 1; 0 0 0];
    S2=S12+S12';
    S3=0.5*[1 0 -1; 0 0 0; -1 0 1];
    S4=(1/24)*[2 1 1; 1 2 1; 1 1 2];
%     Ng=[1-xi-eta; xi; eta];
elseif r==2
        ne=6;  % ne = nodos por elemento
        S1=[ 3 1 0 -4 0 0; 1 3 0 -4 0 0; 0 0 0 0 0 0; -4 -4 0 8 0 0; ...
            0 0 0 0 8 -8; 0 0 0 0 -8 8 ]/6;
        S2=[ 6 1 1 -4 0 -4; 1 0 -1 -4 4 0;  1 -1 0 0 4 -4; -4 -4 0 8 -8 8; ...
            0 4 4 -8 8 -8; -4 0 -4 8 -8 8 ]/6;
        S3=[ 3 0 1 0 0 -4; 0 0 0 0 0 0; 1 0 3 0 0 -4; 0 0 0 8 -8 0; ...
            0 0 0 -8 8 0; -4 0 -4 0 0 8]/6;
        S12=[3 0 1 0 0 -4; 1 0 -1 -4 4 0; 0 0 0 0 0 0; -4 0 0 4 -4 4; ...
            0 0 4 -4 4 -4; 0 0 -4 4 -4 4]/6;
        S4=[ 6 -1 -1  0 -4  0; -1  6 -1  0  0 -4; -1 -1  6 -4  0  0; ...
            0  0 -4 32 16 16; -4  0  0 16 32 16; 0 -4  0 16 16 32]/360;
%         Ng=[(1-xi-eta).*(1-2*xi-2*eta); xi.*(2*xi-1); eta.*(2*eta-1); ...
%             4*xi.*(1-xi-eta); 4*xi.*eta; 4*eta.*(1-xi-eta)];
elseif r==3
        ne=10;  % ne = nodos por elemento
        S1=[34    -7     0   -54    27    -3    -3     3     3     0
            -7    34     0    27   -54     3     3    -3    -3     0
            0     0     0     0     0     0     0     0     0     0
            -54    27     0   135  -108     0     0     0     0     0
            27   -54     0  -108   135     0     0     0     0     0
            -3     3     0     0     0   135   -27    27    27  -162
            -3     3     0     0     0   -27   135  -135    27     0
            3    -3     0     0     0    27  -135   135   -27     0
            3    -3     0     0     0    27    27   -27   135  -162
            0     0     0     0     0  -162     0     0  -162   324]/80;

        S2=[68    -7    -7   -51    30    -6    -6    30   -51     0
            -7     0     7    24   -57    57   -24     0     0     0
            -7     7     0     0     0   -24    57   -57    24     0
            -51    24     0   135  -108    27    27   -27   135  -162
            30   -57     0  -108   135  -135    27   -27   -27   162
            -6    57   -24    27  -135   135    54    27    27  -162
            -6   -24    57    27    27    54   135  -135    27  -162
            30     0   -57   -27   -27    27  -135   135  -108   162
            -51     0    24   135   -27    27    27  -108   135  -162
            0     0     0  -162   162  -162  -162   162  -162   324]/80;

        S3=[34     0    -7     3     3    -3    -3    27   -54     0
            0     0     0     0     0     0     0     0     0     0
            -7     0    34    -3    -3     3     3   -54    27     0
            3     0    -3   135   -27    27    27     0     0  -162
            3     0    -3   -27   135  -135    27     0     0     0
            -3     0     3    27  -135   135   -27     0     0     0
            -3     0     3    27    27   -27   135     0     0  -162
            27     0   -54     0     0     0     0   135  -108     0
            -54     0    27     0     0     0     0  -108   135     0
            0     0     0  -162     0     0  -162     0     0   324]/80;

        S12=[ 68     0   -14     6     6    -6    -6    54  -108     0
            -14     0    14    48  -114   114   -48     0     0     0
            0     0     0     0     0     0     0     0     0     0
            -108     0     0   135   -27    27    27   -27   135  -162
            54     0     0  -189   135  -135    27   -27   -27   162
            -6     0   -48    27  -135   135   135    27    27  -162
            -6     0   114    27    27   -27   135  -135    27  -162
            6     0  -114   -27   -27    27  -135   135   -27   162
            6     0    48   135   -27    27    27  -189   135  -162
            0     0     0  -162   162  -162  -162   162  -162   324]/160;

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
else
    error('r debe ser 1, 2, o 3')
end
   x1=z(:,1);y1=z(:,2);x=x1(T(1:ne,:)');y=y1(T(1:ne,:)');
   clear x1 y1;

% Partial derivatives of x and y wrt xi eta
   xxi=x(:,2)-x(:,1); xet=x(:,3)-x(:,1);
   yxi=y(:,2)-y(:,1); yet=y(:,3)-y(:,1);

% The Jacobians
   J=(x(:,2)-x(:,1)).*(y(:,3)-y(:,1))-(x(:,3)-x(:,1)).*(y(:,2)-y(:,1));
   
   % Parciales de xi y eta con respecto a x e y
   xix=(y(:,3)-y(:,1))./J; xiy=(x(:,1)-x(:,3))./J;
   etx=(y(:,1)-y(:,2))./J; ety=(x(:,2)-x(:,1))./J;
%   J=J(Is); xix=xis(Is); xiy=xiy(Is); etx=etx(Is); ety=ety(Is);

e=ones(ne,1);

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
Ahn=sparse(I,J1,Se(:),nn,nn); Mhn=sparse(I,J1,Me,nn,nn);

