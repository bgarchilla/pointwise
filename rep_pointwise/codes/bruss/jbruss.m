function D=jbruss(t,y,Ah,pAh,Iu,Iv,T,z,J,r,pA,pB)
% Differential with respect to y of function bruss.m
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
nnu=length(Iu); nnv=length(Iv);
u=y(1:nnu); v=y(nnu+1:nnu+nnv);
vv=zeros(nn,1);
if nargin<11; uu=vv; else; uu=pA*ones(nn,1); vv=(pB/pA)*ones(nn,1); end
uu(Iu)=u; vv(Iv)=v;
[~,pu,pv]=u2v(uu,vv,T,z,J,r);
Dp=[pu(Iu,Iu),pv(Iu,Iv);-pu(Iv,Iu),-pv(Iv,Iv)];
D=Ah+Dp;