function f=bruss(t,y,Ah,pAh,Iu,Iv,T,z,J,r,pA,pB)
% right-hand side of a quadratic fem discretiztion of the
% Brusselator
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
vv=zeros(length(z),1);
if nargin<11
    uu=vv;
else
    uu=pA*ones(nn,1);
    vv=(pB/pA)*ones(nn,1);
end
uu(Iu)=u; vv(Iv)=v;
pn=u2v(uu,vv,T,z,J,r);
f=Ah*y+[pn(Iu);-pn(Iv)];
if nargin<11
   f(1:length(Iu))=f(1:length(Iu))+pAh;
else
    f=f+pAh;
end