% to compute the derivatives wrt t of the snapshots
clear all 
% load the_snapshots
load ../../data/output_data/bruss/the_snapshots.mat

% loaded variables are
% Tri 6 x nt with the triangulation for quadratic elements (nt=12800)
%            Column j contains the indexes of the nodes of triangle j
% z   nn x 2 xy coordinates of nodes in the triangulation
%            Row i contais the xy coordinates of node i
% r       degree of piecewise polynomials of the FEM approximation (r=2)
% U,V     nn x 257 matrices with the u and v components of the snapshots
% tp      period of the limit cycle.
% tiempos 1 x 257 with times of the snapshots tiempos=[0:tp/1024:tp];
% nuu, nuv diffusion for components u and v of the brusselator
%         (nuu=nuv=0.002)
% pA, pB  parameters for the Brusselator pA=1; pB=3
% Iu,Iv   indexes of nodes NOT on the Dirichlet boundary for u and v.
% Mhn,Shn mass and stiffness matrix
% J   nt x 1 matrix with the jacobians of the mapping from the reference
%         triangle to each of the triangles of the triangulation.
% IS IE IN IW idexes of nodes on sides south, east,  nort and west of the
%         of the triangulation


% values and matrices required for the computation of derivatives

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

nIu=length(Iu); nIv=length(Iv);
Ahn=[-nuu*Shn-(pB+1)*Mhn,sparse(1,1,0,nn,nn);...
    pB*Mhn, -nuv*Shn];
Ah=[-nuu*Shn(Iu,Iu)-(pB+1)*Mhn(Iu,Iu),sparse(1,1,0,nIu,nIv);...
    pB*Mhn(Iv,Iu), -nuv*Shn(Iv,Iv)];
Mh=[Mhn(Iu,Iu),sparse(1,1,0,nIu,nIv); sparse(1,1,0,nIv,nIu), Mhn(Iv,Iv)];

pAh=pA*(Mhn(Iu,:)*ones(nn,1)); % source of pA in all domain
% pAh=pA*[Ahn(Iu,IW);Ahn(nn+Iv,IW)]*ones(length(IW),1); % u = pA on West boundary
pAh=[pAh;zeros(length(Iv),1)];
IWS=unique([IN;IE]);
pAh=pAh+[Ahn([Iu;nn+Iv],IWS)*(pA*ones(size(IWS)))] +...
    [Ahn([Iu;nn+Iv],nn+IWS)*((pB/pA)*ones(size(IWS)))]; 

% computing the derivatives
disp('Computing time derivatives (this may take several minutes) ...')
tic

Ut=zeros(size(U)); Vt=zeros(size(V)); 
Utt=zeros(size(U)); Vtt=zeros(size(V)); 
Uttt=zeros(size(U)); Vttt=zeros(size(V)); 
Utttt=zeros(size(U)); Vtttt=zeros(size(V)); 
Uttttt=zeros(size(U)); Vttttt=zeros(size(V)); 
F=zeros(length(Iu)+length(Iv),size(U,2));
nU2=size(U,2);
Ft=F; Ftt=Ft; Fttt=Ft; Ftttt=Ft;
for j=1:size(U,2)
    y=[U(Iu,j);V(Iv,j)];
    [ftttt,fttt,ftt,ft,f]=derivs(0,y,Ah,pAh,Iu,Iv,Tri,z,J,r,pA,pB,Mh);
    F(:,j)=f; Ft(:,j)=ft; Ftt(:,j)=ftt; Fttt(:,j)=fttt;  Ftttt(:,j)=ftttt;
end
S=F;
Ut(Iu,:)=S(1:nIu,1:nU2); Vt(Iv,:)=S(nIu+1:nIu+nIv,1:nU2); 
S=Ft;
Utt(Iu,:)=S(1:nIu,1:nU2); Vtt(Iv,:)=S(nIu+1:nIu+nIv,1:nU2); 
S=Ftt;
Uttt(Iu,:)=S(1:nIu,1:nU2); Vttt(Iv,:)=S(nIu+1:nIu+nIv,1:nU2); 
S=Fttt;
Utttt(Iu,:)=S(1:nIu,1:nU2); Vtttt(Iv,:)=S(nIu+1:nIu+nIv,1:nU2); 
S=Ftttt;
Uttttt(Iu,:)=S(1:nIu,1:nU2); Vttttt(Iv,:)=S(nIu+1:nIu+nIv,1:nU2); 

disp('... done')
elapsed_time_computing_derivatives=toc

% save deriv_data tiempos Ut Vt Utt Vtt Uttt Vttt Utttt Vtttt Uttttt Vttttt
save ../../data/output_data/bruss/deriv_data tiempos Ut Vt Utt Vtt Uttt Vttt Utttt Vtttt Uttttt Vttttt