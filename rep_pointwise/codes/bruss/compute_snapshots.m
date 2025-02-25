clear all 
% load the daa
load ../../data/input_data/bruss/data_for_snapshots.mat
% loaded variables are
% Tri 6 x nt with the triangulation for quadratic elements (nt=12800)
% TT  3 x 4*nt with triangulation for ploting the snapshots
%            Column j contains the indexes of the nodes of triangle j
% z   nn x 2 xy coordinates of nodes in the triangulation
%            Row i contais the xy coordinates of node i
% u,v     nn x 1 vectors with the u and v components of the initial values
%         for the snapshots snapshots
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
% TOL   Tolerance for the relative local erros for time integration
%        (TOL=1e-9). Do not change if you intend to replicate results in
%              B. Garcia-Archilla and J. Novo,
%              Pointwise error bounds in POD methods without difference quotients,
%              Journal of Scientific Computing (to appear)
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


r=2; % degree of pw polynomials used
nIu=length(Iu); nIv=length(Iv);
[Shn,Mhn,J]=matrices(Tri,z,r); % stiffnes and mass matrices and Jacobians of triangulation
nn=length(z);

% Ah matrix of the linear part of of the Brusselator (
Ahn=[-nuu*Shn-(pB+1)*Mhn,sparse(1,1,0,nn,nn);...
    pB*Mhn, -nuv*Shn]; % for all nodes
Ah=[-nuu*Shn(Iu,Iu)-(pB+1)*Mhn(Iu,Iu),sparse(1,1,0,nIu,nIv);...
    pB*Mhn(Iv,Iu), -nuv*Shn(Iv,Iv)]; % for active nodes

% Mh mass matrix of the brusselator (only active nodes)
Mh=[Mhn(Iu,Iu),sparse(1,1,0,nIu,nIv); sparse(1,1,0,nIv,nIu), Mhn(Iv,Iv)];
pAh=pA*(Mhn(Iu,:)*ones(nn,1)); % source of pA in all domain (only active nodes)
pAh=[pAh;zeros(length(Iv),1)];
IWS=unique([IN;IE]);
pAh=pAh+[Ahn([Iu;nn+Iv],IWS)*(pA*ones(size(IWS)))] +...
    [Ahn([Iu;nn+Iv],nn+IWS)*((pB/pA)*ones(size(IWS)))]; 


nn=length(z);
y0=[u(Iu);v(Iv)];
Nsn=1024; dt=tp/1024;
tiempos=[0:dt:tp];
fun=@(t,y)bruss(t,y,Ah,pAh,Iu,Iv,Tri,z,J,r,pA,pB);
jfun=@(t,y)jbruss(t,y,Ah,pAh,Iu,Iv,Tri,z,J,r,pA,pB);
misops=odeset('RelTol',TOL,'AbsToL',TOL/1000,'Mass',Mh,'Jacobian',jfun,'Stats','on');
disp('Computing the snapshots (this may take up to one minute) ...')
tic; [T,Y]=ode15s(fun,tiempos,y0,misops); toc;
disp(' ... done')
U=zeros(nn,length(T)); V=U;
U(Iu,:)=Y(:,1:nIu)'; V(Iv,:)=Y(:,nIu+1:nIu+nIv)'; 
U(IWS,:)=pA; V(IWS,:)=(pB/pA);
normaU=sqrt(abs(sum(U.*(Mhn*U)))); normaV=sqrt(abs(sum(V.*(Mhn*V)))); 
norma1U=sqrt(abs(sum(U.*(Shn*U)))); norma1V=sqrt(abs(sum(V.*(Shn*V)))); 

% To view the solution uncomment the following lines;
    % figure(3);
    % for j=1:length(T)
    %     trimesh(TT',z(:,1),z(:,2),U(:,j)); title(num2str(T(j))); axis(aa); drawnow
    % end
    % figure(4); clf;  aa=[0 1 0 1 0 5];
    % for j=1:length(T);
    %     trimesh(TT',z(:,1),z(:,2),V(:,j)); title(num2str(T(j))); axis(aa); drawnow
    % end

 % keeping only 257 snapshots
 U=U(:,1:4:end); V=V(:,1:4:end); tiempos=tiempos(1:4:end);
 save ../../data/output_data/bruss/the_snapshots IS IE IN IW Iu Iv nuu nuv ...
     pA pB Tri z r tp J tiempos U V Shn Mhn
