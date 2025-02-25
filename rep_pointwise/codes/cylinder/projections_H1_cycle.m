clear all
disp('Loading data. Please wait ...')
load ../../data/input_data/cylinder/the_snapshots.mat
load ../../data/output_data/cylinder/pod_basis_H1
load ../../data/output_data/cylinder/deriv_data_L2
disp(' ... done')
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

[~,Mn]=makronvq_skw(T,z,1,nTc,phi,ccphi);
incre=increPOD;

U=U(:,1:incre:end); V=V(:,1:incre:end);
u0=mean(U(:,2:end),2); v0=mean(V(:,2:end),2); 
% U=U(:,2:end); V=V(:,2:end); N=size(U,2); Nkeep=N
% UV=[U;V]; DUV=[UV(:,2:end)-UV(:,1:end-1),UV(:,1)-UV(:,end)]/(tiempos(end)/N);
N=size(U,2)-1; Nkeep=N;
UV=[U;V]; DUV=[UV(:,2:end)-UV(:,1:end-1)]/(tiempos(end)/N);
M=kron(eye(2),Mn);
S=kron(eye(2),Sn);

% r=19;

nPhi=sqrt(abs(sum(Phi.*(M*Phi))));
figure(21); clf; subplot(2,1,1);
semilogy(nPhi,'b-','LineWidth',1.5); 
set(gca,'FontSize',16);
title('$\left\| \varphi^k\right\|$ for $X=H^1_0$','FontSize',20,'Interpreter','LaTex')
axis([0 246 1e-3 1e-1])


rhos=zeros(70,4); err_PD1=zeros(70,1); rats=rhos; mus=rhos;
erres=[11:80]';

disp('Computing L2(L2) projection errors and mu_m for r=11 to r=80 (this may take up to a minute) ...') 

for j=1:70
r=erres(j); gamma_r=norm(ss(r+1:end));
gamma_rr=norm(ss(r+1:end).*nPhi(r+1:end))*sqrt(tiempos(end));
C=Phi(:,1:r)'*(S*DUV);
E=DUV - (Phi(:,1:r)*C);
e=sqrt(abs(sum(E.*(M*E)))); nUV=sqrt(abs(sum(UV.*(M*UV))));
format short e, err_PD1(j)=norm(e)*sqrt(tiempos(end)/N); format short

% computing the mu_m

the_errs=zeros(1,5);

m=2;
dembk=tiempos2; Z=[Utt; Vtt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));

m=3;
dembk=tiempos3; Z=[Uttt; Vttt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));

m=4;
dembk=tiempos4; Z=[Utttt; Vtttt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));

m=5;
dembk=tiempos5; Z=[Uttttt; Vttttt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));

cA=sqrt(2+1/sqrt(2)); 
the_errs=2*the_errs;

format short e, mus(j,:)=the_errs(2:end); format short
end

rats=mus./err_PD1;

[Mrat,iMrat]=max(rats)
disp(' .... done.')
disp('Them maximum values of the overestimation ratios for m=2,3,4,5 are ...')
disp(Mrat)
disp('... and they are reached at the following values of r:')
disp(10+iMrat)
[mrat,imrat]=min(rats);
disp('Them minimum values of the overestimation ratios for m=2,3,4,5 are ...')
disp(mrat)
disp('... and they are reached at the following values of r:')
disp(10+imrat)
mrats=mean(rats(:));
disp(strcat(['The mean overestimation ratio is = ',blanks(1),num2str(mrats)]))


figure(31); clf
semilogy(erres,err_PD1,'b-','LineWidth',1.5); hold on;
semilogy(erres,mus(:,1),'c-','LineWidth',1.5);
semilogy(erres,mus(:,2),'m-','LineWidth',1.5);
semilogy(erres,mus(:,3),'r-','LineWidth',1.5);
semilogy(erres,mus(:,4),'k-','LineWidth',1.5);
set(gca,'FontSize',18);
xlabel('$r$','FontSize',24,'Interpreter','LaTeX')


% r=30;

% r=30; gamma_r=norm(ss(r+1:end));
% gamma_rr=norm(ss(r+1:end).*nPhi(r+1:end))*sqrt(tiempos(end)/N)
% C=Phi(:,1:r)'*(S*(DUV));
% E=DUV - (Phi(:,1:r)*C);
% e=sqrt(abs(sum(E.*(M*E)))); nUV=sqrt(abs(sum(UV.*(M*UV))));
% err_PD1=norm(e)*sqrt(tiempos(end)/N)
% 
% % computing the rhs of 7
% 
% the_errs=zeros(1,5);
% 
% m=2;
% dembk=tiempos2; Z=[Utt; Vtt];
% C=Phi(:,1:r)'*(S*Z); 
% E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
% dt=mean(diff(dembk));
% the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
% 
% m=3;
% dembk=tiempos3; Z=[Uttt; Vttt];
% C=Phi(:,1:r)'*(S*Z); 
% E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
% dt=mean(diff(dembk));
% the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
% 
% m=4;
% dembk=tiempos4; Z=[Utttt; Vtttt];
% C=Phi(:,1:r)'*(S*Z); 
% E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
% dt=mean(diff(dembk));
% the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
% 
% m=5;
% dembk=tiempos5; Z=[Uttttt; Vttttt];
% C=Phi(:,1:r)'*(S*Z); 
% E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
% dt=mean(diff(dembk));
% the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
% 
% cA=sqrt(2+1/sqrt(2)); 
% the_errs=2*the_errs;
% 
% mus=the_errs(2:end)
% rats=mus/err_PD1
% 
