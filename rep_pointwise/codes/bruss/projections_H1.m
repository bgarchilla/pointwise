clear all
disp('Loading data. Please wait ...')
load ../../data/output_data/bruss/the_snapshots
load ../../data/output_data/bruss/pod_basis_H1
load ../../data/output_data/bruss/deriv_data
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

Mn=Mhn;
incre=increPOD;

U=U(:,1:incre:end); V=V(:,1:incre:end);
u0=mean(U(:,2:end),2); v0=mean(V(:,2:end),2); 
U=U(:,2:end); V=V(:,2:end); N=size(U,2); Nkeep=N
UV=[U;V]; DUV=[UV(:,2:end)-UV(:,1:end-1),UV(:,1)-UV(:,end)]/(tiempos(end)/N);
M=kron(eye(2),Mn);
S=kron(eye(2),Sn);

% r=24;

nPhi=sqrt(abs(sum(Phi.*(M*Phi))));
figure(21); clf; subplot(2,1,1);
semilogy(nPhi,'b-','LineWidth',1.5); hold on
plot([0 256],sqrt(2)*[1 1]/pi,'r:','LineWidth',1.5)
set(gca,'FontSize',20,'LineWidth',1.2);
yticks([0.01 0.1])
title('$\left\| \varphi^k\right\|$ for $X=H^1_0$ and constant $C_P$','FontSize',25,'Interpreter','LaTex')
text(124,2e-1,'$C_P$','FontSize',26,'Interpreter','LaTex')
axis([0 256 0.9*min(nPhi) 0.6])

print -depsc ../../figures/bruss/normasphi.eps
disp('Generated figure normasphi.eps. Check in figures/bruss')

disp(' ')
disp(' ')
disp(' ')

r=24; 
disp(strcat(['For r =',blanks(1),num2str(r),blanks(1),'the values are ... ']))

gamma_r=norm(ss(r+1:end));
gamma_rr=norm(ss(r+1:end).*nPhi(r+1:end))*sqrt(tiempos(end)/N);
C=Phi(:,1:r)'*(S*DUV);
E=DUV - (Phi(:,1:r)*C);
e=sqrt(abs(sum(E.*(M*E)))); nUV=sqrt(abs(sum(UV.*(M*UV))));
disp('L2(L2) norm of projection errors of first-order differences:')
format short e, err_PD1=norm(e)*sqrt(tiempos(end)/N), format short

% computing the mu_m

the_errs=zeros(1,5); the_orrs=the_errs;

m=2;
dembk=tiempos; Z=[Utt; Vtt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
e=sqrt(abs(sum(E.*(S*E))));
the_orrs(m)=2*((sqrt(dt)*norm(e))^(1/(m)))*(gamma_r^(1-1/(m)));

m=3;
dembk=tiempos; Z=[Uttt; Vttt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
e=sqrt(abs(sum(E.*(S*E))));
the_orrs(m)=2*((sqrt(dt)*norm(e))^(1/(m)))*(gamma_r^(1-1/(m)));

m=4;
dembk=tiempos; Z=[Utttt; Vtttt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
e=sqrt(abs(sum(E.*(S*E))));
the_orrs(m)=2*((sqrt(dt)*norm(e))^(1/(m)))*(gamma_r^(1-1/(m)));

m=5;
dembk=tiempos; Z=[Uttttt; Vttttt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
e=sqrt(abs(sum(E.*(S*E))));
the_orrs(m)=2*((sqrt(dt)*norm(e))^(1/(m)))*(gamma_r^(1-1/(m)));

cA=sqrt(2+1/sqrt(2)); 
the_errs=2*the_errs;
the_orrs=2*the_orrs;

disp('the mu_m and rho_m for m=2,3,4,5:')
format short e, mus=the_errs(2:end), format short
format short e, rhos=the_orrs(2:end), format short

disp('and their respecive overstimation ratios:')
rats=mus/err_PD1
rots=rhos/err_PD1

disp(' ')
disp(' ')
disp(' ')

% r=30;

r=31;
disp(strcat(['For r =',blanks(1),num2str(r),blanks(1),'the values are ... ']))

gamma_r=norm(ss(r+1:end));
gamma_rr=norm(ss(r+1:end).*nPhi(r+1:end))*sqrt(tiempos(end)/N);
C=Phi(:,1:r)'*(S*(DUV));
E=DUV - (Phi(:,1:r)*C);
e=sqrt(abs(sum(E.*(M*E)))); nUV=sqrt(abs(sum(UV.*(M*UV))));
% err_PD1=norm(e)*sqrt(tiempos(end)/N)
disp('L2(L2) norm of projection errors of first-order differences:')
format short e, err_PD1=norm(e)*sqrt(tiempos(end)/N), format short

% computing the mu_m
the_errs=zeros(1,5); the_orrs=the_errs;

m=2;
dembk=tiempos; Z=[Utt; Vtt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
e=sqrt(abs(sum(E.*(S*E))));
the_orrs(m)=2*((sqrt(dt)*norm(e))^(1/(m)))*(gamma_r^(1-1/(m)));

m=3;
dembk=tiempos; Z=[Uttt; Vttt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
e=sqrt(abs(sum(E.*(S*E))));
the_orrs(m)=2*((sqrt(dt)*norm(e))^(1/(m)))*(gamma_r^(1-1/(m)));

m=4;
dembk=tiempos; Z=[Utttt; Vtttt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
e=sqrt(abs(sum(E.*(S*E))));
the_orrs(m)=2*((sqrt(dt)*norm(e))^(1/(m)))*(gamma_r^(1-1/(m)));

m=5;
dembk=tiempos; Z=[Uttttt; Vttttt];
C=Phi(:,1:r)'*(S*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(m)))*(gamma_rr^(1-1/(m)));
e=sqrt(abs(sum(E.*(S*E))));
the_orrs(m)=2*((sqrt(dt)*norm(e))^(1/(m)))*(gamma_r^(1-1/(m)));

cA=sqrt(2+1/sqrt(2)); 
the_errs=2*the_errs;
the_orrs=2*the_orrs;

disp('the mu_m and rho_m for m=2,3,4,5:')
format short e, mus=the_errs(2:end), format short
format short e, rhos=the_orrs(2:end), format short

disp('and their respecive overstimation ratios:')
rats=mus/err_PD1
rots=rhos/err_PD1

