clear all
disp('Loading data. Please wait ...')
load ../../data/input_data/cylinder/the_snapshots
load ../../data/output_data/cylinder/pod_basis_L2
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

incre=increPOD;

U=U(:,1:incre:end); V=V(:,1:incre:end);
u0=mean(U(:,2:end),2); v0=mean(V(:,2:end),2); 
U=U(:,2:end); V=V(:,2:end); N=size(U,2);
UV=[U;V];
M=kron(eye(2),Mn);

% r=12;
rhs_of_8=zeros(70,4); max_err_L2=zeros(70,1); rats=rhs_of_8;
erres=[11:80]';

%r=24;
disp('Computing maximum projection errors and rhs of (8) for r=11 to r=80 (this may take up to a minute) ...') 
for j=1:70
    r=erres(j);
    format short e, gamma_r=norm(ss(r+1:end)); format short, gamma_rr=gamma_r*sqrt(tiempos(end));
C=Phi(:,1:r)'*(M*(UV-w0));
E=UV - (w0 + Phi(:,1:r)*C);
e=sqrt(abs(sum(E.*(M*E)))); nUV=sqrt(abs(sum(UV.*(M*UV))));
format short e, max_err_L2(j)=max(e); format short

% computing the rhs of 7

the_errs=zeros(1,5);

m=2;
dembk=tiempos2; Z=[Utt; Vtt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

m=3;
dembk=tiempos3; Z=[Uttt; Vttt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

m=4;
dembk=tiempos4; Z=[Utttt; Vtttt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

m=5;
dembk=tiempos5; Z=[Uttttt; Vttttt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

cA=sqrt(2+1/sqrt(2)); 
the_errs=sqrt(2)*cA*the_errs + sqrt(2)*gamma_r;
rhs_of_8(j,:)=the_errs(2:end);
rats(j,:)=rhs_of_8(j,:)/max_err_L2(j);

% format short e, rhs_of_7=the_errs(2:end), format short
end

[Mrat,iMrat]=max(rats);
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
% rats=rhs_of_8/max_err_L2
mrats=mean(rats(:));
disp(strcat(['The mean overestimation ratio is = ',blanks(1),num2str(mrats)]))


figure(51); clf
semilogy(erres,max_err_L2,'b-','LineWidth',1.5); hold on;
semilogy(erres,rhs_of_8(:,1),'c-','LineWidth',1.5);
semilogy(erres,rhs_of_8(:,2),'m-','LineWidth',1.5);
semilogy(erres,rhs_of_8(:,3),'r-','LineWidth',1.5);
semilogy(erres,rhs_of_8(:,4),'k-','LineWidth',1.5);
set(gca,'FontSize',18)



