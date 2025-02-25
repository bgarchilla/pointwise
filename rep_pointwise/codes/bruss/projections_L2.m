clear all
disp('Loading data. Please wait ...')
load ../../data/output_data/bruss/the_snapshots
load ../../data/output_data/bruss/pod_basis_L2
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

incre=increPOD;

U=U(:,1:incre:end); V=V(:,1:incre:end);
u0=mean(U(:,2:end),2); v0=mean(V(:,2:end),2); 
U=U(:,2:end); V=V(:,2:end); N=size(U,2);
UV=[U;V];
M=kron(eye(2),Mn);

% r=24;

r=24;
disp(strcat(['For r =',blanks(1),num2str(r),blanks(1),'the values are ... ']))

gamma_r=norm(ss(r+1:end)), gamma_rr=gamma_r*sqrt(tiempos(end));
C=Phi(:,1:r)'*(M*(UV-w0));
E=UV - (w0 + Phi(:,1:r)*C);
e=sqrt(abs(sum(E.*(M*E)))); nUV=sqrt(abs(sum(UV.*(M*UV))));
max_err_L2=max(e)

% computing the rhs of 7

the_errs=zeros(1,5);

m=2;
dembk=tiempos; Z=[Utt; Vtt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

m=3;
dembk=tiempos; Z=[Uttt; Vttt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

m=4;
dembk=tiempos; Z=[Utttt; Vtttt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

m=5;
dembk=tiempos; Z=[Uttttt; Vttttt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

cA=sqrt(2+1/sqrt(2)); 
the_errs=sqrt(2)*cA*the_errs + sqrt(2)*gamma_r;

format short e, rhs_of_8=the_errs(2:end), format short
disp(' ... and the overestimatin ratios are ... ')
rats=rhs_of_8/max_err_L2

disp(' ')
disp(' ')
disp(' ')

% r=31;

r=31;
disp(strcat(['For r =',blanks(1),num2str(r),blanks(1),'the values are ... ']))

gamma_r=norm(ss(r+1:end)), gamma_rr=gamma_r*sqrt(tiempos(end));
C=Phi(:,1:r)'*(M*(UV-w0));
E=UV - (w0 + Phi(:,1:r)*C);
e=sqrt(abs(sum(E.*(M*E)))); nUV=sqrt(abs(sum(UV.*(M*UV))));
max_err_L2=max(e)

% computing the rhs of 7

the_errs=zeros(1,5);

m=2;
dembk=tiempos; Z=[Utt; Vtt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

m=3;
dembk=tiempos; Z=[Uttt; Vttt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

m=4;
dembk=tiempos; Z=[Utttt; Vtttt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

m=5;
dembk=tiempos; Z=[Uttttt; Vttttt];
C=Phi(:,1:r)'*(M*Z); 
E=Z-Phi(:,1:r)*C; e=sqrt(abs(sum(E.*(M*E))));
dt=mean(diff(dembk));
the_errs(m)=((sqrt(dt)*norm(e))^(1/(2*m)))*(gamma_rr^(1-1/(2*m)));

cA=sqrt(2+1/sqrt(2)); 
the_errs=sqrt(2)*cA*the_errs + sqrt(2)*gamma_r;

rhs_of_8=the_errs(2:end);
format short e, rhs_of_8=the_errs(2:end), format short

disp(' ... and the overestimatin ratios are ... ')
rats=rhs_of_8/max_err_L2

