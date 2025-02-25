clear all
% load the_snapshots
load ../../data/output_data/bruss/the_snapshots.mat

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

incre=1;
U=U(:,1:incre:end); V=V(:,1:incre:end);
u0=mean(U(:,2:end),2); v0=mean(V(:,2:end),2); 
U=U(:,2:end); V=V(:,2:end); N=size(U,2);
w0=[u0;v0];

UV=[U;V]-w0;
I=Iu; nI=length(I); % note that Iu=Iv;
Sn=Shn;

disp('computing chol of S ...')
[R,iflag,P]=chol(Sn(I,I));
disp('... done')
AA=[R*(P'*UV(I,:));R*(P'*UV(length(z)+I,:))];
sizeAA=size(AA)
disp('computing SVD ...')
[Wus,SS,Wvs]=svd(full(AA),0);
disp('... done')
ss=diag(SS)/sqrt(N);
figure(3);clf; semilogy(ss,'b-')
k10=find(ss>1e-10);
Wu=zeros(2*length(z),size(Wus,2));
Wu(I,:)=P*(R\Wus(1:nI,:));Wu(length(z)+I,:)=P*(R\Wus(1+nI:2*nI,:)); clear Wus Wvs R
Phi=Wu;

increPOD=incre;
save ../../data/output_data/bruss/pod_basis_H1 I Sn w0 Phi ss increPOD
