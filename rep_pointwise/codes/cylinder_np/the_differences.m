% load the_snapshots
load ../../data/input_data/cylinder/the_snapshots.mat
% variables are
%  U V P snapshots of the two components of the velocity and the pressure
%  Ut Vt time derivatives of the snapshots
%  T 6 x nt matrix. Column j contains the indexes of the nodes of triangle j
%  z nn x 2 matrix. Row i contains the xy coordinates of node i
%  tp period of the orbit
%  tiempos 1 x 1025 matrix with the times of the snapshots
%  epsilon the viscosity
%  gamma the grad-div parameter
%  nTc idexes of the triangles whos side 2-3 is curved
%  phi @phiiso map fot the triangles with curved side
%  ccphi extra argument for phi
%  is ns x 3 matrix with the boundary sides. Row i contains nodes of sidei.
%  ldir  sides 1 to ldir correspond to Dirichlet boundary condition.
%  TT 3 x 4*nt matrix to plot with command trimesh the quadratic velocity.
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

Utkeep=Ut; tkeep=tiempos; Vtkeep=Vt;
incre=1; % the inrement. Must be 2,4,8, etc. 
         % Results with incre=1 are too much affected by round-off error.
tiempos=tiempos(1:incre:end); Ut=Ut(:,1:incre:end); Vt=Vt(:,1:incre:end);

% Removing the first snapshot due to periodicity.
tiempos=tiempos(2:end); Ut=Ut(:,2:end); Vt=Vt(:,2:end);
dt=mean(diff(tiempos));

Utt=(Ut(:,2:end)-Ut(:,1:end-1))/dt; Vtt=(Vt(:,2:end)-Vt(:,1:end-1))/dt;
tiempos2=(tiempos(2:end)+tiempos(1:end-1))/2;

Uttt=(Utt(:,2:end)-Utt(:,1:end-1))/dt; Vttt=(Vtt(:,2:end)-Vtt(:,1:end-1))/dt;
tiempos3=tiempos(2:end-1);

Utttt=(Uttt(:,2:end)-Uttt(:,1:end-1))/dt; Vtttt=(Vttt(:,2:end)-Vttt(:,1:end-1))/dt;
tiempos4=(tiempos3(2:end)+tiempos3(1:end-1))/2;

Uttttt=(Utttt(:,2:end)-Utttt(:,1:end-1))/dt; Vttttt=(Vtttt(:,2:end)-Vtttt(:,1:end-1))/dt;
tiempos5=tiempos3(2:end-1);

[An,Mn,Dx,Dy,Mp,Sp,GD11,GD12,GD22,J,Se,Me,Xi,Eta]=makronvq_skw(T,z,epsilon,nTc,phi,ccphi);

nUt=sqrt(abs(sum(Ut.*(Mn*Ut)))); nVt=sqrt(abs(sum(Vt.*(Mn*Vt))));
nUtt=sqrt(abs(sum(Utt.*(Mn*Utt)))); nVtt=sqrt(abs(sum(Vtt.*(Mn*Vtt))));
nUttt=sqrt(abs(sum(Uttt.*(Mn*Uttt)))); nVttt=sqrt(abs(sum(Vttt.*(Mn*Vttt))));
nUtttt=sqrt(abs(sum(Utttt.*(Mn*Utttt)))); nVtttt=sqrt(abs(sum(Vtttt.*(Mn*Vtttt))));
nUttttt=sqrt(abs(sum(Uttttt.*(Mn*Uttttt)))); nVttttt=sqrt(abs(sum(Vttttt.*(Mn*Vttttt))));

figure(11); clf
semilogy(tiempos,sqrt(nUt.^2+nVt.^2),'b-','LineWidth',1.5);
semilogy(tiempos2,sqrt(nUtt.^2+nVtt.^2),'c-','LineWidth',1.5);
semilogy(tiempos3,sqrt(nUttt.^2+nVttt.^2),'m-','LineWidth',1.5)
semilogy(tiempos4,sqrt(nUtttt.^2+nVtttt.^2),'r-','LineWidth',1.5)
semilogy(tiempos5,sqrt(nUttttt.^2+nVttttt.^2),'k-','LineWidth',1.5)

Xt=Utkeep; Yt=Vtkeep; temps=tkeep;
incre=2*incre;
temps=temps(1:incre:end); Xt=Xt(:,1:incre:end); Yt=Yt(:,1:incre:end);

% Removing the first snapshot due to periodicity.
temps=temps(2:end); Xt=Xt(:,2:end); Yt=Yt(:,2:end);
dt=mean(diff(temps));


Xtt=(Xt(:,2:end)-Xt(:,1:end-1))/dt; Ytt=(Yt(:,2:end)-Yt(:,1:end-1))/dt;
temps2=(temps(2:end)+temps(1:end-1))/2;

Xttt=(Xtt(:,2:end)-Xtt(:,1:end-1))/dt; Yttt=(Ytt(:,2:end)-Ytt(:,1:end-1))/dt;
temps3=temps(2:end-1);

Xtttt=(Xttt(:,2:end)-Xttt(:,1:end-1))/dt; Ytttt=(Yttt(:,2:end)-Yttt(:,1:end-1))/dt;
temps4=(temps3(2:end)+temps3(1:end-1))/2;

Xttttt=(Xtttt(:,2:end)-Xtttt(:,1:end-1))/dt; Yttttt=(Ytttt(:,2:end)-Ytttt(:,1:end-1))/dt;
temps5=temps3(2:end-1);

a=0.5*(tiempos2(2:2:end-1)+tiempos2(3:2:end)) - temps2;
Eutt=0.5*(Utt(:,2:2:end-1)+Utt(:,3:2:end))-Xtt; eutt=sqrt(abs(sum(Eutt.*(Mn*Eutt))));
Evtt=0.5*(Vtt(:,2:2:end-1)+Vtt(:,3:2:end))-Ytt; evtt=sqrt(abs(sum(Evtt.*(Mn*Evtt))));
nXtt=sqrt(abs(sum(Xtt.*(Mn*Xtt)))); nYtt=sqrt(abs(sum(Ytt.*(Mn*Ytt)))); 
elerr_tt=abs(sqrt(eutt.^2+evtt.^2))./sqrt(nXtt.^2+nYtt.^2);
ee_tt=[norm(a),max(elerr_tt)]

a=tiempos3(3:2:end)-temps3;
Euttt=Uttt(:,3:2:end)-Xttt; euttt=sqrt(abs(sum(Euttt.*(Mn*Euttt))));
Evttt=Vttt(:,3:2:end)-Yttt; evttt=sqrt(abs(sum(Evttt.*(Mn*Evttt))));
nXttt=sqrt(abs(sum(Xttt.*(Mn*Xttt)))); nYttt=sqrt(abs(sum(Yttt.*(Mn*Yttt)))); 
elerr_ttt=abs(sqrt(euttt.^2+evttt.^2))./sqrt(nXttt.^2+nYttt.^2);
ee_ttt=[norm(a),max(elerr_ttt)]

a=0.5*(tiempos4(3:2:end-1)+tiempos4(4:2:end)) - temps4;
Eutttt=0.5*(Utttt(:,3:2:end-1)+Utttt(:,4:2:end))-Xtttt; eutttt=sqrt(abs(sum(Eutttt.*(Mn*Eutttt))));
Evtttt=0.5*(Vtttt(:,3:2:end-1)+Vtttt(:,4:2:end))-Ytttt; evtttt=sqrt(abs(sum(Evtttt.*(Mn*Evtttt))));
nXtttt=sqrt(abs(sum(Xtttt.*(Mn*Xtttt)))); nYtttt=sqrt(abs(sum(Ytttt.*(Mn*Ytttt)))); 
elerr_tttt=abs(sqrt(eutttt.^2+evtttt.^2))./sqrt(nXtttt.^2+nYtttt.^2);
ee_tttt=[norm(a),max(elerr_tttt)]


a=tiempos5(4:2:end-1)-temps5;
Euttttt=Uttttt(:,4:2:end-1)-Xttttt; euttttt=sqrt(abs(sum(Euttttt.*(Mn*Euttttt))));
Evttttt=Vttttt(:,4:2:end-1)-Yttttt; evttttt=sqrt(abs(sum(Evttttt.*(Mn*Evttttt))));
nXttttt=sqrt(abs(sum(Xttttt.*(Mn*Xttttt)))); nYttttt=sqrt(abs(sum(Yttttt.*(Mn*Yttttt)))); 
elerr_ttttt=abs(sqrt(euttttt.^2+evttttt.^2))./sqrt(nXttttt.^2+nYttttt.^2);
ee_ttttt=[norm(a),max(elerr_ttttt)]

format short e, max_norms_derivs=[max(sqrt(nUtt.^2+nVtt.^2));max(sqrt(nUttt.^2+nVttt.^2));...
    max(sqrt(nUtttt.^2+nVtttt.^2));max(sqrt(nUttttt.^2+nVttttt.^2))]', format short

save ../../data/output_data/cylinder_np/deriv_data_L2 max_norms_derivs Utt Vtt Uttt Vttt ...
    Utttt Vtttt Uttttt Vttttt tiempos2 tiempos3 tiempos4 tiempos5