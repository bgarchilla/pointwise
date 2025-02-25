clear all
projections_L2_cycle
fnts=23;
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

figure(51); clf
semilogy(erres,max_err_L2,'b-','LineWidth',1.5); hold on;
semilogy(erres,rhs_of_8(:,1),'c-','LineWidth',1.5);
semilogy(erres,rhs_of_8(:,2),'m-','LineWidth',1.5);
semilogy(erres,rhs_of_8(:,3),'r-','LineWidth',1.5);
semilogy(erres,rhs_of_8(:,4),'k-','LineWidth',1.5);
set(gca,'FontSize',18)
xlabel('$r$','FontSize',fnts,'Interpreter','LaTex')
title('$X=L^2$','FontSize',fnts,'Interpreter','LaTex')
legend('$\max \| (I-P_X^r)u^n\|$','rhs of~(8), $m=2$','rhs of~(8), $m=3$',...
    'rhs of~(8), $m=4$','rhs of~(8), $m=5$','FontSize',fnts,'Interpreter','LaTex')
print -depsc ../../figures/cylinder_np/cyl_np_L2.eps
disp('Generated figure cyl_np_L2.eps. Check figures/cylinder_np')

clear all
 projections_H1_cycle
 fnts=23;

figure(61); clf
semilogy(erres,err_PD1,'b-','LineWidth',1.5); hold on;
semilogy(erres,mus(:,1),'c-','LineWidth',1.5);
semilogy(erres,mus(:,2),'m-','LineWidth',1.5);
semilogy(erres,mus(:,3),'r-','LineWidth',1.5);
semilogy(erres,mus(:,4),'k-','LineWidth',1.5);
set(gca,'FontSize',18)
xlabel('$r$','FontSize',fnts,'Interpreter','LaTex')
title('$X=H^1$','FontSize',fnts,'Interpreter','LaTex')
legend('lhs of (12)','$\mu_2$','$\mu_3$',...
    '$\mu_4$','$\mu_5$','FontSize',fnts,'Interpreter','LaTex')
print -depsc ../../figures/cylinder_np/cyl_np_H1.eps
disp('Generated figure cyl_np_H1.eps. Check figures/cylinder_np')

