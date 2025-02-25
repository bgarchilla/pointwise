function [x,y,xix,xiy,etx,ety,J,xxi,xet,yxi,yet]=phia(xi,eta,z1,z2,c)
% z1 and z2 nIc x 3
% xi and eta 1 x nxi
%r=sqrt((z1-c(:,1)*ones(size(z1))).^2+(z2-c(:,2)*ones(size(z2))).^2);
% the values of the basis functions and their derivatives on xi and eta
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

cx=xi; cy=eta;
Nc=[(1-cx-cy).*(1-2*cx-2*cy); cx.*(2*cx-1); cy.*(2*cy-1); ...
         4*cx.*(1-cx-cy); 4*cx.*cy; 4*cy.*(1-cx-cy)];
Ncxi=[-(1-2*cx-2*cy)-2*(1-cx-cy); 4*cx-1; zeros(size(cx));...
          4*(1-cx-cy)-4*cx; 4*cy; -4*cy];
Ncet=[-(1-2*cx-2*cy)-2*(1-cx-cy); zeros(size(cx)); 4*cy-1;...
          -4*cx; 4*cx; 4*(1-cx-cy)-4*cy];
clear cx cy

r=sqrt((z1-c(:,1)).^2+(z2-c(:,2)).^2);
I=[2 3];
% I=find(r>1-10000*eps);
% if length(I)<2; I=[]; end
% if ~isempty(I)
%     j=find(r<=1-10000*eps);
%     k=I(1); m=I(2);
%     E=[1 -1 -1; 0 1 0; 0 0 1]; % the coordinates of the three vertices;
%     D=[1;0]*E(k,:) + [0; 1]*E(m,:);
% else
     j=1; k=2; m=3;
% end
    x=z1*Nc;
    y=z2*Nc;
% if isempty(I)
%     J=(z1(2)-z1(1))*(z2(3)-z2(1)) - (z1(3)-z1(1))*(z2(2)-z2(1));
% %     xxi=(z1(k)-z1(j))*ones(size(xi)); xet=(z1(m)-z1(j))*ones(size(xi));
% %     yxi=(z2(k)-z2(j))*ones(size(xi)); yet=(z2(m)-z2(j))*ones(size(xi));
%     xix=(z2(3)-z2(1))/J; xiy=-(z1(3)-z1(1))/J;
%     etx=-(z2(2)-z2(1))/J; ety=(z1(2)-z1(1))/J;
%     J=J*ones(size(xi));
%     xix=xix*ones(size(xi)); xiy=xiy*ones(size(xi)); 
%     etx=etx*ones(size(xi)); ety=ety*ones(size(xi));
% else
    xxi=z1*Ncxi; xet=z1*Ncet;
    yxi=z2*Ncxi; yet=z2*Ncet;
%     % transform of xi and eta so that the arc goes between vertices 2 and 3
%     xiorig=xi; etaorig=eta;
%     xi = D(1,1) + D(1,2:3)*[xiorig; etaorig]; eta = D(2,1) + D(2,2:3)*[xiorig; etaorig];
     J=xxi.*yet - xet.*yxi;
     xix=yet./J; xiy=-xet./J; etx=-yxi./J; ety=xxi./J;
end
