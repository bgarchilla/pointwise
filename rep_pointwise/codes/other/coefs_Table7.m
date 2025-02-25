% The coefficients c_m para j=2,3,m+1, in the estimate
% || f' || <= c_m || f ||^((j-1)/j)|| f^(j) ||^(1/j)
%
%   Writen by Bosco Garcia-Archilla (last modified: January 2025).
%
%   This code comes with no guarantee or warranty of any kind.

m=6; % change number at will


cA=sqrt(2+2/3);
cA=sqrt(2)
cA1=(1+cA^(4/3))^(3/4);
cB1=(2+(2*cA*cA1)^2)^(1/2);

ch=zeros(1,m-1);
c=zeros(1,m);
elprod=cB1; c(2)=elprod;
j=1;
ch(j)=((1+cB1^(2*(j+1)))^((j+1)/(2*(j+2)*j)))*(elprod^((j+1)/((j+2)*j)));
c(3)=c(2)*(ch(1)^(1/(j+1)));
for j=2:m-2
    elprod=elprod*(ch(j-1)^(1/j));
    ch(j)=((1+ch(j-1)^(2*(j+1)))^((j+1)/(2*(j+2)*j)))*(elprod^((j+1)/((j+2)*j)));
    c(j+2)=c(j+1)*(ch(j)^(1/(j+1)));
end
disp('The values in Table 7 (for m=2,3,4,5,6) are the following ones:')
format short e, cm=c(2:m), format short