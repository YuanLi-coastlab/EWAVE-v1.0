%  calculate the bed slope based on Long et al. 2006

function [p] = bed_slope_functions(be,dx)
% nx=length(be);
% for ix=3:nx-2;
%      p(ix)=(4/3)*(be(ix+1,2)-be(ix-1,2))./(2*dx)-(1/3)*(be(ix+2,2)-be(ix-2,2))./(4*dx);
% end
% p(1)=(4/3)*(be(2,2)-be(1,2))./dx-(1/3)*(be(3,2)-be(1,2))./(2*dx);
% p(2)=(4/3)*(be(3,2)-be(1,2))./(2*dx)-(1/3)*(be(4,2)-be(1,2))./(3*dx);
% p(nx-1)=(4/3)*(be(nx,2)-be(nx-2,2))./(2*dx)-(1/3)*(be(nx,2)-be(nx-3,2))./(3*dx);
% p(nx)=(4/3)*(be(nx,2)-be(nx-1,2))./(dx)-(1/3)*(be(nx,2)-be(nx-2,2))./(2*dx);

nx=length(be);
for ix=2:nx-1
     p(ix)=(be(ix+1,2)-be(ix-1,2))/(2*dx);
end
p(1)=(be(2,2)-be(1,2))/(dx);
p(nx)=(be(nx,2)-be(nx-1,2))/(dx);






