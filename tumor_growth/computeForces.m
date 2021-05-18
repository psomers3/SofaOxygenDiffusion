function [F_rep F_adh F] = computeForces1111(cell1,cell2)
%computes the forces between two cells
poisson=0.5;
E=1*10^(3);
alpha_star=2*pi*3.1*10^5;


R_star = cell1.radius*cell2.radius/(cell1.radius+cell2.radius);

if cell1.pos==cell2.pos
    d_ij=0.0001;
else
    d_ij=cell2.pos-cell1.pos;
end
h_ij=cell2.radius+cell1.radius-norm(d_ij);
if h_ij<0
    h_ij=0;
end
F_rep=2/3*(E/(1-poisson^2))*h_ij*R_star^0.5;%4/3*E_star*(R_star^(0.5))*(h_ij^(1.5));

F_adh=alpha_star*(cell1.radius-0.25*h_ij)*h_ij;
F=(F_rep-F_adh)*d_ij/norm(d_ij);
end

