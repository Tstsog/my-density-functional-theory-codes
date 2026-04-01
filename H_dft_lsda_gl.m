% This matlab code computes the ground state energy of hydrogen (H) atom by
% solving the Kohn-Sham (KS) equation in the local-spin-density-approximation
% (LSDA) of the density-functional theory (dft) [1]. The exchange and correlation contribution (E_xc) is computed in the Gunnarsson-Lundqvist (GL) approximation [2]. 
%
% The Hamiltonian is approximated with the Legendre pseudo-spectral method
% [3]. The atomic unit (au) is used in the calculation.
%
% Refs: 
% [1] R.G. Parr and W. Yang, Density-Functional Theory of Atoms and Molecules (Oxford University Press, Oxford, 1989);
% [2] O. Gunnarsson and B. I. Lundqvist, ``Exchange and correlation in atoms, molecules, and solids by the spin-density-functional formalism", Phys. Rev. B \textbf{13}, 4274 (1976)
% [3] Ts Tsogbayar and M Horbatsch 2013 J. Phys. B: At. Mol. Opt. Phys. 46 085004; Tsogbayar Tsednee and Danny L Yeager 2017 Chinese Phys. B 26 083101;  
% https://yorkspace.library.yorku.ca/xmlui/bitstream/handle/10315/28160/Tsednee_Tsogbayar_2014_PhD.pdf?sequence=2
% 
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% March 31, 2026 & Institute of Physics and Technology, Mongolian Academy of Sciences
%
function [] = H_dft_lsda_gl
%
clc;
clear H_dft_lsda_gl
format short 
%
Z = 1.;   % charge of magnesium nucleus; 
N = 2*128;   % number of grid point; you may change it
rb = 0.0; % beginning point of coordinate r 
re = 30.; % ending point of coordnate r; you may change it
% unit used here is atomic unit
%%%
itermax = 500; tol = 10^(-8);
alf = 0.75; 
%
a = rb;
b = re;
%
[r,wr,D] = legDC2(N,a,b);
wr = wr(2:N);
r = r(2:N);
D2 = (2/(b-a))^2*D^2;
D2 = D2(2:N,2:N);
%
%[En_1, V1, rho_1s] = H0_eig_val(D2,wr,r,Z,eig_ind,ell)
[En_1s, V1s, rho_1s] = H0_eig_val(D2,wr,r,Z,1,0.);
[En_2s, V2s, rho_2s] = H0_eig_val(D2,wr,r,Z,2,0.);
%
%
gamma = 0.297;
Ry = 0.5; % Ry = 0.5au
alpha = (4/(9*pi))^(1/3);
c_p = 0.0666;
c_f = 0.0406;
r_p = 11.4; 
r_f = 15.9;
%
rho_alpha = (rho_1s ); 
rho_total = rho_alpha; 
rho_1s_old = rho_total;
%
for iter = 1:itermax
    %
    iter
    rho_1s = rho_1s_old;
    xi = (rho_alpha )./rho_1s;    
    %
    f_rhs_11 = -4.*pi.*r.*(rho_1s); %
%
    u_11 = D2\f_rhs_11;
    u_11 = u_11 + Z.*r/(b-a); % for 1s state or s-states
    %u_11 = [0; u_11; 0];
    Vpois_11 = u_11./r;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% O. GUNNARSSON AND B. I. LUNDQVIST 
    rs = (3./(4.*pi.*rho_1s)).^(1/3);
    x = rs/11.4;
    %
    beta_rs = 1. + 0.0545.*rs.*log(1. + (1./x));
    delta_rs = 1 - 0.036.*rs - (1.36.*rs)./(1. + 10.*rs);
    mu_x_p = -2./(pi.*alpha.*rs);
    V_xc_alpha = mu_x_p.* (beta_rs + (1/3).*delta_rs.*xi./(1 + gamma.*xi)).*Ry;
    V_xc_beta = mu_x_p.* (beta_rs - (1/3).*delta_rs.*xi./(1 - gamma.*xi)).*Ry;
%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %
    [En_1s_alpha, V1_1s, rho_1s_alpha] = H_KS_eig_val(D2,wr,r,Z, 1.,0., Vpois_11, V_xc_alpha);
    [En_1s_beta, V1_1s, rho_1s_beta] = H_KS_eig_val(D2,wr,r,Z, 1.,0., Vpois_11, V_xc_beta);
    %
    rho_alpha = rho_1s_alpha; 
    rho_beta  = rho_1s_beta; 
    rho_total = rho_alpha + 0.*rho_beta;
    rho_new = rho_total;
    rho_n = alf.*rho_new + (1.-alf).*rho_1s;
    %
    if (abs(rho_n - rho_1s_old) < tol)
        break
    end
    %
    rho_1s_old = rho_n; 
%  
end
%
[En_1s_alpha, En_1s_beta] * 27.211386;
%
%%%
x_p = rs./r_p;
x_f = rs./r_f;
eps_x_p = -(3./(2.*pi.*alpha.*rs)).*Ry;
eps_x_f = 2^(1/3).*eps_x_p;
%
eps_xc_p = eps_x_p - c_p.*( (1.+x_p.^3).*log(1.+(1./x_p)) + 0.5.*x_p - x_p.^2 - (1/3)).*Ry;
eps_xc_f = eps_x_f - c_f.*( (1.+x_f.^3).*log(1.+(1./x_f)) + 0.5.*x_f - x_f.^2 - (1/3)).*Ry;
f_xi = ( (1.+xi).^(4/3) + (1.-xi).^(4/3) - 2 )./(2^(4/3)-2);
%
eps_xc = eps_xc_p + (eps_xc_f - eps_xc_p).*f_xi;      %
%
V_xc_alpha = sum(wr.*rho_alpha.*V_xc_alpha.*r.*r.*4.*pi);
V_xc_beta = sum(wr.*rho_beta.*V_xc_beta.*r.*r.*4.*pi) ;
V_xc = V_xc_alpha + 0.*V_xc_beta;
%
Vee = sum(wr.*rho_1s_old.*Vpois_11.*r.*r.*4.*pi);
En_xc = sum(wr.*rho_1s_old.*eps_xc.*r.*r.*4.*pi); 
%%%
En_total_gr_au = En_1s_alpha - 0.5*Vee + (En_xc - V_xc); % 
%
[En_1s_alpha, real(En_total_gr_au) ] % -0.2834   -0.4920
% [-0.2834   -0.4920]

%%%
return
end

%%%
function [En_1, V1, rho_1s] = H0_eig_val(D2,wr,r,Z,eig_ind,ell)
%
%%%
H_ham = -0.5.*D2 + diag(ell*(ell+1)./(2.*r.^2)) - diag(Z./r);
[Vec,En] = eig(H_ham);                                     % Eigenvalue problem
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
%[En(1),En(2),En(3),En(4),En(5)];
%           
Vec = Vec(:,ij);                       % The unnormalized eigenfunctions
V1 = Vec(:,eig_ind);                         % The unnormalized eigenfunction for the ground state,
% here 1 is for the ground state, 2,3... are corresponded
% for the excited states
%V1 = [0.,;V1,;0.];
n_c = sum(wr.*conj(V1).*V1.*4.*pi);
V1 = 1./sqrt(n_c).*V1;
R1s = V1./r;
rho_1s = R1s.*R1s;  

%plot(r,V1)
En_1 = En(eig_ind);
%En_2 = En(2);
%En_3 = En(3);

%%%
return
end
%%%
%%%
function [En_1, V1, rho_1s] = H_KS_eig_val(D2,wr,r,Z,eig_ind,ell,Vpois_11,  V_xc)
%
%%%
H_ham = -0.5.*D2 + diag(ell*(ell+1)./(2.*r.^2)) - diag(Z./r) + diag(Vpois_11) + diag( V_xc);
[Vec,En] = eig(H_ham);                                     % Eigenvalue problem
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
%[En(1),En(2),En(3),En(4),En(5)];
%           
Vec = Vec(:,ij);                       % The unnormalized eigenfunctions
V1 = Vec(:,eig_ind);                         % The unnormalized eigenfunction for the ground state,
% here 1 is for the ground state, 2,3... are corresponded
% for the excited states
%V1 = [0.,;V1,;0.];
n_c = sum(wr.*conj(V1).*V1.*4.*pi);
V1 = 1./sqrt(n_c).*V1;
R1s = V1./r;
rho_1s = R1s.*R1s; % times 2   !!!!!!!

%plot(r,V1)
En_1 = En(eig_ind);
%En_2 = En(2);
%En_3 = En(3);
%%%
return
end
%%%
%%%
        function [xi,w,D]=legDC2(N,a,b)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % legDc.m
            %
            % Computes the Legendre differentiation matrix with collocation at the
            % Legendre-Gauss-Lobatto nodes.
            %
            % Reference:
            %   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
            %   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
            %
            % Written by Greg von Winckel - 05/26/2004
            % Contact: gregvw@chtm.unm.edu
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % Truncation + 1
            N1=N+1;
            
            % CGL nodes
            xc=cos(pi*(0:N)/N)';
            
            % Uniform nodes
            xu=linspace(-1,1,N1)';
            
            % Make a close first guess to reduce iterations
            if N<3
                x=xc;
            else
                x=xc+sin(pi*xu)./(4*N);
            end
            
            % The Legendre Vandermonde Matrix
            P=zeros(N1,N1);
            
            % Compute P_(N) using the recursion relation
            % Compute its first and second derivatives and
            % update x using the Newton-Raphson method.
            
            xold=2;
            while max(abs(x-xold))>eps
                
                xold=x;
                
                P(:,1)=1;    P(:,2)=x;
                
                for k=2:N
                    P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
                end
                
                x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
            end
            
            X=repmat(x,1,N1);
            Xdiff=X-X'+eye(N1);
            
            L=repmat(P(:,N1),1,N1);
            L(1:(N1+1):N1*N1)=1;
            D=(L./(Xdiff.*L'));
            D(1:(N1+1):N1*N1)=0;
            D(1)=(N1*N)/4;
            D(N1*N1)=-(N1*N)/4;
            
            % Linear map from[-1,1] to [a,b]
            xi=(a*(1-x)+b*(1+x))/2;        % added by Tsogbayar Tsednee
            
            % Compute the weights
            w=(b-a)./(N*N1*P(:,N1).^2);    % added by Tsogbayar Tsednee
            
        end

