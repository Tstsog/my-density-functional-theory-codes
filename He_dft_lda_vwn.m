% This matlab code computes the ground state energy of helium (He) atom by
% solving the Kohn-Sham (KS) equation in the local-density-approximation
% (LDA)  of the density-functional theory (dft) [1]. The exchange and correlation contribution (E_xc) is computed in the Vosko, Wilk, and
% Nusair (VWN) approximation [2]. The obtained energy value is compared
% with available data [3]. 
%
% The Hamiltonian is approximated with the Legendre pseudo-spectral method
% [4]. The atomic unit (au) is used in the calculation.
%
% Refs: 
% [1] R.G. Parr and W. Yang, Density-Functional Theory of Atoms and Molecules (Oxford University Press, Oxford, 1989);
% [2] S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980); S.H. Vosko and L. Wilk, Phys. Rev. B 22, 3812 (1980).;
% [3] https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations.
% [4] Ts Tsogbayar and M Horbatsch 2013 J. Phys. B: At. Mol. Opt. Phys. 46 085004; Tsogbayar Tsednee and Danny L Yeager 2017 Chinese Phys. B 26 083101;  
% https://yorkspace.library.yorku.ca/xmlui/bitstream/handle/10315/28160/Tsednee_Tsogbayar_2014_PhD.pdf?sequence=2
% 
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% July 24, 2023 & University of North Dakota 
%
function [] = He_dft_lda_vwn
%
clear; clc;
format long
%
clear; clc;
Z = 2.;   % charge of helium nucleus; 
N = 4*64;   % number of grid point; you may change it
rb = 0.0; % beginning point of coordinate r 
re = 30.; % ending point of coordnate r; you may change it
% unit used here is atomic unit
%%%
itermax = 100; tol = 10^(-8);
alf = 0.95; 
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
[En_1s, rho_1s] = H0_eig_val(D2,wr,r,Z,1,1,0); % ns = 1 & 1s orbital
[En_2s, rho_2s] = H0_eig_val(D2,wr,r,Z,2,2,0); % ns = 2 & 2s orbital
%
[En_1s, En_2s];
%
rho_old = 2.*rho_1s + 0.*rho_2s; 
%
x_alpha = 3./4;
%
ac1 = 0.0621814/2;
bc1 = 3.72744;
cc1 = 12.9352;
x0 = -0.10498;
%
n_xc = 1.; % if n_ex = 0, X-lapha; n_ex = 1, LDA-VWN;
%
for iter = 1:itermax
    %
    iter
    rho_ = rho_old;
    %
    f_rhs_11 = -4.*pi.*r.*(rho_./Z); %
%
    u_11 = D2\f_rhs_11;
    u_11 = u_11 + r./(b-a); %  for spherical symmetry
    %u_11 = [0; u_11; 0];
    Vpois_11 = u_11./r ;
    %
    V_ex = -(3./pi).^(1/3).*rho_.^(1/3);
    eps_e = x_alpha.*V_ex;
    %%%
    rs = (3./(4.*pi.*rho_)).^(1/3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Vosko-Wilk-Nusair correlation & VWM5
    x = rs.^(1/2);
    Xx = x.^2 + bc1.*x + cc1;
    Xx0 = x0^2 + bc1.*x0 + cc1;
    Q = sqrt(4.*cc1 - bc1^2);
    %
    t2 = Q./(2.*x + bc1);
    t1 = log((x-x0).^2./(Xx)) + (2.*(bc1+2.*x0)./Q).*atan(t2);
    t = log(x.^2./Xx) + (2.*bc1./Q).*atan(t2) - (bc1.*x0.*t1)./Xx0;
    eps_c = ac1.*t;
    V_corr = eps_c - (ac1./3).* (cc1.*(x - x0) - bc1.*x.*x0)./((x-x0).*(rs+bc1.*x+cc1));        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %[En, rho_ks]    = H_KS_eig_val(D2,wr,r,Z, ns, wf_ind, ell, Vpois_11,V_ex, V_corr)
    [En_1s, rho_1s_] = H_KS_eig_val(D2,wr,r,Z, 1., 1., 0., Vpois_11,V_ex, n_xc.*V_corr); % KS: 1s orbital
    [En_2s, rho_2s_] = H_KS_eig_val(D2,wr,r,Z, 2., 2., 0., Vpois_11,V_ex, n_xc.*V_corr); % KS: 2s orbital    
    %
    rho_new = 2.*rho_1s_ + 0.*rho_2s_;
    %
    rho_n = alf.*rho_new + (1.-alf).*rho_;
    %
    if (abs(rho_n - rho_old) < tol)
        break
    end
    %
    rho_old = rho_n; 

end
%%%
[En_1s, En_2s] ;% 1s and 2s orbital energies
%
Vee = sum(wr.*rho_old.*Vpois_11.*r.*r.*4.*pi);
%
E_x = sum(wr.*rho_old.*eps_e.*r.*r.*4.*pi) ;
V_x = sum(wr.*rho_old.*V_ex.*r.*r.*4.*pi) ;
V_c = sum(wr.*rho_old.*V_corr.*r.*r.*4.*pi) ;
E_c = sum(wr.*rho_old.*eps_c.*r.*r.*4.*pi);
%
En_xc = (E_x + n_xc.*E_c ) - (V_x + n_xc.*V_c);
%
En_total_gr_au = (2.*En_1s) - 0.5*Z.*Vee + En_xc ;
%
[En_1s, En_total_gr_au] 
% [-0.570424722963091  -2.834835624055025]  vs 
% [-0.570425,          -2.834836] from https://www.nist.gov/pml/atomic-reference-data-electronic-structure-calculations/atomic-reference-data-electronic-7-1
 
%%%

plot(r, rho_n.*4*pi.*r.^2, 'b')
ylabel('4\pi r^{2}\rho')
xlabel('r (au)')

%%%
return
end

%%%
function [En, rho] = H0_eig_val(D2,wr,r,Z,ns,wf_ind,ell)
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
V1 = Vec(:,wf_ind);                         % The unnormalized eigenfunction for the ground state,
%V1 = [0.,;V1,;0.];
n_c = sum(wr.*V1.*V1.*4.*pi);
V1 = 1./sqrt(n_c).*V1;
R = V1./r;
rho = R.*R; 
%
En = En(ns);
%%%
return
end
%%%
%%%
function [En, rho_ks] = H_KS_eig_val(D2,wr,r,Z, ns, wf_ind, ell, Vpois_11,V_ex, V_corr)
%
%%%
H_ham = -0.5.*D2 + diag(ell*(ell+1)./(2.*r.^2)) - diag(Z./r) + Z.*diag(Vpois_11) + diag(V_ex) + diag(V_corr);
[Vec,En] = eig(H_ham);                                     % Eigenvalue problem
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
%[En(1),En(2),En(3),En(4),En(5)];
%           
Vec = Vec(:,ij);                       % The unnormalized eigenfunctions
V1 = Vec(:,wf_ind);                         % The unnormalized eigenfunction for the ground state,
%V1 = [0.,;V1,;0.];
n_c = sum(wr.*V1.*V1.*4.*pi);
V1 = 1./sqrt(n_c).*V1;
R = V1./r;
rho_ks = R.*R; 
%
En = En(ns);
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

