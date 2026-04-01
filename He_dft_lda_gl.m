% This matlab code computes the ground state energy of helium (He) atom by
% solving the Kohn-Sham (KS) equation in the local-density-approximation
% (LDA) of the density-functional theory (dft) [1]. The exchange and correlation contribution (E_xc) is computed in the Gunnarsson-Lundqvist (GL) approximation [2]. 
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
function [] = He_dft_lda_gl
%
clc;
clear He_dft_lda_gl
format short 
%
Z = 2.;   % charge of Helium nucleus; 
N = 4*64;   % number of grid point; you may change it
rb = 0.0; % beginning point of coordinate r 
re = 30.; % ending point of coordnate r; you may change it
% unit used here is atomic unit
%%%
itermax = 100; tol = 10^(-8);
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
[En_1, En_2, En_3, V1, rho_1s] = H0_eig_val(D2,wr,r,Z,1,0.);
[En_1, En_2, En_3, V2, rho_2s] = H0_eig_val(D2,wr,r,Z,2,0.);
%
[En_1, En_2,En_3]
%
rho_1s_old = 2.*rho_1s + 0.*rho_2s; 
%
x_alpha = 3./4;
%
n_xc = 1.;  % if n_ex = 0, X-lapha; n_ex = 1, LDA-GL;
%
for iter = 1:itermax
    %
    iter
    rho_1s = rho_1s_old;
    %
    f_rhs_11 = -4.*pi.*r.*(rho_1s); %
%
    u_11 = D2\f_rhs_11;
    u_11 = u_11 + Z.*r/(b-a); % for 1s state or s-states
    %u_11 = [0; u_11; 0];
    Vpois_11 = u_11./r;
    %
    V_ex = -(3./pi).^(1/3).*rho_1s.^(1/3);
    eps_e = x_alpha.*V_ex;
    %%%
    rc = (3./(4.*pi.*rho_1s)).^(1/3);
    xx = rc./11.4;
    eps_c = -0.0666.*( (1.+xx.^3).*log(1.+(1./xx)) + 0.5.*xx - xx.^2 - (1./3) ) *(1/2);
    V_corr = -0.0666.*log(1. + (1./xx)) *(1/2);
    %%%          
    %
    [En_1_1s, En_2_1s, En_3_1s, V1_1s, rho_1s_] = H_KS_eig_val(D2,wr,r,Z, 1.,0., Vpois_11,V_ex, n_xc.*V_corr); % 1s orbital 
    [En_1_2s, En_2_2s, En_3_2s, V1_2s, rho_2s_] = H_KS_eig_val(D2,wr,r,Z, 2.,0., Vpois_11,V_ex, n_xc.*V_corr); % 2s orbital
    [En_1_2p, En_2_2p, En_3_2p, V1_2p, rho_2p_] = H_KS_eig_val(D2,wr,r,Z, 1.,1., Vpois_11,V_ex, n_xc.*V_corr); % 2p orbital - virtial    
    %
    rho_new = 2.*rho_1s_ + 0.*rho_2s_ + 0.*rho_2p_;
    %
    rho_n = alf.*rho_new + (1.-alf).*rho_1s;
    %
    if (abs(rho_n - rho_1s_old) < tol)
        break
    end
    %
    rho_1s_old = rho_n; 

end
%%%

[En_1_1s, En_2_2s, En_1_2p]; % % orbital energy eigenvalues 

%
%
Vee = sum(wr.*rho_1s_old.*Vpois_11.*r.*r.*4.*pi);
%
E_x = sum(wr.*rho_1s_old.*eps_e.*r.*r.*4.*pi) ;
V_x = sum(wr.*rho_1s_old.*V_ex.*r.*r.*4.*pi) ;
V_c = sum(wr.*rho_1s_old.*V_corr.*r.*r.*4.*pi) ;
E_c = sum(wr.*rho_1s_old.*eps_c.*r.*r.*4.*pi);
%
En_xc = (E_x + n_xc.*E_c ) - (V_x + n_xc.*V_c);
%
En_total_gr_au = (2.*En_1_1s) - 0.5*Vee + En_xc ;
%
[En_1_1s, En_total_gr_au] % -0.2416   -0.4537


%%%
return
end

%%%
function [En_1, En_2, En_3, V1, rho_1s] = H0_eig_val(D2,wr,r,Z,eig_ind,ell)
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
n_c = sum(wr.*V1.*V1.*4.*pi);
V1 = 1./sqrt(n_c).*V1;
R1s = V1./r;
rho_1s = R1s.*R1s; % times 2   !!!!!!!

%plot(r,V1)
En_1 = En(1);
En_2 = En(2);
En_3 = En(3);

%%%
return
end
%%%
%%%
function [En_1, En_2, En_3, V1, rho_1s] = H_KS_eig_val(D2,wr,r,Z,eig_ind,ell,Vpois_11,V_ex,V_corr)
%
%%%
H_ham = -0.5.*D2 + diag(ell*(ell+1)./(2.*r.^2)) - diag(Z./r) + diag(Vpois_11) + 1.*diag(V_ex) + 1.*diag(V_corr);
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
n_c = sum(wr.*V1.*V1.*4.*pi);
V1 = 1./sqrt(n_c).*V1;
R1s = V1./r;
rho_1s = R1s.*R1s; % times 2   !!!!!!!

%plot(r,V1)
En_1 = En(1);
En_2 = En(2);
En_3 = En(3);
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

