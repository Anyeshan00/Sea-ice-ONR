function [ R_TM, R_TE] = planewave_response_TGRS(f,eps, mu,sigma, d, kp)

%% initalize variables
w= 2*pi*f;
e0= 8.854e-12;
u0= pi*4e-7;
eps= eps*e0;
mu= mu*u0;
len= length(eps);

%% defining k and kz, in compatibility with the main function
%Ref[10] : Chew, Waves and Fields in Inhomogeneous Media.

%%%%%%%%% 
% NOTE: THE MANUSCRIPT IS LEFT WITH A TYPO,   
% ### WE HAVE USED epsilon_complex = epslion + 1i*sigma/w; IN THE MANUSCRIPT IT IS
% WRITTEN epsilon_complex = epslion - 1i*sigma/w; PLEASE CORRECT EQ. 8(a),

eps_c = eps+ 1i*sigma(1:end)/w;
k2= w*sqrt(mu.*eps_c);
kz= sqrt(k2.^2- kp^2);

%% determining R_TM and R_TE
for i= len-1:-1:2
    h(i)= d(i-1)-d(i); %in this convention,
    % h is always a negative quantity in the code and
    % in the equations (8) and (18) of the manuscript.
end

% using formulations of Ref[10], equation (2.1.13),(2.1.14),(2.1.24)
% determing R_TM and  R_TE using recursive computation bottom up
R= (eps_c(len)*kz(len-1) - eps_c(len-1)*kz(len)) / (eps_c(len)*kz(len-1)+ eps_c(len-1)*kz(len));
for i=len-1: -1 : 2
    r= (eps_c(i)*kz(i-1)- eps_c(i-1)*kz(i) )/ (eps_c(i)*kz(i-1)+ eps_c(i-1)*kz(i));
    num= r+ R*exp (-2i*kz(i)*h(i));
    den= 1+ r*R*exp(-2i*kz(i)*h(i));
    R= num/den;
end
R_TM= R; % equation 8a

R= (mu(len)*kz(len-1) - mu(len-1)*kz(len) )/ (mu(len)*kz(len-1)+ mu(len-1)*kz(len));
for i=len-1: -1 : 2
    r= (mu(i)*kz(i-1)- mu(i-1)*kz(i) )/ (mu(i)*kz(i-1)+ mu(i-1)*kz(i));
    num= r+ R*exp (-2i*kz(i)*h(i));
    den= 1+ r*R*exp(-2i*kz(i)*h(i));
    R= num/den;
end
R_TE= R; % equation 8b


end

