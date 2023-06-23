function HED_radial_EzHz_TGRS

% does the radial sweep in the near field, for fig 6.
clc
clear, close all
tic,
%% initializations
f= 1e6; % the operating frequency
w= 2*pi*f; 
lambda= 3e8/f; % operating wavelength
e0= 8.854e-12; % permittivity of free space
u0= pi*4e-7; % permeability of free space

e1= 3.2; % the dielectric constant of the substance, assuming glacier
% e1= 10, assuming  ground  
epsr= [ 1 e1 1]*e0; % the dielectric constants of the three media 
mu= [ 1 1 1]*u0;
LT1= 5; % the loss tangent of the substance, also try for 50, 500
losst= [ 0 LT1 5000]; % losst defined at 1 MHz
sigma= losst.*epsr*w;

I = 1; % the driving current
len=1; % length of the dipole
z=2; % receiver height = 2 m
d_array = [1, 2.5, 5]; % the thickness of the substance under test

for dd= 1: length(d_array)
d= d_array(dd); 
    disp(['calculating radial sweep for thickness = ' num2str(d) 'm']) 
    k= w*sqrt(u0.*e0); % the propagation constant in  the  source (air) medim
    
    constETM= 1i*I*len/(4*pi*w*e0);
    constHTM= I*len/(4*pi);
    constETE= I*len*w*u0/(4*pi);
    constHTE= 1i*I*len/(4*pi);
    
    count=1; %number of data sample marker.
    
    delrho= 0.2; % the 
    rho_init=2; % The starting point of the radial sweep
    rho_end=20; % the end point of the radial  sweep
    
    for rho= (rho_init:delrho:rho_end)
        %% definition of the integral path (Shifat) Along the sommerfeld path , with krhoi and delkrho as function of rho;
        delkrho= 0.05/rho; %seems to converge well
        tol= 1e-6; % tolerance of convergence
        
        %% Determining Hz
        krho= -.01i/rho;
        sumIntegrandHz= 0;
        err= Inf;
        flag=3;
        datapoint=0;
        while( err >tol)        
            temp= sumIntegrandHz;
            %term 1
            kz= (sqrt(k.^2- krho^2));
            J= besselj(1, krho*rho);
            [~, R_TE]= planewave_response_TGRS(f,epsr/e0, mu/u0,sigma,[0 d], krho);            
            term1= Hz(krho,kz(1),J,R_TE,z);
            %term 2
            kz= sqrt(k.^2- (krho+delkrho)^2);
            J= besselj(1, (krho+delkrho)*rho);
            [~, R_TE]= planewave_response_TGRS(f,epsr/e0, mu/u0,sigma,[0 d], krho+delkrho);
            term2= Hz(krho+delkrho,kz(1),J,R_TE,z);
            %passing in simpsons routine
            [sum, flag]= simpsons_routine(datapoint,term1, term2, flag);
            sumIntegrandHz= sumIntegrandHz+sum;
            err= abs((temp-sumIntegrandHz)/sumIntegrandHz);
            krho= krho+delkrho;
            datapoint= datapoint+1;
        end
        Hz_simpsons(count)=constHTE*sumIntegrandHz*delkrho/3;    
        %% Determining Ez
            krho= -.01i/rho;
            sumIntegrandEz= 0;
            err= Inf;
            flag=3;
            datapoint=0;
            while( err >tol)
                temp= sumIntegrandEz;
                %term 1
                kz= sqrt(k.^2- krho^2);
                kz= (sqrt(k.^2- krho^2));
                J= besselj(1, krho*rho);
                [R_TM,~]=planewave_response_TGRS(f,epsr/e0, mu/u0,sigma,[0 d], krho);
                term1= Ez(krho,kz(1),J,R_TM,z);     
                %term 2
                kz= sqrt(k.^2- (krho+delkrho)^2);
                J= besselj(1, (krho+delkrho)*rho);
                [R_TM,~]=planewave_response_TGRS(f,epsr/e0, mu/u0,sigma,[0 d], krho+delkrho);
                term2=  Ez(krho+delkrho,kz(1),J,R_TM,z);      
                %passing in simpsons routine
                [sum, flag]= simpsons_routine(datapoint,term1, term2, flag);
                sumIntegrandEz= sumIntegrandEz+sum;
                err= abs((temp-sumIntegrandEz)/sumIntegrandEz);
                krho= krho+delkrho;
                datapoint= datapoint+1;
            end
            Ez_simpsons(count)=constETM*sumIntegrandEz*delkrho/3;       
            count= count+1 ;
    end
    
    figure(1)
    semilogy(rho_init:delrho:rho_end,abs(Hz_simpsons),'linewidth', 2);
    title('Hz Simpsons')
    grid on
    hold on
    figure(2)
    semilogy(rho_init:delrho:rho_end,abs(Ez_simpsons),'linewidth', 2);
    title('Ez Simpsons')
    grid on
    hold on
end
legend([num2str(transpose(d_array))]);
disp (['completed. time elapsed= ' num2str(toc) 's'])
end

function y=Hz (krho,kz,J,R_TE,z)
y=(krho^2/kz)*(1+R_TE)*J*exp(1i*kz*z); % equation 3d
end

function y= Ez(krho,kz,J,R_TM,z) % equation 3a, not presented in  the paper
y= krho^2 *(1-R_TM)*J* exp(1i*kz*z);
end

function [sum, flag]= simpsons_routine(pointval,term1, term2, flag)
if (pointval==0);
    sum= term1+term2;
else
    sum= flag*term1 +term2;
    if flag==3
        flag=1;
    elseif flag==1
        flag=3;
    else
        disp ('Wrong flagging')
    end
end
end

