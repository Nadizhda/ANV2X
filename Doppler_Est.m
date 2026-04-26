% Function for the calculation of the Doppler shift based on the reference
% paper :J. Cai, W. Song and Z. Li, “Doppler Spread Estimation for Mobile OFDM Systems in Rayleigh Fading Channels,”
% in IEEE Transactions on Consumer Electronics, vol. 49, no. 4, pp. 973-977, Nov. 2003 
% m is the maximum delay spread
% Ts is the symbol duration without guard band
% noise_var is the noise power 
% N is the number of subcarriers(number of OFDM symbols, N_IFFT size)
% mode_symb is the input data, and it can be the modulated symbols (BPSK, QPSK, 16QAM).

 function Delta_f = Doppler_Est(N,noise_var,Ts,mode_symb,tau_max,f,L_c)
% N=512;
% mode_symb= Chann_signal(1:N+N/4);


% L_c= N/4; % the duration of the gurad time is L*Ts = 4.6  micro second ceil(L_c);


T_ds = Ts/(N+L_c); % useful symbol 

M =  ceil(tau_max/T_ds(1)) ;

[s,m]= size(mode_symb);

% Considering all the symbols are experinceing the same channel in the
% considered time. 
if s~=1
    non_demsi = reshape(mode_symb, 1,[]);
else
    non_demsi =mode_symb;
end

N_nonzer=[];
F_rho=[];
E_rho=[];
Rh_re=[];
Rh_re=[];
c = 3e8;
A_r=[];

A_r=[];
% D =Doppler_Est(N,noise_var,Ts,non_demsi,tau_max);
for i=1:s

%     F_rho=[];
%     for j =(N - L):N
        rho_re=[];
        alpha_re=[];
        rho_im=[];
        alpha_im=[];
        l = zeros(1,M(1));l1=zeros(1,M(1));
        for k=1:M
            l(k)= floor(((i-1)*(N+L_c- M)+ k));
            l1(k) = floor((((i-1)*(N+L_c-M))+ (k+N)));
            rho_re(k) = real(mode_symb(l(k))) * real(mode_symb(l1(k)));
            rho_im(k) = imag(mode_symb(l(k))) * imag(mode_symb(l1(k)));
            alpha_re(k) = (real(mode_symb(l(k))) )^2;
            alpha_im(k) = (imag(mode_symb(l1(k))))^2 ;
        end
%         l,l1
        Rh_re(i) = (1/M) * sum(rho_re);
        Rh_im(i) = (1/M) * sum(rho_im);
        A_r(i) = (1/M) *sum(alpha_re);
        A_im(i) = (1/M) *sum(alpha_im);
  
        F_rho(i) = 0.5*((Rh_re(i)/A_r(i))+ (Rh_im(i)/A_im(i)));
%     end
    N_nonzer(i) = (nonzeros(F_rho(i)));
    E_rho(i) = sum(N_nonzer(i),'all')/L_c;
    
end

snr = real(non_demsi)/ real(noise_var);
% 
snr_avg = sum(snr)/length(snr);
% % %  This value should be in the range of [-0.5,1]
n = 0;  % Order of the Bessel function
% The average  Doppler shift is given by:
clear mean
if s==1
    avgD_shift= E_rho;
else
    avgD_shift = mean(E_rho);
end

bessel_value = abs ((1+(1/snr_avg)) * avgD_shift);
% abs(bessel_value);
% % % Example usage
n = 0;  % Order of the Bessel function
x = inverse_bessel(0,bessel_value );
f_max = x/(2*pi*Ts); % this to be at the N^t symbol
f_min  = x/(2*pi*N*Ts);
% 
% disp(['Inverse Bessel function:', num2str(x)]);

Delta_f=[f_max, f_min]; 
% 
% f_max1=  x/ (2*pi*Ts); % here we should follows the normalized figure as given in the reference book (page 109)
% 
V_max = c*f_max./f ;
V_min = c*f_min./f ;
