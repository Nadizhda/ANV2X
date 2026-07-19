% first adjust the code to run the python code in MATLAb
% 1- import sys
% 2- sys.excutable
% %%% add the path to the MATLAB codes 
%  3- pyenv('Version', '.../python')
% The system model for the adaptive subcarrier spacing and Tg using
% Jakes model for signal envelope simulation of Rayleigh fading
% clc

close all;

% clear all;

%% The number of pathes

N= 18;

M_path = 0.5*(N/2-1);% N/2 should be an odd integer.

Wn(M_path) = 0;

alpha(M_path) = 0;

riamp=[];

rqamp=[];

rialpha(1:2,:) = 0; 

%% Doppler Parameters

f_c= 5.9; 

% Realtive vecloities causing the doppler shift

v = [20,40,100,160,240]; % added to consider the second iteration 

% Doppler shift for sub 6GHz (Theoritical) 

fD = v*1e3*f_c/((3*10^8)*3600);


Wd = 2*pi*fD; % 
%% number of samples based on the symbol duration
%% FFT parameters 

N_SCS = [1024,512,256,128] ;

TB_s = [66.67e-6,33.33e-6,16.67e-6, 8.335e-6] ; % the length of the OFDM modulated symbol.

T_das = TB_s./N_SCS; %the length of the digital modulated symbol. (Not important)

Lc = N_SCS/4; % cyclic prefix : initial assumption

Gaud_int=[4.69e-6, 2.34e-6,1.17e-6,0.95e-6] ;

Ts = TB_s+ Gaud_int; % Total OFDM symbol interval 

N_g = Gaud_int./T_das; % The number of bins in the gurad time (cyclic prefix) 

Rat = floor(TB_s./Gaud_int); % multipliers of the OFDM symbol duration compared to the guard time


%% Channel parameters 
Eo = 1;
SNR = 25; % The target SNR

delt_fc = 0.00001;

scs = 15e3*2.^(0:3); % mmWave (2:3)// Sub 6 GHz (0:2)

noise_var = 10.^(-SNR/10);

N0 = 1;

N_iter=100;

A1=[];

time = 100; % ms

tau_max= 1e-6; % mean delay spread ; ignore the name  randi([1,600],1,1) * 1e-9;% The symbol duration. 

Bw =10e6; 


%% Transmission process


result1=[];

% SINR=zeros(length(fD),:);

Value=[];

BER=[];

Aenvl=zeros ([length(fD),length(N_SCS),1280]);

PDP=[];

tau=[];

An=[];

TF= 10;

MCS = 2;

Aenvl=[];



pd = makedist('normal','mu',tau_max);

tau = random(pd,[1,M_path]);

% PDP = exp(-tau/tau_max)* tau_max;

Mu0=0;

An=[];

R_set =  1; %get_RSet (Mu0, sim_time);

sinr_R=zeros(length(fD),R_set,1280); % sub 6GHz

Delta_f=[];

 Aenvl= zeros(length(fD),R_set,1280 ); % sub 6GHz


for i= 1: length(fD)
    % i
    % Rx_ve = [abs(Tx_ve(1) - v(i)),0,0];

    % for ii=1:length(N_SCS) % smulation time 

    % data_symb = [];

    % modulate the data bits (A)

    Data=(randi([0 1],1, 200*8)); % 7 RB

    % data_symb=mod_BPSK(Data); For the use of digital modulation 
    data_symb =  Data 
    
    L = length(data_symb);
    
    mode_symb=[];

    % initi_SC = N_SCS(1);
    % detrmine the number of resource set during the simulation time

   % initiate the design of the resource grid and the number of
   % resource sets withn the simulation time
   % The resource grid is instantaneous
   % fr is the resource set index
   % AvEnv= zeros()
   for fr = 0 : R_set
      % Non-Adaptive
       if fr== 0
           N_SC= N_SCS(1); % NFFT size
           scs_n=scs(1); % subcarrier spacing
           L_c= Lc(1); % number of frequency samples in the guard interval
           Ts_d=Ts(1); % The symbol duration without guard band
           Mu= Mu0 ; % consider any initila value of the SCS
           Tg = Gaud_int(1); % duration of guard interval s
           [mode_symb,~]=   OFDM_test(data_symb,N_SC,L_c); 
    
       else
          % adaptive 
          scs_n =An(i,fr,2);

          Tg = An(i,fr,6); % [Mod_data, ~] =   OFDM_test(data_symb,scs_n,Tg); 

          [ mode_symb,~,N_SC,scs_n,L_c,T_ds,Ts_n] = OFDM_adaptive(data_symb, scs_n,Tg);

           Mu= find(scs==scs_n)-1 ;
           Ts_d= TB_s(Mu+1)+ T_ds;
           % % Non adaptive 
           % N_SC= N_SCS(2); % NFFT size
           % scs_n=scs(2); % subcarrier spacing
           % L_c= Lc(2); % number of instance of the guard interval
           % Ts_d=Ts(2); % The symbol duration without guard band
           % Mu= Mu0 ; % consider any initila value of the SCS
           % Tg = T_das(Mu+1); % duration of guard interval s
           % [mode_symb,~]=   OFDM_test(data_symb,N_SC,L_c); 
       end

      
        output = Get_ResGrid (Mu,TF,Bw); % The resource grid is adjusted based on the slot size (1/2^mu)
         N_RB = output(1); % Numbe of resource blocks 
         N_TS = 10;  % Numbe of subframes in one frame

        % prepare the packet for transmission.

        time= (fr*N_TS+1:1/2^Mu:(fr+1)*N_TS);
        
        % Update the resource allocation
        
        for ts =1:length(time)

            valu = discretize(ts, [1,2])==1 ; 

            if valu ==1

                % Use the Contol Channel Channel 

                [CD(i,fr+1,:),~] =  Update_Contrdata (MCS,Mu,int16(L_c),int8(fD(i)),0); % Control Data 

                Rb_indx = randi([1,N_RB],1);

                Ts_indx = randi([3,time(end)-2],1);
        
                % Check the design of the resource set and update the N_Ts
             elseif ts>=13

                % send the feed back channel the SINR level
                % according to the CQI
                % and the numerology and the limit of the
                % Tg
                [CD(i,fr+1,:),~] =  Update_Contrdata (MCS,Mu,int16(L_c),int8(fD(i)),0);

            else

                % Check the new N_Ts

                Rb_indx = randi([1,N_RB],1);

                Ts_indx = randi([3,time(end)-2],1);

                SCS_out = scs_n;

            end

        end

         [sz,sz2] = size(mode_symb);

         t = 0 :Ts_d: Ts_d*sz*sz2;
 
        % Check if there is fed back about the SCS for the previose
        % transmission 
    
        riamp=[];
    
        rqamp=[];

        rialpha=[];

        ber=[];

        meen=0;
        %%%%%%%%%% Start transmission using Jakes spectrum under Montecarlo Simulations
        for k=N_iter:
           
            SINR2=[];
    
            for n= 1:M_path

               % The input here is modulated using OFDM signal 

                for tt = 1:length(t)-1
                  
                    Wn(n) = Wd(i)'.*cos(2*pi*n/N);% Doppler shift at each wave coming from multipath
    
                    alpha(n) = 2*pi*n/M_path; % phase shift due to the scatter object

                    % (tt*scs/2) Is the central frequency of the
                    % bandwidth part with respect to the subcarrier
                    % frequency

                    riamp(n,tt,:) = 0.2*cos(alpha(n))*cos(Wn(n)*(((t(tt) +(Ts_indx/(1/2^Mu)))*1e-3)-tau(n))).*real(mode_symb(tt) ) *...
                        cos(2*pi*(f_c+(tt*scs_n/2))*(((t(tt)+ +(Ts_indx/(1/(2^Mu))))*1e-3)-tau(n))); %The real value of the received signal (should be the modulated symbols)
    
                    rqamp(n,tt,:) = 0.2*sin(alpha(n))*cos(Wn(n)*(((t(tt) +(Ts_indx/(1/2^Mu)))*1e-3)-tau(n))).* imag(mode_symb(tt)) *...
                        sin (2*pi*(f_c+(tt*scs_n/2))*(((t(tt) +(Ts_indx/(1/(2^Mu))))*1e-3)-tau(n))); % The imaginary value of the recieved signal (should be the modulated symbols)
    
                    rialpha(1,tt)= sqrt(2)*cos(Wd(i)*t(tt));
                end
            end  

            Chann_signal=[];

            % the values of the ri and alpha

            ri = sum(riamp)+rialpha(1,:);
    
            rq = sum(rqamp );
            %%%%%%%%%% End transmission using Jakes spectrum under Montecarlo Simulations
            % The recieved signal at OFDM signal 
            Rcvd =[];
            
            Rcvd = ri+1j*rq; 
        
            % Add the more fading to the signal 
            % shadwing and blockage and other fading things to be updated here as an additional fading          
            h=1/sqrt(2)*(randn(sz,sz2)+1j*randn(sz,sz2));
        
            %One tap channel estimation 
            h_x=[];

            h_x(:,:)=  reshape(h,1,[]); %(:)*ones(1,sz2); 

            H_x= OFDM_test(h_x,N_SC,L_c);

            Noise=[];
                
            Noise = noise_var*sqrt(N0/2)*(randn(1,sz*sz2)+ 1j*randn(1,sz*sz2));
            
            N_noise = FFT(Noise, N_SC); %% FFT for the noise based on the selected numerology.

            h_x= reshape(h_x, 1,[]);

            Chann_signal = Rcvd (:,1:end,:).*h_x +  noise_var*sqrt(N0/2)*(randn(1,length(h_x))+ 1j*randn(1,length(h_x)));
              

            % In this place, I have to estimate the correct doppler shift
          

            % Doppler Effect function based on Reference [12]
             
             Delta_f(i,fr+1,:) = Doppler_Est(N_SC,noise_var,Ts_d,Chann_signal,tau_max,f_c,L_c);
            
            % Delta_f(isnan(Delta_f))=fD(i) ;
            % Direct estimation if Refrence [12] is not adopted. 
            % Delta_f(i,fr+1,:) = [abs(randn(1))* fD(i),0]; % Only based on random number upperbounded
    
            % check the optimization problem with respect to the current setings
            % after writing the data to a text file
            
            % Then apply the SINR function to evaluate the change in the SINR

            Chann_signal1= reshape(Chann_signal, [sz,sz2]);
            
            % OFDM demodulation 
            rceiv_symbol= OFDM_demo(Chann_signal1,L,N_SCS(Mu+1),L_c);

            % adaptive scheme use M+1=j
     
            % Digital Demodulation 
            
            Data_bits = BPSK_dem(rceiv_symbol);

            r = sqrt(real(Chann_signal).^2+imag(Chann_signal).^2 );

            meen = meen+ sum(r)/(length(t)-1);
            
        end
     %%%%%%%%% ANV2X %%%%%%%%%%
     I=[];
        delt_f(i)= fD(i)+ Delta_f(i,fr+1,1);
        epsilon(i,fr+1)= delt_f(i)/scs_n;
        b1= exp(1j*pi*epsilon(i,fr+1)*(N_SC-1)/N_SC) ;
        for lk=1:length(Chann_signal)
            Ik=0;
            for ik=1:length (Chann_signal)
                c1= sin(pi*(ik-lk+(delt_f(i)/scs_n)))/(N_SC*sin(pi*(ik-lk+(delt_f(i)/N_SC*scs_n))));
                if lk~=ik
                    Ik =Ik +abs(Chann_signal(ik)*Chann_signal(ik))*exp(1j*(ik-lk)*(N_SC-1)/N_SC)*c1;
                end
            end

           I(lk)=b1*Ik;

           a1= sin(pi*epsilon(i,fr+1))/(N_SC*sin(pi*epsilon(i,fr+1)/N_SC));
 
           R_rx(lk)= a1*b1* (Chann_signal(lk));
 
           sinr_R(i,fr+1,lk)= abs(R_rx(lk)/ (N_noise(1,lk)+ I(lk))); % Correct version Don't change it
           % bitsz_sinr(i,fr+1,lk)= log2(1+abs(sinr_R(i,fr+1,lk)));

        end
        % clear mean;
        Rate_1(i,fr+1,1:lk) = Rate_adaption(sinr_R(i,fr+1,1:lk)); % Rate mapping based on Table IV 
        Rate(i,fr+1) = mean(Rate_adaption(sinr_R(i,fr+1,1:lk)));

         N = length(Noise);
        
         % %  for SCIPY implementation 

         % Delta_f(i,fr+1,2) = Delta_f(i,fr+1,2)+ fD(i); % Calculation based on the reference paper
         for m=1:length(N_SCS)

           y0= (Gaud_int(m)/max(tau))* max(tau);

           rms_d = rms(Chann_signal);

           %  Writing the data into a text to be used in the python file for solving the SE maximization based on the resulting above data. 

           write_text(Bw,scs(m),2,5, N_SCS(m),min(scs),max(scs),3,5,4,32,min(N_SCS),max(N_SCS),Delta_f(i,fr+1,1) ,max(tau), N,y0,mode_symb, Noise,h_x,data_symb,rms_d);

           % SCS and T_g as an optimization variables for the frequency
           % based oprimization function 

            B(m,:) = double(pyrunfile("ANV2X_SCIPY.py","f_result"); %% Update the python file name 

            Result(i, fr+1,m,:) = [B(m,1), B(m,6)];
        end
        
        scs_ad= find (B(:,1)==max(B(:,1)));

        if length(scs_ad)~=1
            An(i, fr+1,:) = B(1,:) ;

        else

            An(i, fr+1,:) = B(scs_ad,:) ;
        end

        TSINR(i,fr+1,1:length(SINR)) = abs(SINR/N_iter);
    
        
        Value(i,fr+1,:) = An(i,fr+1,:)  ; 
        
        AvEnv(i,fr+1,1:length(r)) = 10*log10(r)-10*log10(meen/N_iter);

        % AR(i,j,fr+1) =  mean( Rate_adaption(r));% mean(Rate_adaption(abs(AvEnv(i,j,fr+1,1:length(r)))));
        AR1(i,fr+1) =  mean( Rate_adaption(TSINR(i,:))); % Rate adaption 

        % single numerology 
        AvEnv1=[];

        AvEnv1 = reshape(AvEnv(i,fr+1,:),1,[]);

        s= length(nonzeros( AvEnv1));

        Aenvl (i,fr+1,1:N) = AvEnv1(1:N)./ Noise ;
       
        % Rate_nonad(i,j,fr+1) = sum (scs(j) *log2(1+(TSINR(i,j,fr+1,1:s)/(N_iter*sz2)))); 

       % rau = Number of useful carriers*subcarrier spacing
       % ach_rate (i,ii)= rau * (log2(1+ real(abs(i,ii,:)))); to combare
       %  with the obtained 
       % result1(i,ii) = An(ii,:);

       MU_matrix(i,fr+1) = Mu;


        j=j+1;
       snr_li(i,fr+1,1:sz2) = SINR2/(sz2*N_iter);

  end

  
end
