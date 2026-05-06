% Get Resource grid 
% initial setting for the test onle
%  Mu,TF is the numerology value and the frame duraion respectively
% Mu, TF shold be announced by the cluster head. 
% This function should be implemented in the clutster memebers 
% To be implemented based on fuzzy logic system 
% aver_Rvel is teh average relative velocity 
function output = Get_ResGrid (Mu,TF,Bw)
%  C_ID will be used to get the resource grid configuration based on the
%  average veloity and the density of vehilces 
%   Mu = 0; 
%   TF = 10; % ms
  SCS = 15e3 * 2^Mu;
%   Bw = SCS*12;
  T_Sym = 1/SCS;
  
  Ts = 1/ (2^(Mu));
%   N_Ts = TF*1e3/ Ts; 
   N_Ts = TF/ Ts; 
  
  if Mu == 0
      N_RB = 50;
  elseif Mu == 1
      N_RB = 24;
  elseif Mu==2
      N_RB = 10;
  else
      N_RB = ceil(Bw/(SCS*12));
  end
  
  output  = [N_RB, N_Ts] ;
      
end
