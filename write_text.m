% write the output to the texxt file that will be used as an input to the
% python code
% num_SC is the current N_size IFFT
% SCS is teh current subcarier spacing
% Mod_Sy : should be used to validate what is the bext SCS with N_IFFT size
% that gives the maximum spectral effeiciency with respect to the current
% received signal. 
function  write_text(Bw,SCS0,alpha0, beta0, N_IFFT0,SCS_min,SCS_max,alphamin, alpamax, betamin, betamax,N_min, N_max,f_d,tau_max,N,y0,tra_symbol,N_noise,h_x,Mod_Sy,rms_d)

  f=fopen('text.txt','w');
  fprintf(f,'%E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E %E',Bw,SCS0,alpha0, beta0, N_IFFT0,SCS_min,SCS_max,alphamin, alpamax, betamin, betamax,N_min, N_max,f_d,tau_max,N,y0,rms_d);
  fprintf(f,' %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f %6.5f ',real(tra_symbol),...
  imag(tra_symbol),real(N_noise),imag(N_noise),real(h_x),imag(h_x),real(Mod_Sy),imag(Mod_Sy),rms_d);
  fgetl(f);
  fclose(f);
end

