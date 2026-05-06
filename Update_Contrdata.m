%  Evaluate the control data for between the Tx  and the Rx
function [Cont_data,L]  = Update_Contrdata (MCS,Mu,L_c,Fd,var)
%  [CD(i,j,:),~] =  Update_Contrdata (MCS,Mu,int32(L_c),int32(fD(i)),0)
    if var ==0
        Tg = 256;
        
        con_data = [int2bit(MCS,2); int2bit(Mu,2); int2bit(L_c,8)]';
        cont_mes = mod_BPSK (con_data);
        Cont_data  = [MCS,Mu,Tg,Fd]; % This should be updated only at the initaition of the transmission
    elseif var ==1
        % Rx side
        con_data = [int2bit(MCS,2); int2bit(Mu,2); int2bit(int(Tg),4),int2bit((Fd),log2(Fd))]';

        cont_mes = mod_BPSK (con_data);
        Cont_data  = [MCS,Mu,Tg,Fd];
    
    end
    
    L = length(cont_mes);

end