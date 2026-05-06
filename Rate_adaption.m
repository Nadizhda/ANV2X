%  Discrete Rate adaption 
function  [SE, CQI_u, BER] = Rate_adaption(SINR)
% SINR=[-7.5717,    5.0329,    1.9138,    4.2761,   -1.9799,    3.9392];
    CQI= 0:15;
    SNR= [-1.889, -0.817, 0.954, 2.984,4.889, 7.39, 8.898, 11.02, 13.32, 14.68, 16.62, 18.91,21.58, 24.88, 29.32];
    Rate = [0,0.1532,0.2344, 0.3770, 0.6016, 0.8770, 1.1758, 1.4766,1.9141, 2.4036, 2.7305, 3.3223, 3.9023, 4.5234,5.1152,5.5547];
    SE=[];
    CQI_u=[];
    BER=[];
    for j = 1:length(SINR)
        a=SINR(j);
        for i =1:length(SNR)
            if i>0 && i<length(SNR)
                if a < SNR(1)
                    se = Rate(1);
                    cqi=CQI(1);
                elseif SNR(i)<a && a <= SNR(i+1)
                    se = Rate(i+1);
                    cqi=CQI(i+1);
                elseif a> SNR(end)
                    se =Rate(end);
                    cqi=CQI(end);
                end
             end

        end


        SE(j) =se;
        CQI_u(j) = cqi;
    end
   % SE_av = sum(SE)/length(SE);
end