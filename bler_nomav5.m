%start_date:22.09.2022
%last update: 22.09.2022

%goal: average latency vs transmit SNR
clc; 
clear all;
close all;

%% initial parameters
transmit_snrdb_vec = linspace(10,20,20);
noisepower = 0.1;
u_i = 5;
v_i = 2;
mu = 0.7;
n = 1000;%number of bits
m = 100; %blocklength
v = rand(1,n);
K = 3;%nbusers
sumc = zeros(K,1); 
delta = rand(n,n);
w_j = 0.7;
P_j = 1;
alpha_j = 0.9;
T_sym = 0.007;
final_analy_latency = zeros(length(transmit_snrdb_vec),1);
final_simula_latency = zeros(length(transmit_snrdb_vec),1);

%% analytical delay 
for idx =1 : length(transmit_snrdb_vec)
    
    ENodB = transmit_snrdb_vec(idx);%snr
    
    for k = 1:K  %user index
      P_k(k) = 10^(ENodB/10)*noisepower;%tx power
      for s = 1:m %symbol index %for each symbol
        for zeta = s-1:1:s+1 %neighbor symbol index
             if(zeta>1)
                sumc(k) = (1-delta(k,zeta)+1) + sumc(k);
             else
                sumc(k) = (1- delta(k,zeta+1)); 
             end
        end
   
        if(zeta+1==s)%check here
          deltaval = 1;
        else
          deltaval = 0;
        end

        clear delta_mat;
        clear reverse_delta_mat;
        %time offsets between users 
        %delta_mat: rows -> user index, columns-> symbol index %time offset with
        delta_mat = zeros(K,m);%user max: K, symbol max: m

        %the total interference from both s-1th, s+1th symbols to the desired sth 
        %delta_mat(1,:) = zeros(1,m);%check here
        time_offset = rand(K-1,m);
        delta_mat(1:K-1,:) = time_offset.*ones(K-1,m);%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

        reverse_delta_mat = zeros(K,m);
        rev_time_offset = rand(K-1,m);
        reverse_delta_mat(2:K,:) = rev_time_offset.*ones(K-1,m);%A1, A2

        del_offset = delta_mat(k,s);
        num1 = P_k(k);
        den1 = ((1-deltaval)*del_offset*P_k(k)/2+w_j*P_j*alpha_j^2 ...
        +noisepower^2);
        bler_ther_p1 = (2*pi*(2^(2*n/m)-1)/m)^(-1/2);
        bler_ther_p2 = (u_i - v_i);
        bler_ther_p3 = (num1/den1)^2*(exp(-v_i/(num1/den1))/v_i^2 ...
            - exp(-u_i/(num1/den1))/u_i^2);%check this

        bler_ther(k) = bler_ther_p1*(bler_ther_p2+bler_ther_p3);

        analy_latency(k,:) = mu*(m*T_sym+v)/(1 - bler_ther(k)); 

        avg_analy_latency(idx,k) = mean(analy_latency(k,:)); %bler of each user     
      end
    end
    final_analy_latency(idx) = mean(avg_analy_latency(idx,k),2);
end
final_analy_latency
%% simulation delay 
%for reps = 1:10%check here 
%nbrandom simulations 
for idx =1 : length(transmit_snrdb_vec)
%% 
for k = 1:K  
%generate a sequence of 1s and 0s with 0.8% of 1s randomly distributed 
%1s = error block, 0s = no error block (correct)
numberOfElements = m; % Whatever.
percentageOfOnes = abs(mean(bler_ther(k)))*10; % Whatever.
numberOfOnes = round(numberOfElements * percentageOfOnes / 100);
% Make initial signal with proper number of 0's and 1's.
blckstream = [ones(1, numberOfOnes), zeros(1, numberOfElements - numberOfOnes)];
% Scramble them up with randperm
blckstream = blckstream(randperm(length(blckstream)));
% Count them just to prove it
numOnes = sum(blckstream);

%% latency 
% number of retransmissions required
% until the bler is < 10
% goal: retransmissions required
nbretransmissions = 0;
success = 0;
blockidx = 2;

%count the place where the first 1 appears in the signal stream
while (success == 0)
   if blckstream(blockidx-1)== 0
       success = 1;
   elseif blckstream(blockidx-1)== 1 & blckstream(blockidx)== 0
       nbretransmissions = nbretransmissions + 1;
       success = 1;
   else
       nbretransmissions = nbretransmissions + 1; 
       blockidx = blockidx + 1; 
       success = 0;
       if (blockidx>=length(blckstream))
            %initialize block stream 
            blockidx = 2;%initialize block stream index
            %generate a sequence of 1s and 0s with 0.8% of 1s randomly distributed 
            %1s = error block, 0s = no error block (correct)
            numberOfElements = m; % Whatever.
            percentageOfOnes = abs(bler_ther(k))*10; % Whatever.
            numberOfOnes = round(numberOfElements * percentageOfOnes / 100);
            % Make initial signal with proper number of 0's and 1's.
            blckstream = [ones(1, numberOfOnes), zeros(1, numberOfElements - numberOfOnes)];
            % Scramble them up with randperm
            blckstream = blckstream(randperm(length(blckstream)));
            % Count them just to prove it
            %numOnes = sum(blckstream);
        end       
   end
end 

delay_transmission = T_sym*m;
sim_avg_latency(idx,k) = nbretransmissions*delay_transmission;
end%end k=1:K
final_simula_latency(idx)= mean(sim_avg_latency(idx),2);
end

%% plots
grid on;
%plot(transmit_snrdb_vec,(final_simula_latency),'--r')
hold on 
plot(transmit_snrdb_vec,(final_analy_latency),'-b')
xlabel('SNR') 
ylabel('avg latency')