clear all;
close all;
clc;

%goal: block error rate based latency simulation comparision with
%analytical
%% initial parameters
bler = 0.8;

%% 
%generate a sequence of 1s and 0s with 0.8% of 1s randomly distributed 
%1s = error block, 0s = no error block (correct)
numberOfElements = 10; % Whatever.
percentageOfOnes = 50; % Whatever.
numberOfOnes = round(numberOfElements * percentageOfOnes / 100);
% Make initial signal with proper number of 0's and 1's.
blckstream = [ones(1, numberOfOnes), zeros(1, numberOfElements - numberOfOnes)];
% Scramble them up with randperm
blckstream = blckstream(randperm(length(blckstream)));
% Count them just to prove it
numOnes = sum(blckstream);
blckstream

%% latency 
% number of retransmissions required
%until the bler is < 10
%goal: retransmissions required
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
            %renew transmission blockstream
            numberOfElements = 10; % Whatever.
            percentageOfOnes = 50; % Whatever.
            numberOfOnes = round(numberOfElements * percentageOfOnes / 100);
            %Make initial signal with proper number of 0's and 1's.
            blckstream = [ones(1, numberOfOnes), zeros(1, numberOfElements - numberOfOnes)];
            % Scramble them up with randperm
            blckstream = blckstream(randperm(length(blckstream)));
        end       
   end
end 

n = 100;%number of bits
m = 10; %blocklength
T_sym = 0.007;
delay_transmission = T_sym*m;
fprintf('nbretransmission %f\n',nbretransmissions);
sim_avg_latency = nbretransmissions*delay_transmission;
fprintf('sim delay %f\n',sim_avg_latency);

%% analytical delay 
v = rand(1,n);
u_i = 0.5;
v_i = 0.1;
P_k = 1;
zeta = 3;
s = 1;
delta_offset = 0.5;
mu = 10;
%n = 10;%total nbits

bler_ther_p1 = (2*pi*(2^(2*n/m)-1)/m)^(-1/2);
delta = rand(n,n);
w_j = 0.1;
P_j = 0.01;
alpha_j = 0.9;
sigma = 0.1;
num1 = P_k^2;
K = 3;%nbusers
sumc = zeros(K,1); 

for k = 1:K
for zeta = s-1:1:s+1
 if(zeta>1)
    sumc(k) = (1-delta(k,zeta)+1)+sumc(k);
 else
    sumc(k) = (1- delta(k,zeta+1)); 
 end
end
end

if(zeta<s)
  deltaval = 0;
else
  deltaval = 1;
end

clear delta_mat;
clear reverse_delta_mat;
%time offsets between users 
%delta_mat: rows -> user index, columns-> symbol index %time offset with
delta_mat = zeros(K,K);

%the total interference from both s-1th, s+1th symbols to the desired sth 
delta_mat(1,:) = zeros(K,1);
time_offset = rand(1);
delta_mat(2:K,:) = time_offset.*ones(K-1,K);%B1,... Bn, C1....,Cn, ....... %Z1,....Zn

reverse_delta_mat(K,:) = zeros(K,1);
reverse_delta_mat(1:K-1,:) = time_offset*ones(K-1,K);%A1, A2

del_offset = delta_mat(k,s);

den1 = ((1-deltaval+1)*del_offset*P_k/2+w_j*P_j*alpha_j^2 ...
+sigma^2)^2;

bler_ther_p2 = (u_i - v_i)+num1/den1;
bler_ther_p3 = (exp(-u_i/den1)/u_i^2 - exp(-v_i/den1)/v_i^2);

bler_ther = bler_ther_p1*bler_ther_p2*bler_ther_p3;

analy_avg_latency  = mu*(m*T_sym+v)/(1 - bler_ther);
fprintf('analytical delay %f\n',mean(analy_avg_latency));
