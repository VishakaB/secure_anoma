%example 
m  = 10;
u_i = 1;%lower sinr level
v_i = 0.1; %higher sinr level

Nk = 100;%number of bits transmitted by each user
eta_k = 0.2;%interference created by asynchronous neighbor symbols 
sigma = 0.1;%noise power
p_k = 2; %D2D transmission power

term1 = p_k^2/(eta_k+sigma^2);
te1 = -v_i/(p_k^2/(eta_k+sigma^2));
term2 = exp(te1)/v_i^2;
term3 = -u_i/(p_k^2/(eta_k+sigma^2));

bler_th = [2*pi*(2^(2*Nk/m)-1)/m]^(-1/2)*(u_i - v_i + term1.*(term2-term3))

