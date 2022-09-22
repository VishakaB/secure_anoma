EbNodB = 5;%snr 
Nerrs  = 0; Nblocks = 10;
k = 1;%message size
sigma = 1;%noise variance
n = 3;%number of users
for i = 1: Nblocks
    msg = randi([0,1],1,k);
    %encoding 
    cword = [msg msg msg];
    s = 1 - 2*cword;
    r = s+sigma*randn(1,n);
end
