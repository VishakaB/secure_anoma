%start date: Sep 3 2022
%last update: 07/09/2022

%Title: Optimization for URLLC performance in A-NOMA assisted D2D communication 
%goal1: average BLER metric analysis analytical
%goal2: comparison with simulation results

%goal2: latency metric analysis
%goal3: secrecy capacity metric analysis
%goal4: optimization of ursllc metrics

%% comparison analytical vs simulation

%analytical
close all;
clear all;
clc;

d1 = 100;
eta = 4;
p1 = 100;
users=2;            % Number of Users

%------------------ Generating data for User1 -------------------------------
N=10;                               % Number of Bits for  data_user1
data_user1= rand(1,N)>0.5;          % Generation of data for user1
data_user1bpsk = (2*data_user1-1);  % BPSK modulation 0 -> -1; 1 -> 0 

%creating rayleigh fading channel 
gain1=sqrt(d1^-eta)*sqrt(p1/2)*[randn(1,N).^2 + j*randn(1,N).^2];% Gain for Tap1

%passing data through channel 
data_channel1=data_user1bpsk.*gain1;

%simulation

%% optimization 
%maximize number of decoded symbols 
%under BLER and latency thresholds
