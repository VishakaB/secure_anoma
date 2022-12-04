%start_date:22.09.2022
%last update: 22.09.2022

%goal: average latency vs transmit SNR
clc; 
clear all;
close all;

%% initial parameters
transmit_snrdb_vec = linspace(10,50,20);
noisepower = 0.1;

u_i = 
v_i = 