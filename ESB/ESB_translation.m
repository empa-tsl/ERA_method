% Preparing Earth system boundaries  for Monte Carlo simulations

clc
clear all

load('ATP.mat');
load('Land_appropriable.mat');

n_runs = 100000;
ESB = MCrand2(xlsread('ESB_translation.xlsx','PB','D24:G35'),n_runs);
ESB_distribution = xlsread('ESB_translation.xlsx','PB','G24:G35');
ESB(8,:,:) = Land_appropriable;
ESB(12,:,:) = E_total;

save('ESB.mat','ESB','ESB_distribution')

