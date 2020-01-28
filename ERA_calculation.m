% calculting ERA for materials
tic
clc
clear all

%%%%%%%%%%%%%%%%%%
% Inupt parameter:
%%%%%%%%%%%%%%%%%%
sector_name = 'metals';      % change sector name here
segment_number = 5;          % number of segment [1,13]
P_v = 0.01;                  % probability of violation
%%%%%%%%%%%%%%%%%%

addpath('functions/');

%%%%%%%%%%%%%%%%%%
% reading input files

filename = ['UI_SoP/ERA_' sector_name '.xlsx'];      
figurename = ['display/ERA_' sector_name];

P_cut = P_v/2;              % cut off probability for quantile calculation

load('SoSOS_IxI/results/materials/SoSOS_S_wdc_materials_2011.mat')
load('ESB/ESB.mat')
load('SoSOS_IxI/results/materials/index_t_s.mat')

SoP_cur = xlsread(filename,'SoP','H5:H104');

UI_raw = xlsread(filename,'UI','G5:BB104'); % Unit impacts, imput data
UI_raw(isnan(UI_raw))=0;

n_ESB = size(UI_raw,2)/4;
n_resources = size(UI_raw,1);
n_runs = size(ESB,3);

LCI_uncertainty = MCrand2(xlsread('files/LCI_uncertainty.xlsx','LCI_uncertainty','E5:H30'),n_runs);

%%%%%%%%%%%%%%%%%%
% building calculation matrices
% Unit impacts (UI)
UI=zeros(n_resources,n_ESB,n_runs);
for i=1:n_ESB
    UI(:,i,:) = MCrand2(UI_raw(:,(1+4*(i-1)):(4+4*(i-1))),n_runs);
end
for i=1:n_runs
    UI(:,:,i) = UI(:,:,i) * diag(LCI_uncertainty(:,:,i));
end

% Share of safe operating space for sector (SoSOS)
SoSOS_k = SoSOS_S_wdc(:,segment_number,:);
for i=1:n_runs
    for j=1:n_ESB
        if SoSOS_k(j,:,i)==0
            SoSOS_k(j,:,i)=0.00001;    % to avoid restrictions for emissions not accounted in Exiobase, but Ecoinvent
        end
    end
end

% segment boundary (SB)
SB = quantile(SoSOS_k,0.5,3) .* ESB;     % excluding uncertainty from SoSOS
%SB = SoSOS_k .* ESB;                    % including uncertainty from SoSOS

SB_limit = quantile(SB,P_cut,3);
UI_limit = quantile(UI,1-P_cut,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERA calculation with current production share (mode in SoP)
UI_k_cur=zeros(1,n_ESB,n_runs);
for i=1:n_runs
        UI_k_cur(:,:,i) = SoP_cur' * UI(:,:,i);
end
% calculate ERA for current SoP
ERA_k_cur =  SB_limit ./ (quantile(UI_k_cur,1-P_cut,3))';
[ERA_limit_cur, limiting_boundary_cur] = min(ERA_k_cur);

[ERA_total_cur, P_v_check_cur, count_cur] = ERA_adjustment_to_P_v(ERA_limit_cur, UI_k_cur, SB, P_v, n_ESB, n_runs);
[~,limiting_boundary_cur_rev] = max(P_v_check_cur);

ERA_cur = SoP_cur * ERA_total_cur;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save data
xlswrite(filename,ERA_cur,'result','G5')

toc   