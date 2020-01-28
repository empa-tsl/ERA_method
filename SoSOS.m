% calcuation of share of operating space (SoSOS) from Exiobase v.3.4
% A) direct emissions per resource segment 
% B) direct + supply chain emission allocated to the resource segment
% (according to Cabernard et.al. 2019)

clc
clear all
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0) preperation of data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameter
target_segment = 'materials';
year = 2011;



filename = ['SoSOS_IxI/SoSOS_' target_segment '.xlsx'];
filename2 = 'SoSOS_IxI/Files/Impact_assessment.xlsx';

addpath('SoSOS_IxI')
addpath('functions')
addpath('ESB')
load('ESB.mat');
[n_ESB,~,n_runs] = size(ESB);


% Read in coefficient matrix, final demand, environmental extensions,
% household emissions
datapath = ['SoSOS_IxI/Files/IOT_' num2str(year) '_ixi/'];
A = dlmread([datapath 'A.txt'],'\t',3,2); % A = coefficient matrix with 7987 rows x 7987 columns (7987 = 163 sectors x 49 regions)
Y = dlmread([datapath 'Y.txt'],'\t',3,2); % Y = final demand with 7987 rows x 343 columns (343 = 49 regions x 7 final demand categories)
datapath = ['SoSOS_IxI/Files/IOT_' num2str(year) '_ixi/satellite/'];
Ext = dlmread([datapath 'F.txt'],'\t',2,1); % Ext = environmental and social extensions (rows) for 7987 sector-region combinations (columns)
Ext_hh = dlmread([datapath 'F_hh.txt'],'\t',2,1); % Ext_hh = environmental and social extensions (rows) for households (columns = 49 regions x 7 final demand categories)

n_flows = size(Ext,1);

% Calculate Leontief Inverse (L)
I = eye(size(A));
L = inv(I-A);

% Derive Total output of each region-sector combination to satisfy the global
% final demand (TotalOut)
X = L*Y;
TotalOut =sum(X,2);


% selecting target industries (t)
industries_selected = xlsread(filename,'dimensions','C3:C165');
industries_selected(isnan(industries_selected))=0;
n_industries = size(industries_selected,1);
n_target = sum(industries_selected,1);
% building a vector with indix numbers of selected industries
k=1;
index_t_s = zeros(1, n_target);
for i=1:n_industries
    if industries_selected(i,1)==1
        index_t_s(1,k) = i;
        k=k+1;
    end
end
% Index of target regions (default: all regions on the globe)
index_t_r = [1:49];
% derive the index of non-target sector region combinations (O)
index_all_s = [1:163];
index_all_r = [1:49];

index_o_s = setdiff(index_all_s, index_t_s);
index_o_r = setdiff(index_all_r, index_t_r);
% target index for whole range (industry * region)
z = 0;
for n = index_t_r
    for i = index_t_s
        q_wdc = i + (163 * (n-1)) ;
        z= z+1;
        index_t(z)= q_wdc;
    end
end
n_t = length(index_t);
n_country = size(index_t_r,2);

index_all = [1:7987];
index_o = setdiff(index_all, index_t);
n_o = length(index_o);



% segment definition
S = xlsread(filename,'S','C5:P200');
S(isnan(S))=0;
n_segment = size(S,2);
S_t = zeros(n_target+2,n_segment);
S_t(1,end)=1;                   % final consumption
S_t(2,end-1)=1;                 % Rest of economy
S_t(3:end,:) = S(index_t_s,:);  % target industries



% impact assessment on PB
% impact characterization
Q_raw = xlsread(filename2,'Q','D5:AY1108');
Q_raw(isnan(Q_raw))=0;
% normalized uncertainty based on inventory flows (mean = 1)
addpath('SoSOS_IxI/functions/')
LCI_uncertainty = MCrand2(xlsread(filename2,'LCI_uncertainty','E5:H1108'),n_runs);

Q = zeros(n_flows,n_ESB,n_runs);
for i=1:n_ESB
    Q(:,i,:) = MCrand2(Q_raw(:,(1+4*(i-1)):(4+4*(i-1))),n_runs);
end
for i=1:n_runs
    Q(:,:,i) = diag(LCI_uncertainty(:,:,i)) * Q(:,:,i);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A) direct flows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flows
% final demand (household)
f_hh = sum(Ext_hh,2);
% total global flows
f_total_global = sum(Ext,2) + f_hh;

% flows from the target sector for all regions
f_T_direct_r = Ext(:,index_t);
% flows for target sector, global
f_T_direct = zeros(size(Ext,1),n_target);
for i=1:n_country
    for j=1:n_target
        f_T_direct(:,j) = f_T_direct(:,j) + f_T_direct_r(:,(j+n_target*(i-1))); 
    end
end
% flows from non-targets (t)
f_O_direct = sum(Ext(:,index_o),2);

f_direct = zeros(n_flows,n_target+2);
f_direct(:,1) = f_hh;
f_direct(:,2) = f_O_direct;
f_direct(:,3:end) = f_T_direct;

% impacts
q_direct = zeros(n_ESB,(n_target+2),n_runs);
q_S_direct = zeros(n_ESB,n_segment,n_runs);
for i=1:n_runs
    q_direct(:,:,i) = Q(:,:,i)' * f_direct;
    q_S_direct(:,:,i) = q_direct(:,:,i) * S_t;
end
q_direct_total = sum(q_direct,2);

% SoSOS
SoSOS_S_direct = zeros(n_ESB,(n_segment),n_runs);
for i=1:(n_segment)
   SoSOS_S_direct(:,i,:)=q_S_direct(:,i,:) ./ q_direct_total;
end

SoSOS_S_q_direct(:,:) = quantile(SoSOS_S_direct,0.5,3);   % average value for display
SoSOS_S_q_total_direct = sum(SoSOS_S_q_direct,2);
for i=1:n_ESB
    SoSOS_S_q_direct(i,:) = SoSOS_S_q_direct(i,:)/SoSOS_S_q_total_direct(i);     %scaling to 1 for balance
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B) direct + supply chain flows without double counting (wdc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Derive the total output without double counting (X_t_wdc) according to
% the method of Dente et al. 2018, Cabernard et.al. 2019
L_oo_dash = inv(I(index_o,index_o)-A(index_o,index_o));

X_t_wdc_C(:,1:343) = Y(index_t,:) + A(index_t,index_o) * L_oo_dash * Y(index_o,:); %link of X_t_wdc between target (rows) and final demand (column)

% Calculation of impact coefficients of the economy 
D = Ext ./ repmat(TotalOut',length(sum(Ext,2)),1);
D(isinf(D))=0 ; 
D(isnan(D))=0 ;  
D = D .* (D>=0); 
% calculating all flows for wdc
f_T_wdc_r = D * L(:,index_t) * diag(sum(X_t_wdc_C,2));
TotalOut_t = TotalOut(index_t);
Ext_t = Ext(:,index_t);
for i=1:n_t
    if TotalOut_t(i) == 0
        f_T_wdc_r(:,i) = f_T_wdc_r(:,i) + Ext_t(:,i);
    end
end
f_global_wdc = f_total_global;
f_O_wdc = f_global_wdc - sum(f_T_wdc_r,2) - f_hh;
% deviation from direct emission through conversion into coefficients
%f_total_check = (f_total_global - f_global_wdc) ./ f_total_global;


f_T_wdc = zeros(size(Ext,1),n_target);
for i=1:n_country
    for j=1:n_target
        f_T_wdc(:,j) = f_T_wdc(:,j) + f_T_wdc_r(:,(j+n_target*(i-1))); 
    end
end

f_wdc = zeros(n_flows,n_target+2);
f_wdc(:,1) = f_hh;
f_wdc(:,2) = f_O_wdc;
f_wdc(:,3:end) = f_T_wdc;


% impacts
q_wdc = zeros(n_ESB,(n_target+2),n_runs);
q_S_wdc = zeros(n_ESB,n_segment,n_runs);
for i=1:n_runs
    q_wdc(:,:,i) = Q(:,:,i)' * f_wdc;
    q_S_wdc(:,:,i) = q_wdc(:,:,i) * S_t;
end
q_wdc_total = sum(q_wdc,2);


% share of total impact
SoSOS_wdc = zeros(n_ESB,n_target+2,n_runs);
for i=1:n_target+2
   SoSOS_wdc(:,i,:)=q_wdc(:,i,:) ./ q_wdc_total;
end
SoSOS_S_wdc = zeros(n_ESB,n_segment,n_runs);
for i=1:n_segment
   SoSOS_S_wdc(:,i,:)=q_S_wdc(:,i,:) ./ q_wdc_total;
end

SoSOS_S_q_wdc(:,:) = quantile(SoSOS_S_wdc,0.5,3);   % average value for display
SoSOS_S_q_total = sum(SoSOS_S_q_wdc,2);
for i=1:n_ESB
    SoSOS_S_q_wdc(i,:) = SoSOS_S_q_wdc(i,:)/SoSOS_S_q_total(i);     %scaling to 1 for balance
end

% re-calculate total impacts per segment from total
q_S_wdc_re = SoSOS_S_wdc .* q_direct_total;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save and display results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save data
datapath = ['SoSOS_IxI/results/' target_segment '/'];
    save([datapath 'f_total_global_' target_segment '_' num2str(year)],'f_total_global');
    save([datapath 'f_wdc_' target_segment '_' num2str(year)],'f_wdc');
    save([datapath 'f_direct_' target_segment '_' num2str(year)],'f_direct');
    save([datapath 'S_t'],'S_t');
    save([datapath 'index_t_s'],'index_t_s');
    save([datapath 'SoSOS_S_wdc_' target_segment '_' num2str(year) '.mat'],'SoSOS_S_wdc');
    save([datapath 'SoSOS_wdc_' target_segment '_' num2str(year) '.mat'],'SoSOS_wdc');
    save([datapath 'SoSOS_S_direct_' target_segment '_' num2str(year) '.mat'],'SoSOS_S_direct');
    save([datapath 'q_S_direct' target_segment '_' num2str(year) '.mat'],'q_S_direct');
    save([datapath 'q_S_wdc_re_' target_segment '_' num2str(year) '.mat'],'q_S_wdc_re');

% save to xls
i=1;
if i==1
    xlswrite(filename,f_total_global,'f','E5')
    xlswrite(filename,f_wdc * S_t,'f_S','F5')
    xlswrite(filename,f_direct * S_t,'f_S_direct','F5')
    xlswrite(filename,index_t_s','dimensions','O5')

    xlswrite(filename,quantile(q_direct_total,0.01,3),'f_S_direct','E1111');
    xlswrite(filename,quantile(q_direct_total,0.5,3),'f_S_direct','E1124');
    xlswrite(filename,quantile(q_direct_total,0.99,3),'f_S_direct','E1137');

    xlswrite(filename,quantile(q_S_wdc_re,0.01,3),'f_S','F1111');
    xlswrite(filename,quantile(q_S_wdc_re,0.5,3),'f_S','F1124');
    xlswrite(filename,quantile(q_S_wdc_re,0.99,3),'f_S','F1137');

    xlswrite(filename,quantile(q_S_direct,0.01,3),'f_S_direct','F1111');
    xlswrite(filename,quantile(q_S_direct,0.5,3),'f_S_direct','F1124');
    xlswrite(filename,quantile(q_S_direct,0.99,3),'f_S_direct','F1137');
    
    xlswrite(filename,quantile(SoSOS_S_direct,0.01,3),'SoSOS','F6');
    xlswrite(filename,quantile(SoSOS_S_direct,0.5,3),'SoSOS','F19');
    xlswrite(filename,quantile(SoSOS_S_direct,0.99,3),'SoSOS','F32');
    
    xlswrite(filename,quantile(SoSOS_S_wdc,0.01,3),'SoSOS','F48');
    xlswrite(filename,quantile(SoSOS_S_wdc,0.5,3),'SoSOS','F61');
    xlswrite(filename,quantile(SoSOS_S_wdc,0.99,3),'SoSOS','F74');

end

toc