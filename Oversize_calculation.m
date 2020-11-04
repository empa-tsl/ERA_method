% calcuation of oversize factor omega from Exiobase v.3.4
% (based on the approach from Cabernard et.al. 2019 and Dente et al. 2018)
% description can be found in Desing et al. 2020 https://doi.org/10.1017/sus.2020.26

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


% Read in coefficient matrix, final demand, environmental extensions,
% household emissions
datapath = ['SoSOS_IxI/Files/IOT_' num2str(year) '_ixi/'];
A = dlmread([datapath 'A.txt'],'\t',3,2); % A = coefficient matrix with 7987 rows x 7987 columns (7987 = 163 sectors x 49 regions)
Y = dlmread([datapath 'Y.txt'],'\t',3,2); % Y = final demand with 7987 rows x 343 columns (343 = 49 regions x 7 final demand categories)


% Calculate Leontief Inverse (L)
I = eye(size(A));
L = inv(I-A);

% Derive Total output of each region-sector combination to satisfy the global
% final demand (TotalOut)
X = L*Y;
TotalOut =sum(X,2);


% selecting target industries (t)
industries_selected = xlsread(filename,'dimensions','C3:C165');
[~,industries_names,~] = xlsread(filename,'dimensions','B3:B165');
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direct + supply chain output without double counting (wdc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Derive the total output without double counting (X_t_wdc) according to
% the method of Dente et al. 2018, Cabernard et.al. 2019
L_oo_dash = inv(I(index_o,index_o)-A(index_o,index_o));

X_t_wdc_C(:,1:343) = Y(index_t,:) + A(index_t,index_o) * L_oo_dash * Y(index_o,:); %link of X_t_wdc between target (rows) and final demand (column)


X_t = X(index_t,:);

X_t_wdc_GLO = zeros(n_target,1);
X_t_GLO = zeros(n_target,1);
for i=1:n_country
    for j=1:n_target
        X_t_wdc_GLO(j,1) = X_t_wdc_GLO(j,1) + sum(X_t_wdc_C((j+n_target*(i-1)),:),2); 
        X_t_GLO(j,1) = X_t_GLO(j,1) + sum(X_t((j+n_target*(i-1)),:),2); 
    end
end


Omega_r = sum(X_t,2) ./ sum(X_t_wdc_C,2);


Omega = X_t_GLO ./ X_t_wdc_GLO;




%%%%%
% writing resultss
%%%%%

xlswrite(filename,index_t_s','Omega','A5');
xlswrite(filename,industries_names(index_t_s),'Omega','B5');
xlswrite(filename,Omega,'Omega','C5');


toc