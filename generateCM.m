% load your structrual connectivity matrix here
% the dimension of the example matrix here is 420*86*86
sc = load('example_sc.mat');

%% Below are five different communication models(CMs) used in the paper, you can find the detail of the description of the CMs in the paper

% Search Information
SIsc = []; % array to save the transformed SC 
% loop all subjects and compute SIsc for each subject
for k = 1:420 % you can change 420 to the number of subjects in your dataset 
    weight = reshape(sc.data(k,:,:),86,86);
    length = - log10(weight / max(max(weight))); % map the weight matrix to length matrix
    SI = search_information(weight, length);
    SI = (SI + SI.') / 2; % make the matrix symmetric
    SIsc = [SIsc; {SI}];  
end

% Shortest Path Efficiency
SPEsc = [];
for k = 1:420
    weight = reshape(sc.data(k,:,:),86,86);
    length = - log10(weight / max(max(weight)));
    [SPL,hops,Pmat] = distance_wei_floyd(length);
    SPE = 1 ./ SPL;
    SPE = (SPE + SPE.') / 2;
    SPEsc = [SPEsc; {SPE}];
end

% Navigation Efficiency
% note that navigation efficiency requires distance matrix which is the
% eculidean distance between regions, you need to calculate it for each subject by yourself
% and store it in Distance variable
NEsc = [];
Distance = load('example_distance.mat');
Distance = Distance.Distance;
for k = 1:420
    weight = reshape(sc.data(k,:,:),86,86);
    length = - log10(weight / max(max(weight)));
    D = Distance{k}; % this is your distance matrix 
    [sr, PL_bin, PL_wei, PL_dis, paths] = navigation_wu(length, D);
    PL_eff = 1 ./ PL_dis;
    PL_eff = (PL_eff + PL_eff.') / 2;
    NEsc = [NEsc; {PL_eff}];
end

% Diffusion Efficiency
DEsc = [];
for k = 1:420
   [GEdiff,Ediff] = diffusion_efficiency(reshape(sc.data(k,:,:),86,86));
   Ediff = (Ediff + Ediff.') / 2;
   DEsc = [DEsc; {Ediff}];
end

% Communicability
CMYsc = [];
for k = 1:420
    weight = reshape(sc.data(k,:,:),86,86);
    degree = diag(sum(weight,2));
    CMY = exp(degree^(-1/2) * weight * degree^(-1/2));
    CMY = (CMY + CMY.')./2;
    CMYsc = [CMYsc;{CMY}];
end
