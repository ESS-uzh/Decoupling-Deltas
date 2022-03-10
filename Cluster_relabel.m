%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Re-labelling of ecosystem service clusters                          %%
%  Code to re-label the output of an ensemble of cluster algorithms    %%
%  Reader et al., Communications Earth & Environment, 2022             %%
%  Re-labelling maximises the consistency between cluster algorithms   %%
%  To limit the number of permutations, the script looks for           %%
%  maximal consistency when fixing the cluster labels for one method   %%
%  Required input: The file 'DEScluster.xlsx'                          %%
%  Last update: March 03, 2022                                         %%
%  Questions regarding the code: martin.reader@geo.uzh.ch or           %%
%  maarten.eppinga@geo.uzh.ch                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear the current workspace
clear 

%% Load Data and define parameters describing the data structure
[Data,~,~] = xlsread('DEScluster.xlsx');              % Load the dataset
NoClus     = length(unique(Data));             % Number of clusters to be considered for relabelling
NoVar      = length(Data);                     % Number of variables included in the cluster analysis
NoMeth     = numel(Data)/NoVar;                % Number of methods used to cluster the data
OrgVal     = 1:NoClus;                         % Set of original labels assigned to clusters

%% Create matrices to enable consideration of the possible re-labelling options
AllCombs          = perms(1:NoClus);           % Create the possible cluster orders
AllOrders         = perms(1:NoMeth);           % Create the possible method orders
[NoCombs,~]       = size(AllCombs);            % Count number of cluster orders
[NoOrders,~]      = size(AllOrders);           % Count the number of method orders
ComboMatrix       = nan(NoClus,NoMeth);        % Create matrix format for labelling
ComboMatricesBest = nan(NoClus,NoMeth,NoMeth); % Create matrix to store label output
AllData           = nan(NoVar,NoMeth,NoMeth);  % Create matrix to store cluster designations

%% Start re-labelling procedure
for h = NoMeth:-1:1                            % Select one clustering method to fix
    BestFit       = 0;                         % First attempt will be chosen as best fit, and subsequently replaced
    DataMh        = Data;                      % Create temporary data file
    ComboMatrixh  = ComboMatrix;               % Create temporary matrix of label combinations
    X = find(AllOrders(:,1)==h);               % Select order combinations starting with specific method
    for i=1:length(X)                          % For all combinations, 
        ComboMatrixh(:,h) = 1:NoClus;          % Set the first method's labels fixed (using the raw data labels)
        DataTemp = DataMh;                     % Store data temporarily, to be overwritten by better fits
        DataMhC  = DataMh;                     % Store data temporarily, to be overwritten by better fits
        ComboMatrixhC = ComboMatrix;           % Store data temporarily, to be overwritten by better fits
        for j = 2:NoMeth                       % For all other methods (except the first one fixed above),
            BestFitC = 0;                      % First attempt will be chosen as best fit, and subsequently replaced
            for r = 1:NoCombs                  % For all possible cluster orders, perform the re-labelling:
                DataTemp(:,AllOrders(X(i),j)) = changem(DataMhC(:,AllOrders(X(i),j)),AllCombs(r,:),OrgVal);
                Freqs = nan(NoVar,NoClus);     % Create dataframe to evaluate fit
                for z = 1:NoClus               % Measure consistency of labelling for all clusters
                    Freqs(:,z) =sum(DataTemp(:, AllOrders(X(i),1:j))==z,2);
                end
                
                TotalFitC = sum(max(Freqs,[],2)); % Calculate total level of consistency among clusters so far
                if  TotalFitC > BestFitC          % If consistency is higher than previous best fit, keep it
                    BestFitC = TotalFitC;         %
                    BestSubC = AllCombs(r,:);     %          
                end
            end
            ComboMatrixhC(:,AllOrders(X(i),j))  = BestSubC; % Store most consistent labelling for clusters considered so far 
            DataMhC(:,AllOrders(X(i),j))        = changem(DataMh(:,AllOrders(X(i),j)),BestSubC,OrgVal);
        end                                       %
        
        Freqs = nan(NoVar,NoClus);                % Calculate the total degree of consistency among all clusters
        for z = 1:NoClus                          %
            Freqs(:,z) = sum(DataMhC==z,2);       %
        end                                       %
        
        TotalFit = sum(max(Freqs,[],2));          % If this particular cluster order provides higher consistency, keep it 
        if TotalFit > BestFit                     %
                BestFit      = TotalFit;          %
                DataMh       = DataMhC;           %
                ComboMatrixh = ComboMatrixhC;     %  
        end
    end
    
    ComboMatricesBest(:,:,h) = ComboMatrixh;      % Store the most consistent ordering, when method h is used first
    AllData(:,:,h)           = DataMh;            % Store the data when labelled accordingly
end

for qq = 1:NoMeth                                          % Out of all combinations:                            
    Freqs = nan(length(sum(AllData(:,:,qq)==z,2)),NoClus); % Determine the most consistent fit
    for z = 1:NoClus                                       %
        Freqs(:,z) = sum(AllData(:,:,qq)==z,2);            % 
    end                                                    %
    TotalFit(qq) = sum(max(Freqs,[],2));                   %
end                                                        %
Opt = find(TotalFit == max(TotalFit));                     %

CodingSystem = ComboMatricesBest(:,:,Opt);                 % Creates a matrix with the re-labelling system
Transformed  = AllData(:,:,min(Opt));                      % Creates a matrix with the re-labelled clusters

%% Create a visualization of the outcome:
figure;subplot(1,2,1); imagesc(Data,[1 NoClus]);           % Shows the original labelling system
title('Original cluster labels')                           % of the input data
subplot(1,2,2); imagesc(Transformed,[1 NoClus])            % Shows the re-labelled system
title('Results of cluster re-labelling')                   % with maximum consistency between methods (i.e. rows)
%% Save the re-labelled result as an Excel file:
xlswrite('Relabelled.xls', Transformed)