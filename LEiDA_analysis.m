%% LEiDA_analysis.m
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
%
% This function clusters and analyses BOLD data using LEiDA
% Comparing between 2 groups ('low' and 'high')
%
% 1 - Compute the leading eigenvector V1 and saves them in LEiDA_data.mat
% 2 - Clusters the Leading Eigenvectors into recurrent PL states
% 3 - Analyzes the percentage of occupancy and lifetime of each PL state
% 4 - Calculates statistical signigifance between groups
%   - Saves Clusters and statistics into LEiDA_results.mat
% 5 - Plots PL states and group differences for a given clustering solution
%
%
% Code from Joana Cabral March 2018
% joana.cabral@psych.ox.ac.uk
% adapted by alonsomartinez.sonsoles@gmail.com March 2020

clc
clearvars

%% 1 - Compute the Leading Eigenvectors from the BOLD datasets

if ~exist('LEiDA_data.mat','file')
    
    % USER: Provide time series for each subject
    load timeseries_data.mat ts 
    % ts is a cell of length=n_Subjects [number of subjects]. In each cell
    % timeseries are stored in a matrix of size=N_areas*Tmax [number of
    % regions*volumes]
    [n_Subjects, ~]=size(ts);
    [N_areas, Tmax]=size(ts{1,1});
    TR=2;

    disp('Processing the eigenvectors from BOLD data') 
    % Preallocate variables to save PL patterns and associated information
    V1_all   = zeros((Tmax-2)*n_Subjects,N_areas); % All leading eigenvectors
    t_all=0; % Index of time (starts at 0 and will be updated until n_Subjects*(Tmax-2))
    Time_all= zeros((Tmax-2)*n_Subjects,1); % Vector that links each frame to a subject

    % Bandpass filter settings
    fnq=1/(2*TR);                 % Nyquist frequency
    flp = 0.01;                   % lowpass frequency of filter (Hz)
    fhi = 0.1;                    % highpass
    Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
    k=9;                          % 9nd order butterworth filter %USER: adapt according to imaging data
    [bfilt,afilt]=butter(k,Wn);   % construct the filter
    clear fnq flp fhi Wn k

    for s=1:n_Subjects

        % USER: Adapt to load here the BOLD matrix (NxT) from each scan
        BOLD = ts{s,1};

        % Get the BOLD phase using the Hilbert transform
        Phase_BOLD=zeros(N_areas,Tmax);
        for seed=1:N_areas
            ts = demean(detrend(BOLD(seed,:)));
            signal_filt =filtfilt(bfilt,afilt,ts);
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end


        % Slide over time discarding the first and last epochs
        for t=2:Tmax-1

            % Calculate the Instantaneous BOLD PL (Phase-Locking)
            iFC=zeros(N_areas);
            for n=1:N_areas
                for p=1:N_areas
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end

            % Get the leading eigenvector
            [V1,~]=eigs(iFC,1);
            if sum(V1)>0
                V1=-V1;
            end

            % Save V1 from all frames in all fMRI sessions
            t_all=t_all+1; % Update time
            V1_all(t_all,:)=V1;
            Time_all(t_all)=s;
        end
    end
    % save
%     save LEiDA_data.mat N_areas TR Tmax Time_all V1_all...
%         n_Subjects Index_low n_Low Index_high n_High
else
    % Load the Leading Eigenvectors from the previous step
    load LEiDA_data.mat N_areas TR Tmax Time_all V1_all...
        n_Subjects Index_low n_Low Index_high n_High
end


%%  2- Cluster the Leading Eigenvectors

disp('Clustering the eigenvectors into')
% V1_all is a matrix containing all the eigenvectors:
% Collumns: N_areas are brain areas (variables)
% Rows: (Tmax-2)*n_Subjects are all time points (independent observations)

% USER: Set maximum/minimum number of clusters
% There is no fixed number of PL states that the brain can display
% Keep the range small for the first trials
% Extend the range depending on the hypothesis of each work
maxk=14;
mink=4;
rangeK=mink:maxk;

% Set the parameters for Kmeans clustering
Kmeans_results=cell(size(rangeK));
for k=1:length(rangeK)
    disp(['- ' num2str(rangeK(k)) ' PL states'])
    [IDX, C, SUMD, D]=kmeans(V1_all,rangeK(k),'Distance','cosine',...
        'Replicates',100,'MaxIter',1000,'Display','off');
    [~, ind_sort]=sort(hist(IDX,1:rangeK(k)),'descend');
    [~,idx_sort]=sort(ind_sort,'ascend');
    Kmeans_results{k}.IDX=idx_sort(IDX);   % Cluster time course - numeric collumn vectors
    Kmeans_results{k}.C=C(ind_sort,:);     % Cluster centroids (PL patterns) V
    Kmeans_results{k}.SUMD=SUMD(ind_sort); % Within-cluster sums of point-to-centroid distances
    Kmeans_results{k}.D=D(:,ind_sort);     % Distance from each point to every centroid 
end

%% 3 - Calculate %occupancy and lifetimes of each PL state c.

% Preallocate variables
P=zeros(n_Subjects,maxk-mink+1,maxk);
LT=zeros(n_Subjects,maxk-mink+1,maxk);

for k=1:length(rangeK) % Cluster solution k
    for s=1:n_Subjects

        % Select the time points representing this subject
        T=(Time_all==s);
        Ctime=Kmeans_results{k}.IDX(T);

        for c=1:rangeK(k)
            % Percentage of occupancy
            P(s,k,c)=mean(Ctime==c)*100;

            % Mean Lifetime
            Ctime_bin=Ctime==c;

            % Detect switches in and out of this state
            a=find(diff(Ctime_bin)==1);
            b=find(diff(Ctime_bin)==-1);

            % We discard the cases where state sarts or ends ON
            if length(b)>length(a)
                b(1)=[];
            elseif length(a)>length(b)
                a(end)=[];
            elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                b(1)=[];
                a(end)=[];
            end
            if ~isempty(a) && ~isempty(b)
                C_Durations=b-a;
            else
                C_Durations=0;
            end
            LT(s,k,c)=mean(C_Durations)*TR;
        end
    end
end

%% 4 - Statistical significance between groups

% Preallocate variables
disp('Testing differences in State Occupancy between Low and High MDI')
P_pval=zeros(length(rangeK),max(rangeK));
P_pval_bh=zeros(length(rangeK),max(rangeK));
for k=1:length(rangeK)
    disp(['Now running statistics for ' num2str(rangeK(k)) ' PL states'])
    for c=1:rangeK(k)
        % Compare Probabilities
        a=squeeze(P(Index_low,k,c))';  % Vector containing Prob of c in Low MDI
        b=squeeze(P(Index_high,k,c))';  % Vector containing Prob of c in High MDI
        stats=permutation_htest2_np_gender([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest',Gender);
        P_pval(k,c)=min(stats.pvals);
    end
    % FDR corrections
    [~,~,~, P_pval_bh(k,1:rangeK(k))]=fdr_bh(P_pval(k,1:rangeK(k)),0.05,'pdep','yes');   
end

% Preallocate variables
disp('Testing differences in State Lifetimes between between Low and High MDI')
LT_pval=zeros(length(rangeK),max(rangeK));
LT_pval_bh=zeros(length(rangeK),max(rangeK));
for k=1:length(rangeK)
    disp(['Now running statistics for ' num2str(rangeK(k)) ' PL states'])
    for c=1:rangeK(k)
        % Compare Probabilities
        a=squeeze(LT(Index_low,k,c))';  % Vector containing Prob of c in Controls
        b=squeeze(LT(Index_high,k,c))';  % Vector containing Prob of c in Patients
        stats=permutation_htest2_np_gender([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest',Gender);
        LT_pval(k,c)=min(stats.pvals);
    end
     % FDR corrections
     [~,~,~, LT_pval_bh(k,1:rangeK(k))]=fdr_bh(LT_pval(k,1:rangeK(k)),0.05,'pdep','yes'); 
end

% Save results
% save LEiDA_results.mat V1_all Time_all Kmeans_results P LT P_pval LT_pval...
%     P_pval_bh LT_pval_bh rangeK maxk mink n_Low n_High N_areas TR Tmax...
%     n_Subjects Index_low Index_high


%% 5 - Plot PL patterns for a cluster solution of k=9

disp('%%% PLOTS %%%%')

load LEiDA_results.mat
load power223_symm.mat Order netID netName roiName
color_net ={'skyblue','skyblue','purple','pink','red', 'grey','blue','yellow','black','brown','teal','green','tan','lightgrey'};

K = 9;
load mymap.mat mymap % colors for 9 clusters
Clusters=Kmeans_results{rangeK==K};
k=find(rangeK==K);
V=Clusters.C;

close all
figure('Name','Repertoire_of_9_PL_states','color','w')
colormap(jet)
location = 1:2:K*7*2;

for c=1:K
    
    % Pannel A - Plots the PL states over the cortex
    subplot(7,K*2,[location(c),location(c)+1])
    plot_nodes_in_cortex_new_223(V(c,:),mymap(c,:))
    title(['PL state ' num2str(c)],'fontsize',10);
   
    % Pannel B - Plots V1 in matrix format
    subplot(7,K*2,[location(c+K),location(c+K)+1])
    FC_V=V(c,Order)'*V(c,Order);
    li=max(abs(FC_V(:)));
    imagesc(FC_V,[-li li])
    axis square
    set(gca,'YTickLabel',[],'XTickLabel',[])
    title('PL matrix','FontWeight','normal')
    ylabel('Brain region #')
    xlabel('Brain region #')
   
    % Pannel C - Plots V1 in vector format
    subplot(7,K*2,[location(c+2*K):location(c+2*K)+1, location(c+3*K):location(c+3*K)+1, location(c+4*K):location(c+4*K)+1])
    Vc=V(c,Order);
    Vc=Vc/max(abs(Vc));
    hold on
    for nn = 1:14   
        barh((Vc.*(netID==nn)'),'FaceColor',rgb(sprintf('%s',color_net{nn})),'EdgeColor','none','Barwidth',.5)
        hold on
    end
    ylim([0 N_areas+1])
    xlim([-1 1])
    set(gca,'YTickLabel',[])
    if c==1
     ylabel('brain region #')
    end
    pos = get(gca,'Position');
    pos(3)=0.04;
    pos(1)=pos(1)+0.0176;
    set(gca,'Position',pos);
    box on
   
    % Pannel D - Plots the Percentage of occupancy and lifetime of each state in each group
    for CAS = 1:2
        if CAS == 1
            CASE = P;
            P_case = P_pval;
            Pbh_case = P_pval_bh;
            CASELABEL = 'occupancy(%)'; 
        elseif CAS == 2
            CASE = LT;
            P_case = LT_pval;
            Pbh_case = LT_pval_bh;
            CASELABEL = 'lifetime(s)';
        end
        subplot(7,K*2,location(c+(CAS+4)*K):location(c+(CAS+4)*K)+1)
        Low=squeeze(CASE(Index_low,k,c));
        High=squeeze(CASE(Index_high,k,c));
        bar([1 3],[mean(Low) mean(High)],'EdgeColor','w','FaceColor',[.8 .8 .8]);hold on
        scatter(ones(n_Low,1),Low,15,'.','MarkerFaceColor',[0 .7 .7],'jitter', 'on', 'jitterAmount', 0.3);hold on;
        scatter(ones(n_High,1)*3,High,15,'.','MarkerFaceColor',[1 0.5 0.2],'jitter', 'on', 'jitterAmount', 0.3); hold on
        set(gca,'XTick',[1 3],'XTickLabel',{'Low';'High'}); hold on;

        % Error bar containing the standard error of the mean
        errorbar([1 3],[mean(Low) mean(High)],...
                 [std(Low)/sqrt(numel(Low)) std(High)/sqrt(numel(High))],...
                  'LineStyle','none','Color','k','linewidth',1)
        lim = get(gca,'YTick');
        lim = lim(end);
        if Pbh_case(k,c)<=0.05 && Pbh_case(k,c)>0
                text(1.8,lim,'**','fontsize',12)
        end
        if c==1
             ylabel(CASELABEL,'fontsize',12);
        end
        if CAS == 1
             set(gca,'XTick',[1 3],'XTickLabel',{''});
        end
        box off        
    end
end
xlabel('MDI')
%savefig('Figure_Repertoire_of_9_PL_states.fig')