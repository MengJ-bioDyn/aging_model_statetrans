% model for "Divergent Aging in Isogenic cells"
% state transition model with memory effect for cell 
% UCSD BioDynamics Lab, 2018

% make sure you download shadedErrorBar function from https://github.com/raacampbell/shadedErrorBar

% clc
clear all 

             
     para = [ 1   1,    0.015 -0.04, 0.008  -0.03, 0.012 -0.019,  ...  % passive 0-0 (-0.035 1.089), 0-1', 0-1, 0-2
              0   0,    0.022 0.22,  0.0125 -0.05, ...                 % passive 1'-0 (-0.033, 0.84),  1'-1',  1'-1,
             -0.02 0.4, -0.01  0.25, 0  0, ...                         % 1-0, 1-1', passive 1-1 (0.03 0.35)
              0.008 0.75];                                             % 2-2  (2-0 is 1-p22)
         
     para_memo =[0  0,  0.023  0,      0.009 -0.09,  0.004 -0.02, ...  % passive 0-0 (-0.036 1.11), 0-1', 0-1, 0-2
                 0  0,  0.005  -0.05,  0.007 -0.04,  0.036 0.1];       % passive 0-0 (-0.048 0.99), 0-1', 0-1, 0-2
       

         
a0 = para(1);
b0 = para(2);

a01s = para(3);        
b01s = para(4);       % p0-1s

a01 = para(5);        
b01 = para(6);        % p0-1

a02 = para(7);        % p0-2
b02 = para(8);         

                     
a1s0 = para(9);        
b1s0 = para(10);        

a1s1s = para(11);        
b1s1s = para(12);       % p1s-1s

a1s1 = para(13);        
b1s1 = para(14);        % p1s-1 

%
a10 = para(15);
b10 = para(16);        % p1-0 

a11s = para(17);
b11s = para(18);       % p1-1s

a11 = para(19);
b11 = para(20);        


a22 = para(21);
b22 = para(22);        % p22


% with history of S1/S1'
a0_w1 = para_memo(1);
b0_w1 = para_memo(2);

a01s_w1 = para_memo(3);        
b01s_w1 = para_memo(4);       

a01_w1 = para_memo(5);        
b01_w1 = para_memo(6);        

a02_w1 = para_memo(7);        
b02_w1 = para_memo(8); 

% with history of S1/S1'
a0_w2 = para_memo(9);
b0_w2 = para_memo(10);

a01s_w2 = para_memo(11);        
b01s_w2 = para_memo(12);       

a01_w2 = para_memo(13);        
b01_w2 = para_memo(14);        

a02_w2 = para_memo(15);        
b02_w2 = para_memo(16); 


mymap = [ 255 255 255; 
           205 205 205; 
           230 159 0; 
           222 45  38; 
           49,130,189;  ]/255;
                
s0_clr = [205 205 205]./255;
s1p_clr = [230 159 0]./255;
s1_clr = [222 45  38]./255;
s2_clr = [49,130,189]./255;

my_clr =mymap(2:5,:);


% Simulation time (generations)
Tmax = 80;

% experimental measured transition probability of N consecutive S1/S2 transit to death
% fit exp_transitionProb into a polynomial, max n=2
% Maximum number of consecutive s1 to test
max_Nsi = Tmax;

fit_transitionProb_s1 = 0.01*(1:max_Nsi).^2 - 0.036*(1:max_Nsi)+0.15;
fit_transitionProb_s2 = 0.009*(1:max_Nsi).^2 - 0.035*(1:max_Nsi)+0.12;


temp1 = find(fit_transitionProb_s1>1);
temp2 = find(fit_transitionProb_s2>1);

transitionProb_s1 = fit_transitionProb_s1;
transitionProb_s2 = fit_transitionProb_s2;

transitionProb_s1(temp1)=1;
transitionProb_s2(temp2)=1;


% Initial condition
% [s0 s0_a s0_b s1s s1 s2]
x0 = [1  0  0   0   0  0]';
% Simulation time
TPx = 0:Tmax;
%=========================================
% initiate Transition probability
% Transition Prob 0 to 0, 1, 2
TP01s = a01s*(0:Tmax) + b01s; TP01s((TP01s<0))=0;  
TP01 = a01*(0:Tmax) + b01;    TP01((TP01<0))=0;    
TP02 =  a02*(0:Tmax) + b02;   TP02((TP02<0))=0;    
TP00 = 1-TP01s - TP01 - TP02;
    neg_id = find(TP00<0,1,'first');
    if ~isempty(neg_id )
        TP00(neg_id:end)=0; 
        TP01s(neg_id:end)=TP01s(neg_id); 
        TP01(neg_id:end)=TP01(neg_id);
        TP02(neg_id:end)=TP02(neg_id);  
    end

TP01s_w1 = a01s_w1*(0:Tmax) + b01s_w1; TP01s_w1((TP01s_w1<0))=0; 
TP01_w1 = a01_w1*(0:Tmax) + b01_w1;    TP01_w1((TP01_w1<0))=0;   
TP02_w1 =  a02_w1*(0:Tmax) + b02_w1;   TP02_w1((TP02_w1<0))=0;   
TP00_w1 = 1- TP01s_w1 - TP01_w1 - TP02_w1;
neg_id = find(TP00_w1<0,1,'first');
    if ~isempty(neg_id)
        TP00_w1(neg_id:end)=0; 
        TP01s_w1(neg_id:end)=TP01s_w1(neg_id); 
        TP01_w1(neg_id:end)=TP01_w1(neg_id);
        TP02_w1(neg_id:end)=TP02_w1(neg_id);
    end

TP01_w2 = a01_w2*(0:Tmax) + b01_w2;    TP01_w2((TP01_w2<0))=0;   
TP01s_w2 = a01s_w2*(0:Tmax) + b01s_w2; TP01s_w2((TP01s_w2<0))=0; 
TP02_w2 =  a02_w2*(0:Tmax) + b02_w2;   TP02_w2((TP02_w2<0))=0;  TP02_w2((TP02_w2>1))=1; 
TP00_w2 = 1- TP01s_w2 - TP01_w2 - TP02_w2;
    neg_id = find(TP00_w2<0,1,'first');
    if ~isempty(neg_id)
        TP00_w2(neg_id:end)=0; 
        TP01s_w2(neg_id:end)=TP01s_w2(neg_id); 
        TP01_w2(neg_id:end)=TP01_w2(neg_id);
        TP02_w2(neg_id:end)=TP02_w2(neg_id);
    end

% S1s to S0, S1s, S1
TP1s1s = a1s1s*(0:Tmax) + b1s1s;   TP1s1s((TP1s1s<0))=0; 
TP1s1 =  a1s1*(0: Tmax) + b1s1;    TP1s1((TP1s1<0))=0;
TP1s0 =  1- TP1s1s - TP1s1;
    neg_id = find(TP1s0<0,1,'first');
    if ~isempty(neg_id) 
        TP1s0(neg_id:end) = 0;
        TP1s1s(neg_id:end)= TP1s1s(neg_id); 
        TP1s1(neg_id:end)=  TP1s1(neg_id);
    end


% S1 to S0, S1s, S1
TP10 = a10*(0:Tmax) + b10;     TP10((TP10<0))=0;    
TP11s = a11s*(0:Tmax) + b11s;  TP11s(TP11s<0)=0;  
TP11 = 1- TP10 - TP11s;
neg_id = find(TP10==0,1,'first');
    if ~isempty(neg_id) 
        TP10(neg_id:end) = 0;
        TP11s(neg_id:end)= TP11s(neg_id); 
        TP11(neg_id:end)=  TP11(neg_id);
    end

% S2 to S0, S2
TP22 = a22*(0:Tmax) + b22;
Lg_22 = find(TP22>1);
TP22(Lg_22)=1;
TP20 = 1- TP22;

%% Set up simulation (multiple realizations)
% number of runs
runNum=50;

% Number of trajectories to simulate per run (each trajectory can be treated as a cell)
numReals = 205; 

test =1;  % test 1: show each simulation. 0: not show.
deltaT = 1; % generation for step
Treg_lngth = ceil(Tmax/deltaT+1);

% RLS
check_RLS_D1 = nan(runNum, Treg_lngth);
check_RLS_D2 = nan(runNum, Treg_lngth);
check_RLS_all = nan(runNum, Treg_lngth);
% percentage
check_S0N = nan(runNum, Treg_lngth);
check_S1sN = nan(runNum, Treg_lngth);
check_S1N = nan(runNum, Treg_lngth);
check_S2N = nan(runNum, Treg_lngth);
% count number
check_S1count = nan(runNum, Treg_lngth);
check_S2count = nan(runNum, Treg_lngth);
% death percentage
check_S1D = nan(runNum, Treg_lngth);
check_S2D = nan(runNum, Treg_lngth);

% save Cell States
check_cellstate = cell(1,runNum);  % s0, s1s, s1, s2, Death
check_sortage = cell(3,runNum);  % sort_DT1, sort_DT2, sort_all

for JJ =1:runNum
% initialize rand generator
% reset(RandStream.getGlobalStream);

% Initial state and time
t = repmat(0,[1 numReals]);
x = repmat(x0,[1 numReals]);

% Number of steps done
idx = 0;
% Set up arrays to record trajectories

deltaT = 1; % generation for step
X = nan(size(x,1),ceil(Tmax/deltaT+1),size(x,2));
T = nan(size(X,2),size(x,2));
tsampleIdx = 1;
X(:,1,:) = x;
T(1,:) = t;
Treg = (0:size(X,2)-1)*deltaT;

   
      %  p0-0 p0-1s p0-1 p0-2, p0a-0a p0a-1s p0a-1 p0a-2, p0b-0b p0b-1s p0b-1 p0b-2, p1s-0a p1s-1s p1s-1, p1-0a p1-1s p1-1, p2-0b p2-2 
S = [       0   -1    -1   -1,    0      0     0      0,     0      0     0      0,    0      0     0,      0     0     0,    0    0   % s0
            0    0     0    0,    0     -1    -1     -1,     0      0     0      0,    1      0     0,      1     0     0,    0    0   % s0_a %with history of S1/S1'
            0    0     0    0,    0      0     0      0,     0     -1     -1    -1,    0      0     0,      0     0     0,    1    0   % s0_b %with history of S2
            0    1     0    0,    0      1     0      0,     0      1     0      0,    -1     0    -1,      0     1     0,    0    0   % s1s
            0    0     1    0,    0      0     1      0,     0      0     1      0,    0      0     1,     -1    -1     0,    0    0    % s1
            0    0     0    1,    0      0     0      1,     0      0     0      1,    0      0     0,      0     0     0,   -1    0 ]; % s2

% how long each traj has been in state1
s1_deltaT_sum  = zeros(1, numReals);
s2_deltaT_sum  = zeros(1, numReals);

% Simulation
while any(t<Tmax)

updateThese = t<Tmax;

p00 = TP00(t+1);
p01s = TP01s(t+1);
p01 = TP01(t+1);
p02 = TP02(t+1);

p00_w1 = TP00_w1(t+1);
p01s_w1 = TP01s_w1(t+1);
p01_w1 = TP01_w1(t+1);
p02_w1 = TP02_w1(t+1);

p00_w2 = TP00_w2(t+1);
p01s_w2 = TP01s_w2(t+1);
p01_w2 = TP01_w2(t+1);
p02_w2 = TP02_w2(t+1);

p1s0 = TP1s0(t+1);
p1s1s = TP1s1s(t+1);
p1s1 = TP1s1(t+1);

p10 = TP10(t+1);
p11s = TP11s(t+1);
p11 = TP11(t+1);

p20 = TP20(t+1);
p22 = TP22(t+1);


rates = [p00; p01s; p01; p02;  p00_w1; p01s_w1; p01_w1; p02_w1;  p00_w2; p01s_w2; p01_w2; p02_w2; ...
         p1s0; p1s1s; p1s1; p10; p11s; p11; p20; p22];

drawP = rand(1, numReals);
reaction = zeros(1, numReals);


state_id = cell(1,6);
csrates = cell(1,6);

for k=1:6
      % 1:4, 5:8, 9:12, 13:15, 16:18, 19:20
      if k<4
        csrates{k} = cumsum(rates(1+4*(k-1):4*k,:),1)./repmat(sum(rates(1+4*(k-1):4*k,:),1),[size(rates(1+4*(k-1):4*k,:),1) 1]); % 4 x numReals
      elseif k<6
        csrates{k} = cumsum(rates(13+3*(k-4):13+3*(k-4)+2,:),1)./repmat(sum(rates(13+3*(k-4):13+3*(k-4)+2,:),1), ...
              [size(rates(13+3*(k-4):13+3*(k-4)+2,:),1) 1]); % 3 x numReals
      elseif k==6
        csrates{k} = cumsum(rates(19:20,:),1)./repmat(sum(rates(19:20,:),1), ...
              [size(rates(19:20,:),1) 1]);
      end

      
    state_id{k} = find(x(k,:)==1);

    chooseReactMat{k} = repmat(drawP(state_id{k}), [size(csrates{k},1) 1])  <=  csrates{k}(:,state_id{k});
    indcs = repmat((1:size(chooseReactMat{k},1))',[1 size(chooseReactMat{k},2)]);
    indcs(~chooseReactMat{k}) = inf;
    if k<4
            reaction(state_id{k}) = min(indcs) + 4*(k-1);
    elseif k==4
            reaction(state_id{k}) = min(indcs)+12;
    elseif k==5
            reaction(state_id{k}) = min(indcs)+15;
    elseif k==6
            reaction(state_id{k}) = min(indcs)+18;
    end     
end


% 
% chooseReactMat0 = repmat(drawP(state0_id),[size(csrates_0,1) 1])<=csrates_0(:,state0_id);
% indcs0 = repmat((1:size(chooseReactMat0,1))',[1 size(chooseReactMat0,2)]);
% indcs0(~chooseReactMat0) = inf;
% reaction(state0_id) = min(indcs0,[],1);
% 
% chooseReactMat1s = repmat(drawP(state1s_id),[size(csrates_1s,1) 1])<=csrates_1s(:,state1s_id);
% indcs1s = repmat((1:size(chooseReactMat1s,1))',[1 size(chooseReactMat1s,2)]);
% indcs1s(~chooseReactMat1s) = inf;
% reaction(state1s_id) = min(indcs1s,[],1)+4;
% 
% chooseReactMat1 = repmat(drawP(state1_id),[size(csrates_1,1) 1])<=csrates_1(:,state1_id);
% indcs1 = repmat((1:size(chooseReactMat1,1))',[1 size(chooseReactMat1,2)]);
% indcs1(~chooseReactMat1) = inf;
% reaction(state1_id) = min(indcs1,[],1)+7;
% 
% chooseReactMat2 = repmat(drawP(state2_id),[size(csrates_2,1) 1])<=csrates_2(:,state2_id);
% indcs2 = repmat((1:size(chooseReactMat2,1))',[1 size(chooseReactMat2,2)]);
% indcs2(~chooseReactMat2) = inf;
% reaction(state2_id) = min(indcs2,[],1)+10;

   
deltaX = S(:,reaction(updateThese>0));

% only update those with t<Tmax
x(:,updateThese) = x(:, updateThese) + deltaX;
t(:,updateThese) = t(:,updateThese)+ deltaT;

tsampleIdx = tsampleIdx +1;

% Calculate time for consequtive s1 state for all trajectories.
% and generate random numbers to compare to the probability

s1_countL = double(x(5,:)==1);
s1_deltaT_sum = s1_deltaT_sum + deltaT.*s1_countL; 
s1_deltaT_sum = s1_deltaT_sum.*s1_countL;

% index of trajectories that we need to decide to die or not
s1_index = find(x(5,:)==1);
death_rand = rand(1, length(s1_index));
deathprob = transitionProb_s1(s1_deltaT_sum(s1_index));
death_true = find(death_rand - deathprob<0);
s1_death = s1_index(death_true);

s2_countL = double(x(6,:)==1);
s2_deltaT_sum = s2_deltaT_sum + deltaT.*s2_countL; 
s2_deltaT_sum = s2_deltaT_sum.*s2_countL;

% index of trajectories that we need to decide to die or not
s2_index = find(x(6,:)==1);
death_rand = rand(1, length(s2_index));
deathprob = transitionProb_s2(s2_deltaT_sum(s2_index));
death_true = find(death_rand - deathprob<0);
s2_death = s2_index(death_true);

% store current timestep
X(:, tsampleIdx, updateThese) = x(:, updateThese);
T(tsampleIdx, updateThese) = t(updateThese);


% asign "death" to s1_death trajectories 
% and set X to 2 for death state, t(s1_death) to Tmax, removing them from updateThese pool. 
if ~isempty(s1_death)
    for k=s1_death
        X(5, tsampleIdx+1:end,k) = 2;
        t(k) = Tmax;
    end
end

if ~isempty(s2_death)
    for k=s2_death
        X(6, tsampleIdx+1:end,k) = 2;
        t(k) = Tmax;
    end
end

end

species = cell(1, size(x,1));
for i=1:size(x,1)
    species{i} = squeeze(X(i,:,:));
end
S0 = species{1}';
S0_w1 = species{2}';
S0_w2 = species{3}';
S1s = species{4}';
S1 = species{5}';
S2 = species{6}';

S0_t = S0;
S0_w1_t = S0_w1;
S0_w2_t = S0_w2;
S1s_t = S1s;
S1_t = S1;
S2_t = S2;

S_origin = S1s + 2*S1 +3*S2; % let S1'=1, S1=2, S2=3, for matching colormap in exp data
S_origin_save = S_origin;

% save the states of each run, for conditional prob
for i = 1:numReals
    nan_id0 = find(isnan(S0(i,:)),1,'first');
    nan_id0_w1 = find(isnan(S0_w1(i,:)),1,'first');
    nan_id0_w2 = find(isnan(S0_w2(i,:)),1,'first');

    nan_id1s = find(isnan(S1s(i,:)),1,'first');
    nan_id1 = find(isnan(S1(i,:)),1,'first');
    nan_id2 = find(isnan(S2(i,:)),1,'first');

    if isempty(nan_id1)==0
        DT_tag =1;
        check_end = unique([nan_id0, nan_id0_w1, nan_id0_w2, nan_id1s, nan_id2]);
    elseif isempty(nan_id2)==0
        DT_tag =2;
        check_end = unique([nan_id0, nan_id0_w1, nan_id0_w2, nan_id1s, nan_id1]);
    end
    if length(check_end)>1
        disp(['warning! Run',num2str(numReals),' inconsistent.'])
    end
    
    S_origin_save(i,check_end:end)=-1; % change NaN into -1
end

S_origin_save2 = S_origin_save + 1*(S_origin_save==3); % change S2=3 to S2=4, match the experimental colormap
check_cellstate{JJ} = S_origin_save2;

% sort S based on its lifespan
death_time = zeros(1,numReals);
for i=1:numReals
    nan_id = find(isnan(T(:,i)),1,'first');
    death_time(i) = Treg(nan_id);
    S_origin_save(i,nan_id:end)=-1;
end

[sort_death_time sort_death_id ]= sort(death_time );
S_sorted = S_origin(sort_death_id,:); 

D1_id = find(S1(:,end)==2);
D2_id = find(S2(:,end)==2);
[sort_death_time_D1, sort_D1only_id] = sort(death_time(D1_id));
[sort_death_time_D2, sort_D2only_id] = sort(death_time(D2_id));

sort_DT1_id = D1_id(sort_D1only_id);
sort_DT2_id = D2_id(sort_D2only_id);

% save traj_id order
check_sortage{1,JJ} = sort_DT1_id;
check_sortage{2,JJ} = sort_DT2_id;
check_sortage{3,JJ} = sort_death_id;







%% statistic check



% check 1. state
StateM = nan(7, length(Treg));
StateM(1,:) = sum((S0==1),1)+ sum((S0_w1==1),1)+ sum((S0_w2==1),1);
StateM(2,:) = sum((S1s==1),1);
StateM(3,:) = sum((S1==1),1);
StateM(4,:) = sum((S2==1),1);
StateM(5,:) = sum((S1==2),1); % death 1
StateM(6,:) = sum((S2==2),1); % death 2
StateM(7,:) = StateM(5,:) + StateM(6,:);
all_state = numReals*ones(7, length(Treg));
StateP = StateM./all_state;
if runNum==1
    figure; plot(StateP([1:4,7],:)')
else
    check_S0N(JJ,:) = StateP(1,:);
    check_S1sN(JJ,:)= StateP(2,:);
    check_S1N(JJ,:) = StateP(3,:);
    check_S2N(JJ,:) = StateP(4,:);
    check_S1D(JJ,:) = StateP(5,:);
    check_S2D(JJ,:) = StateP(6,:);
    
end

% check 2. RLS
RLS_1 = 1- StateM(5,:)./max(StateM(5,:));
RLS_2 = 1- StateM(6,:)./max(StateM(6,:));
RLS_all = 1- StateP(7,:);
if runNum==1
    figure; hold on;
    plot(Treg, RLS_1,'ro-', Treg, RLS_2,'bo-',Treg, RLS_all,'ko-');
else
    check_RLS_D1(JJ,:) = RLS_1;
    check_RLS_D2(JJ,:) = RLS_2;
    check_RLS_all(JJ,:) = RLS_all;
end

% check 3. consecutive S1/2 before death.
if runNum==1
    figure; hold on;
    plot(StateM(2,:), 'g+-')
    plot(StateM(3,:), 'r+-')
    plot(StateM(4,:), 'b+-')
else
    check_S1count(JJ,:) = StateM(3,:);
    check_S2count(JJ,:) = StateM(4,:);
end


end


%% check RLS

load 'WT_sanitycheck_5S_pub.mat' % data from WT experiments


simu_RLS_D1_avg = mean(check_RLS_D1, 1);
simu_RLS_D1_std = std(check_RLS_D1, 1);
simu_RLS_D2_avg = mean(check_RLS_D2, 1);
simu_RLS_D2_std = std(check_RLS_D2, 1);
simu_RLS_all_avg = mean(1-check_RLS_all, 1);
simu_RLS_all_std = std((1-check_RLS_all), 1);
simu_S0N_avg = mean(check_S0N,1);
simu_S0N_std = std(check_S0N,1);

figure; 
subplot(2,2,1); hold on; box on;
% WT data
plot(1:length(lifespan_DT1),lifespan_DT1,'.','Color',s1_clr,'markersize',15 ); 
plot(1:length(lifespan_DT2),lifespan_DT2,'.','Color',s2_clr,'markersize',15 ); 
% simulation
shadedErrorBar(0:Tmax, simu_RLS_D1_avg, simu_RLS_D1_std,'lineprops',{'-','Color',s1_clr },'patchSaturation',0.33)
shadedErrorBar(0:Tmax, simu_RLS_D2_avg, simu_RLS_D2_std,'lineprops',{'-','Color',s2_clr },'patchSaturation',0.33)

set(gca,'fontsize',16); xlim([0 60]); ylim([0 1])
title('RLS for Path1 and Path2'); 
xlabel('Division number'); 

ylabel('Survival rate');

subplot(2,2,3); hold on; box on;
% WT data
plot(hist_state([1],:)','.', 'Color',s0_clr ,'markersize',15);
plot(hist_state([5],:)','.','Color','Black','markersize',15);
% simulation
shadedErrorBar(0:Tmax, simu_S0N_avg, simu_S0N_std,'lineprops',{'-','Color',s0_clr },'patchSaturation',0.33)
shadedErrorBar(0:Tmax, simu_RLS_all_avg, simu_RLS_all_std,'lineprops',{'-','Color',[63/255,63/255,63/255]},'patchSaturation',0.33)

set(gca,'fontsize',16); xlim([0 60]); 
title('Population change of S0 and Death (both paths)'); 
xlabel('Division number');
ylabel('Percentage');
% legend('exp RLS DT1','exp RLS DT2','exp S0','exp RLS all','simu RLS DT1','simu RLS DT2','simu S0','simu RLS all')
xlim([0 60]); ylim([-0.01 1.01]); box on;


% check state percentage

simu_S1sN_avg = mean(check_S1sN,1);
simu_S1sN_std = std(check_S1sN,1);
simu_S1N_avg = mean(check_S1N,1);
simu_S1N_std = std(check_S1N,1);
simu_S2N_avg = mean(check_S2N,1);
simu_S2N_std = std(check_S2N,1);

subplot(2,2,2);
hold on; box on;
% WT data
plot(hist_state([2],:)','.','Color',s1p_clr ,'markersize',15)
plot(hist_state([3],:)','.','Color',s1_clr ,'markersize',15)
plot(hist_state([4],:)','.','Color',s2_clr ,'markersize',15)
% simulation
shadedErrorBar(0:Tmax, simu_S1sN_avg, simu_S1sN_std,'lineprops',{'-','Color',s1p_clr },'patchSaturation',0.33)
shadedErrorBar(0:Tmax, simu_S1N_avg, simu_S1N_std,'lineprops',{'-','Color',s1_clr },'patchSaturation',0.33)
shadedErrorBar(0:Tmax, simu_S2N_avg, simu_S2N_std,'lineprops',{'-','Color',s2_clr },'patchSaturation',0.33)

set(gca,'fontsize',16);
xlabel('Division number'); ylabel('Percentage');
title('Population change of S1^\prime, S1 and S2')
xlim([0 60]); ylim([-0.01 0.4]); 
set(gcf, 'position',[50 100 800 640])
       

figure; 
subplot(1,2,1);
imagesc(S_origin(sort_DT1_id,:)); xlim([0 Tmax]);set(gca,'fontsize',14);
title('One set of realization, path 1 subset'); ylabel('single cell trajectory');
subplot(1,2,2);
imagesc(S_origin(sort_DT2_id,:)); xlim([0 Tmax]);
set(gca,'fontsize',14);   colormap(my_clr);
xlabel('time, generations'); 
title('One set of realization, path 2 subset');