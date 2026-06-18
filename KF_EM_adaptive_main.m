clear 
close all




% Load ECG data from MAT file
[file, path] = uigetfile('*.mat','Select ECG mat file');
data = load([path '\' file]);


fs = data.fs;                % Sampling frequency of ECG signal (Hz)

x = data.x;               % Extract stored signal matrix

dt = 1/fs; % time step

max_number_of_iterations = 15; % Maximum Number of EM Iterations 

[warnMsg, warnId] = lastwarn;
if ~isempty(warnMsg)
    lastwarn('')
end

%% adding white gaussian Noise
SNR = 6;   % you can change the Noise SNR level
y = awgn(x,SNR,'measured') ;
x_noisy =y;

[qrs_positions] = pantompkins_qrs(x_noisy,fs);

figure(1),plot(x_noisy,'b'),hold on,plot(qrs_positions,x_noisy(qrs_positions),'*r'),hold off
legend({'Noisy ECG Signal','R Peaks'})
title(file)

ind_Rpeak = qrs_positions;
%% finding Samples belonging to QRS
ind_qrs = [];
for i=1:length(ind_Rpeak)
    ind_qrs = [ind_qrs   ind_Rpeak(i)-round(0.03*fs):ind_Rpeak(i)+round(0.07*fs)];
end
ind_qrs = ind_qrs(ind_qrs>0);
ind_qrs = ind_qrs(ind_qrs<=length(x));
%% finding Samples belonging to P and T waves 
total_ind = 1:length(x);
ind_P_T = total_ind(~ismember(total_ind,ind_qrs));

% modification for colored noise
Psi_pt = 0.001;
Psi_qrs = 0.01;

y_new_qrs = y(ind_qrs);
y_new_pt = y(ind_P_T);
length_y_pt = size(y_new_pt,2);
length_y_qrs = size(y_new_qrs,2);

%% Creating two seperate signals from the original ECG
% one signal (y_new_qrs) is result of concatenation the samples belonging to QRS waves
% the other one (y_new_pt) is result of concatenation the samples belonging to P and T waves
length_y = size(y,2);
for i=1:length_y_pt-1
    y_new_pt(i) = y_new_pt(i+1)-Psi_pt*y_new_pt(i);
end

for i=1:length_y_qrs-1
    y_new_qrs(i) = y_new_qrs(i+1)-Psi_qrs*y_new_qrs(i);
end



% figure(1),plot(y)


R_pt = 15;%%; Initial value for measurement Covariance of y_new_pt
R_qrs = 1.5;%% Initial value for measurement Covariance of y_new_qrs


Q_pt = 1*eye(14); %% Initial value for state Covariance of y_new_pt
Q_qrs = 1*eye(14);%% Initial value for state Covariance of y_new_qrs


A_pt = eye(14)+1e-1*rand;%% Initial value for state transition matrix of y_new_pt
A_pt= triu(A_pt);
A_pt = eye(size(A_pt))+(1-eye(size(A_pt))).*A_pt;

A_qrs = eye(14)+1e-5*rand;%% Initial value for state transition matrix of y_new_qrs
A_qrs= triu(A_qrs);
A_qrs = eye(size(A_qrs))+(1-eye(size(A_qrs))).*A_qrs;


H_pt = zeros(1,14); %% Initial value for measurement  matrix of y_new_pt
H_pt(1,1) = 1;
H_qrs = zeros(1,14);%% Initial value for measurement matrix of y_new_qrs
H_qrs(1,1) = 1;

length_x_pt = size(A_pt,1);
length_x_qrs = size(A_qrs,1);

%% creating Square root versions of matrices necessary for iterative Kalman EM adaptation

B_pt = 1*eye(size(Q_pt));
B_qrs = 1*eye(size(Q_qrs));

H_new_pt = H_pt*A_pt-Psi_pt*H_pt;
H_new_qrs = H_qrs*A_qrs-Psi_qrs*H_qrs;



S_pt = Q_pt*B_pt'*H_pt';
S_qrs = Q_qrs*B_qrs'*H_qrs';

%%%
Q_s_pt = Q_pt.^0.5; % square_root of Q (Q^1/2)
R_s_pt = R_pt.^0.5; % square_root of R (R^1/2)
Q_s_qrs = Q_qrs.^0.5; % square_root of Q (Q^1/2)
R_s_qrs = R_qrs.^0.5; % square_root of R (R^1/2)
%R_new = H*B*Q*B'*H'+ R;
[~,R1_pt] = (qr([Q_s_pt*B_pt'*H_pt';R_s_pt]));
[~,R1_qrs] = (qr([Q_s_qrs*B_qrs'*H_qrs';R_s_qrs]));

R_new_s_pt = R1_pt(1:size(y,1),1:size(y,1))';
R_new_s_qrs = R1_qrs(1:size(y,1),1:size(y,1))';

R_new_pt = R_new_s_pt*R_new_s_pt';
R_new_qrs = R_new_s_qrs*R_new_s_qrs';

inv_R_new_s_pt = inv(R_new_s_pt);
inv_R_new_s_qrs = inv(R_new_s_qrs);

J_pt = B_pt*S_pt/R_new_pt;
J_qrs = B_qrs*S_qrs/R_new_qrs;

%Q_new = B*(Q+S/(R_new)*S')*B';
[~,R1_pt] = (qr([(inv_R_new_s_pt)'*S_pt';Q_s_pt]));
[~,R1_qrs] = (qr([(inv_R_new_s_qrs)'*S_qrs';Q_s_qrs]));

Q_new_s_pt = R1_pt(1:length_x_pt,1:length_x_pt)';
Q_new_s_qrs = R1_qrs(1:length_x_qrs,1:length_x_qrs)';

Q_new_pt = Q_new_s_pt*Q_new_s_pt';
Q_new_qrs = Q_new_s_qrs*Q_new_s_qrs';


A_new_pt = A_pt-J_pt*H_new_pt;
A_new_qrs = A_qrs-J_qrs*H_new_qrs;






%% Initialization of Kalman state vectors and covariance matrices

x_updat_pt = zeros(length_x_pt,length_y_pt);                %(x_k|k)
x_pred_k_k_1_pt= zeros(length_x_pt,length_y_pt);                  %(x_k|k-1)
x_pred_k_plus_k_pt= zeros(length_x_pt,length_y_pt);                  %(x_k+1|k)

x_smoothed_pt= zeros(length_x_pt,length_y_pt);              %(x_k|N)
Ps_pred_k_k_1_pt = cell(1,length_y_pt);                  %(P_k|k-1)^(1/2)
Ps_pred_k_plus_k_pt = cell(1,length_y_pt);               %(P_k+1|k)^(1/2)
Ps_pred_k_plus_k_pt{1,1} = 3*eye(length_x_pt,length_x_pt);  %(P_k+1|k)^(1/2)

Ps_updat_pt =  cell(1,length_y_pt);                      %(P_k|k)
Ps_smoothed_pt = cell(1,length_y_pt);                    %(P_k|N)

P_pred_k_k_1_pt = cell(1,length_y_pt);                  %(P_k|k-1)
P_pred_k_plus_k_pt = cell(1,length_y_pt);               %(P_k+1|k)
P_pred_k_plus_k_pt{1,1} = 3*eye(length_x_pt,length_x_pt);  %(P_k+1|k)

P_updat_pt =  cell(1,length_y_pt);                      %(P_k|k)
P_smoothed_pt = cell(1,length_y_pt);                    %(P_k|N)


P_s0_pt =  eye(length_x_pt);diag([1 .3 .3]);

X0_pt = zeros(length_x_pt,1);

x_k_1_k_2_pt = X0_pt;    %x(k-1)|(k-2)

y_new_k_1_pt = 0*y(:,1); % z_new(k-1)
% square_root_filtering

log_new_pt =0;
log_previous_pt = 0;
log_liklihood_pt = 0;


x_updat_qrs = zeros(length_x_qrs,length_y_qrs);                %(x_k|k)
x_pred_k_k_1_qrs= zeros(length_x_qrs,length_y_qrs);                  %(x_k|k-1)
x_pred_k_plus_k_qrs= zeros(length_x_qrs,length_y_qrs);                  %(x_k+1|k)

x_smoothed_qrs= zeros(length_x_qrs,length_y_qrs);              %(x_k|N)
Ps_pred_k_k_1_qrs = cell(1,length_y_qrs);                  %(P_k|k-1)^(1/2)
Ps_pred_k_plus_k_qrs = cell(1,length_y_qrs);               %(P_k+1|k)^(1/2)
Ps_pred_k_plus_k_qrs{1,1} = 3*eye(length_x_qrs,length_x_qrs);  %(P_k+1|k)^(1/2)

Ps_updat_qrs =  cell(1,length_y_qrs);                      %(P_k|k)
Ps_smoothed_qrs = cell(1,length_y_qrs);                    %(P_k|N)

P_pred_k_k_1_qrs = cell(1,length_y_qrs);                  %(P_k|k-1)
P_pred_k_plus_k_qrs = cell(1,length_y_qrs);               %(P_k+1|k)
P_pred_k_plus_k_qrs{1,1} = 3*eye(length_x_qrs,length_x_qrs);  %(P_k+1|k)

P_updat_qrs =  cell(1,length_y_qrs);                      %(P_k|k)
P_smoothed_qrs = cell(1,length_y_qrs);                    %(P_k|N)


P_s0_qrs =  eye(length_x_qrs);diag([1 .3 .3]);

X0_qrs = zeros(length_x_qrs,1);

x_k_1_k_2_qrs = X0_qrs;    %x(k-1)|(k-2)

y_new_k_1_qrs = 0*y(:,1); % z_new(k-1)
% square_root_filtering



for iteration = 1:max_number_of_iterations
    
    % square_root_filtering
    %     disp('log_new-log_previous')
    %     disp(num2str(log_new-log_previous))
    %     log_liklihood_sum = 0;
    %     inov_vec = zeros(1,300);
    %     denum_vec = zeros(1,300);
    %     lambda = .8;
    
    
    for i=1:length_y_pt
        
        % prediction
        %
        % x(k)|(k-1) = Phi_new*x(k-1)|(k-1)  + J*(z_new(k-1)-H_new*x(k-1)|(k-2))
        
        % P(k)|(k-1) = Phi_new*P(k-1)|(k-1)*Phi_new'+Q_new
        %
        
        
        x_pred_k_k_1_pt(:,i) = A_new_pt*X0_pt+J_pt*(y_new_k_1_pt-H_new_pt*x_k_1_k_2_pt);
        x_k_1_k_2_pt = x_pred_k_k_1_pt(:,i);
        
        
        
        
        % P_pred_k_k_1{i} = A_new*P0*A_new'+Q_new;
        
        
        
        [~,R1_pt] = (qr([P_s0_pt'*A_new_pt';Q_new_s_pt]));
        R1_pt = R1_pt(1:length_x_pt,1:length_x_pt);
        Ps_pred_k_k_1_pt{i} = R1_pt';
        
        P_pred_k_k_1_pt{i} = Ps_pred_k_k_1_pt{i}*Ps_pred_k_k_1_pt{i}';
        
        
        
        
        
        % update
        %
        % K = P(k)|(k-1)*H_new'*inv(H_new*P(k)|(k-1)*H_new'+R_new)
        % x(k)|(k) = x(k)|(k-1) + K*(z_new(k)-H_new*x(k)|(k-1))
        % P(k)|(k) =(I - K*H_new)*P(k)|(k-1)
        
        %K = P_pred_k_k_1{i}*(H_new')/(H_new*P_pred_k_k_1{i}*H_new'+R_new);
        K_pt = ((R1_pt'*R1_pt)*H_new_pt')/(H_new_pt*(R1_pt'*R1_pt)*H_new_pt'+R_new_pt);
        K_update_pt{i} = K_pt;
        
        inov_pt = y_new_pt(:,i)-H_new_pt*x_pred_k_k_1_pt(:,i);
        X0_pt = x_pred_k_k_1_pt(:,i) +K_pt*(inov_pt);
        x_updat_pt(:,i) = X0_pt;
        y_new_k_1_pt = y_new_pt(:,i);
        
        x_pred_k_plus_k_pt(:,i) = A_new_pt*X0_pt+J_pt*(y_new_k_1_pt-H_new_pt*x_k_1_k_2_pt);
        
        %         inov_vec = [inov inov_vec(1:end-1)];
        %         denum_vec = [H_new*(R1'*R1)*H_new'+R_new denum_vec(1:end-1)];
        
        
        % P0 =  (eye(size(Q))-K*H_new)*P_pred_k_k_1{i};
        % P_updat{1,i} = (P0);
        
        R2_pt = qr([R_new_s_pt' zeros(size(y,1),length_x_pt);R1_pt*H_pt' R1_pt]);
        P_s0_pt =  (R2_pt(end-length_x_pt+1:end,end-length_x_pt+1:end))';
        Ps_updat_pt{1,i} = (P_s0_pt);
        
        P_updat_pt{1,i} = Ps_updat_pt{1,i}*Ps_updat_pt{1,i}';
        
        %         if i>300 && iter==25
        %             R_new = lambda*R_new+mean((1-lambda)*(inov_vec.^2)./denum_vec);
        %             R_new_s = R_new^.05;
        %         end
        
        % P_pred_k_plus_k{i} = A_new*P0*A_new'+Q_new;
        [~,R3_pt] = (qr([P_s0_pt'*A_new_pt';Q_new_s_pt]));
        R3_pt = R3_pt(1:length_x_pt,1:length_x_pt);
        Ps_pred_k_plus_k_pt{i} = R3_pt';
        
        P_pred_k_plus_k_pt{i} = Ps_pred_k_plus_k_pt{i}*Ps_pred_k_plus_k_pt{i}';
        
        
        if i<=length_y_qrs
            % prediction
            %
            % x(k)|(k-1) = Phi_new*x(k-1)|(k-1)  + J*(z_new(k-1)-H_new*x(k-1)|(k-2))
            
            % P(k)|(k-1) = Phi_new*P(k-1)|(k-1)*Phi_new'+Q_new
            %
            
            
            x_pred_k_k_1_qrs(:,i) = A_new_qrs*X0_qrs+J_qrs*(y_new_k_1_qrs-H_new_qrs*x_k_1_k_2_qrs);
            x_k_1_k_2_qrs = x_pred_k_k_1_qrs(:,i);
            
            
            
            
            % P_pred_k_k_1{i} = A_new*P0*A_new'+Q_new;
            
            
            
            [~,R1_qrs] = (qr([P_s0_qrs'*A_new_qrs';Q_new_s_qrs]));
            R1_qrs = R1_qrs(1:length_x_qrs,1:length_x_qrs);
            Ps_pred_k_k_1_qrs{i} = R1_qrs';
            
            P_pred_k_k_1_qrs{i} = Ps_pred_k_k_1_qrs{i}*Ps_pred_k_k_1_qrs{i}';
            
            
            
            
            
            % update
            %
            % K = P(k)|(k-1)*H_new'*inv(H_new*P(k)|(k-1)*H_new'+R_new)
            % x(k)|(k) = x(k)|(k-1) + K*(z_new(k)-H_new*x(k)|(k-1))
            % P(k)|(k) =(I - K*H_new)*P(k)|(k-1)
            
            %K = P_pred_k_k_1{i}*(H_new')/(H_new*P_pred_k_k_1{i}*H_new'+R_new);
            K_qrs = ((R1_qrs'*R1_qrs)*H_new_qrs')/(H_new_qrs*(R1_qrs'*R1_qrs)*H_new_qrs'+R_new_qrs);
            K_update_qrs{i} = K_qrs;
            
            inov_qrs = y_new_qrs(:,i)-H_new_qrs*x_pred_k_k_1_qrs(:,i);
            X0_qrs = x_pred_k_k_1_qrs(:,i) +K_qrs*(inov_qrs);
            x_updat_qrs(:,i) = X0_qrs;
            y_new_k_1_qrs = y_new_qrs(:,i);
            
            x_pred_k_plus_k_qrs(:,i) = A_new_qrs*X0_qrs+J_qrs*(y_new_k_1_qrs-H_new_qrs*x_k_1_k_2_qrs);
            
            %         inov_vec = [inov inov_vec(1:end-1)];
            %         denum_vec = [H_new*(R1'*R1)*H_new'+R_new denum_vec(1:end-1)];
            
            
            % P0 =  (eye(size(Q))-K*H_new)*P_pred_k_k_1{i};
            % P_updat{1,i} = (P0);
            
            R2_qrs = qr([R_new_s_qrs' zeros(size(y,1),length_x_qrs);R1_qrs*H_qrs' R1_qrs]);
            P_s0_qrs =  (R2_qrs(end-length_x_qrs+1:end,end-length_x_qrs+1:end))';
            Ps_updat_qrs{1,i} = (P_s0_qrs);
            
            P_updat_qrs{1,i} = Ps_updat_qrs{1,i}*Ps_updat_qrs{1,i}';
            
            %         if i>300 && iter==25
            %             R_new = lambda*R_new+mean((1-lambda)*(inov_vec.^2)./denum_vec);
            %             R_new_s = R_new^.05;
            %         end
            
            % P_pred_k_plus_k{i} = A_new*P0*A_new'+Q_new;
            [~,R3_qrs] = (qr([P_s0_qrs'*A_new_qrs';Q_new_s_qrs]));
            R3_qrs = R3_qrs(1:length_x_qrs,1:length_x_qrs);
            Ps_pred_k_plus_k_qrs{i} = R3_qrs';
            
            P_pred_k_plus_k_qrs{i} = Ps_pred_k_plus_k_qrs{i}*Ps_pred_k_plus_k_qrs{i}';
            
        end
        
    end
    
    %     log_previous = log_liklihood;
    %     log_liklihood = log_liklihood_sum/length_y;
    %     log_new = log_liklihood;
    
    
    P_smoothed_pt{length_y_pt} = Ps_updat_pt{end};
    Sum_S_ti1_T_pt = eps;
    Sum_S_t_t_plus_T_pt = eps;
    
    % square_root_smoothing
    x_k_plus_N_pt = x_updat_pt(:,end); % x_k+1|N
    Ps_k_plus_N_pt = Ps_updat_pt{1,end};
    
    P_smoothed_qrs{length_y_qrs} = Ps_updat_qrs{end};
    Sum_S_ti1_T_qrs = eps;
    Sum_S_t_t_plus_T_qrs = eps;
    
    % square_root_smoothing
    x_k_plus_N_qrs = x_updat_qrs(:,end); % x_k+1|N
    Ps_k_plus_N_qrs = Ps_updat_qrs{1,end};
    
    for i=length_y_pt:-1:1
        K_pt = (Ps_updat_pt{1,i}*Ps_updat_pt{1,i}')*A_new_pt'/(Ps_pred_k_plus_k_pt{i}*Ps_pred_k_plus_k_pt{i}');
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                break;
            end
        K_smoothed_pt{i} = K_pt;
        x_smoothed_pt(:,i) = x_updat_pt(:,i)+K_pt*(x_k_plus_N_pt-x_pred_k_plus_k_pt(:,i));
        
        x_k_plus_N_pt = x_smoothed_pt(:,i);
        
        
        [~,R4_pt] = qr([Ps_updat_pt{1,i}'*A_new_pt'  Ps_updat_pt{1,i}';Q_new_s_pt zeros(length_x_pt,length_x_pt);zeros(length_x_pt,length_x_pt) Ps_k_plus_N_pt'*K_pt']);
        
        R4_pt = (R4_pt(length_x_pt+1:2*length_x_pt,end-length_x_pt+1:end));
        Ps_smoothed_pt{i} = R4_pt';
        Ps_k_plus_N_pt = Ps_smoothed_pt{i};
        
        P_smoothed_pt{i} = Ps_smoothed_pt{i}*Ps_smoothed_pt{i}' ;
        
        if i<=length_y_qrs
            K_qrs = (Ps_updat_qrs{1,i}*Ps_updat_qrs{1,i}')*A_new_qrs'/(Ps_pred_k_plus_k_qrs{i}*Ps_pred_k_plus_k_qrs{i}');
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                break;
            end
            
            K_smoothed_qrs{i} = K_qrs;
            x_smoothed_qrs(:,i) = x_updat_qrs(:,i)+K_qrs*(x_k_plus_N_qrs-x_pred_k_plus_k_qrs(:,i));
            
            x_k_plus_N_qrs = x_smoothed_qrs(:,i);
            
            
            [~,R4_qrs] = qr([Ps_updat_qrs{1,i}'*A_new_qrs'  Ps_updat_qrs{1,i}';Q_new_s_qrs zeros(length_x_qrs,length_x_qrs);zeros(length_x_qrs,length_x_qrs) Ps_k_plus_N_qrs'*K_qrs']);
            
            R4_qrs = (R4_qrs(length_x_qrs+1:2*length_x_qrs,end-length_x_qrs+1:end));
            Ps_smoothed_qrs{i} = R4_qrs';
            Ps_k_plus_N_qrs = Ps_smoothed_qrs{i};
            
            P_smoothed_qrs{i} = Ps_smoothed_qrs{i}*Ps_smoothed_qrs{i}' ;
            
        end
    end
    P_t_ti1_T_pt =cell(1,length_y_pt);
    P_t_ti1_T_pt{length_y_pt} = (eye(length_x_pt)-K_update_pt{length_y_pt}*H_pt)*A_new_pt*P_smoothed_pt{length_y_pt};
    
    Sum_S_t_ti1_T_pt =  eps*eye(length_x_pt);%P_t_ti1_T{length_y} ;
    Sum_S_ti1_T_pt  =  0*eye(length_x_pt);
    Sum_S_t_T_pt = 0*eye(length_x_pt);
    Sum_R_pt = 0;
    Sum_P_k_N_pt = 0; %Sum(P_k|N)
    Sum_y_x_k_N_pt = 0; %Sum(y*(x_k|N)')
    
    P_t_ti1_T_qrs =cell(1,length_y_qrs);
    P_t_ti1_T_qrs{length_y_qrs} = (eye(length_x_qrs)-K_update_qrs{length_y_qrs}*H_qrs)*A_new_qrs*P_smoothed_qrs{length_y_qrs};
    
    Sum_S_t_ti1_T_qrs =  eps*eye(length_x_qrs);%P_t_ti1_T{length_y} ;
    Sum_S_ti1_T_qrs  =  0*eye(length_x_qrs);
    Sum_S_t_T_qrs = 0*eye(length_x_qrs);
    Sum_R_qrs = 0;
    Sum_P_k_N_qrs = 0; %Sum(P_k|N)
    Sum_y_x_k_N_qrs = 0; %Sum(y*(x_k|N)')
    
    
    
    for i=length_y_pt:-1:2
        [warnMsg, warnId] = lastwarn;
        if ~isempty(warnMsg)
            break;
        end
        S_t_T_pt = P_smoothed_pt{i}+x_smoothed_pt(:,i)*x_smoothed_pt(:,i)';%S(t|T)
        S_ti1_T_pt = P_smoothed_pt{i-1}+x_smoothed_pt(:,i-1)*x_smoothed_pt(:,i-1)';%S(t-1|T)
        Sum_S_ti1_T_pt = Sum_S_ti1_T_pt + S_ti1_T_pt;
        Sum_S_t_T_pt = Sum_S_t_T_pt + S_t_T_pt;
        P_t_ti1_T_pt = K_smoothed_pt{i}*P_smoothed_pt{i}; %P(t,t-1|T)
        
        S_t_ti1_T_pt =  P_t_ti1_T_pt+ x_smoothed_pt(:,i)*x_smoothed_pt(:,i-1)';
        Sum_S_t_ti1_T_pt = Sum_S_t_ti1_T_pt + S_t_ti1_T_pt;
        Sum_R_pt = Sum_R_pt + y_new_pt(:,i)* y_new_pt(:,i)'-2*H_pt*x_smoothed_pt(:,i)*y(:,i) + H_pt*S_t_T_pt*H_pt';% (y_new(:,i)-H*x_smoothed(:,i))*(y_new(:,i)-H*x_smoothed(:,i))'-H*P_smoothed{i}*H'
        
        Sum_P_k_N_pt = Sum_P_k_N_pt + P_smoothed_pt{i}; %Sum(P_k|N)
        Sum_y_x_k_N_pt = Sum_y_x_k_N_pt + y_new_pt(:,i)*x_smoothed_pt(:,i)'; %Sum(y*(x_k|N)')
        

        if i<=length_y_qrs
            
            S_t_T_qrs = P_smoothed_qrs{i}+x_smoothed_qrs(:,i)*x_smoothed_qrs(:,i)';%S(t|T)
            S_ti1_T_qrs = P_smoothed_qrs{i-1}+x_smoothed_qrs(:,i-1)*x_smoothed_qrs(:,i-1)';%S(t-1|T)
            Sum_S_ti1_T_qrs = Sum_S_ti1_T_qrs + S_ti1_T_qrs;
            Sum_S_t_T_qrs = Sum_S_t_T_qrs + S_t_T_qrs;
            P_t_ti1_T_qrs = K_smoothed_qrs{i}*P_smoothed_qrs{i}; %P(t,t-1|T)
            
            S_t_ti1_T_qrs =  P_t_ti1_T_qrs+ x_smoothed_qrs(:,i)*x_smoothed_qrs(:,i-1)';
            Sum_S_t_ti1_T_qrs = Sum_S_t_ti1_T_qrs + S_t_ti1_T_qrs;
            Sum_R_qrs = Sum_R_qrs + y_new_qrs(:,i)* y_new_qrs(:,i)'-2*H_qrs*x_smoothed_qrs(:,i)*y(:,i) + H_qrs*S_t_T_qrs*H_qrs';% (y_new(:,i)-H*x_smoothed(:,i))*(y_new(:,i)-H*x_smoothed(:,i))'-H*P_smoothed{i}*H'
            
            Sum_P_k_N_qrs = Sum_P_k_N_qrs + P_smoothed_qrs{i}; %Sum(P_k|N)
            Sum_y_x_k_N_qrs = Sum_y_x_k_N_qrs + y_new_qrs(:,i)*x_smoothed_qrs(:,i)'; %Sum(y*(x_k|N)')
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                break;
            end
            
        end
        
    end
    [warnMsg, warnId] = lastwarn;
    if ~isempty(warnMsg)
        break;
    end
    A_new_EM_pt= Sum_S_t_ti1_T_pt/Sum_S_ti1_T_pt;
    Q_new_EM_pt = (Sum_S_t_T_pt-A_new_pt*Sum_S_ti1_T_pt)/(length_y_pt-1);
    % Q_new = (Sum_S_t_T-A*Sum_S_ti1_T)/(length(x_noisy)-1);
    R_new_EM_pt = Sum_R_pt/(length_y_pt-1);
    H_new_EM_pt = Sum_y_x_k_N_pt/Sum_P_k_N_pt;
    
    A_pt = triu(A_new_EM_pt);%A_new_EM;
    Q_pt = abs(Q_new_EM_pt.*eye(length_x_pt));%
    R_pt = abs(R_new_EM_pt.*eye(size(y,1)));
    Q_s_pt = Q_pt.^0.5; % square_root of Q (Q^1/2)
    R_s_pt = R_pt.^0.5; % square_root of R (R^1/2)
    

    
    
    X0_pt = x_smoothed_pt(:,1);
    P_s0_pt = Ps_smoothed_pt{1};
    
    A_new_pt = A_pt;
    H_new_pt = H_pt*A_new_pt-Psi_pt*H_pt;
    R_new_pt = R_pt;
    Q_new_s_pt = Q_s_pt;
    R_new_s_pt = R_s_pt;
    
    
    A_new_EM_qrs= Sum_S_t_ti1_T_qrs/Sum_S_ti1_T_qrs;
    Q_new_EM_qrs = (Sum_S_t_T_qrs-A_new_qrs*Sum_S_ti1_T_qrs)/(length_y_qrs-1);
    % Q_new = (Sum_S_t_T-A*Sum_S_ti1_T)/(length(x_noisy)-1);
    R_new_EM_qrs = Sum_R_qrs/(length_y_qrs-1);
    H_new_EM_qrs = Sum_y_x_k_N_qrs/Sum_P_k_N_qrs;
    
    A_qrs = triu(A_new_EM_qrs);%A_new_EM;
    Q_qrs = abs(Q_new_EM_qrs.*eye(length_x_qrs));%
    R_qrs = abs(R_new_EM_qrs.*eye(size(y,1)));
    Q_s_qrs = Q_qrs.^0.5; % square_root of Q (Q^1/2)
    R_s_qrs = R_qrs.^0.5; % square_root of R (R^1/2)
    
    X0_qrs = x_smoothed_qrs(:,1);
    P_s0_qrs = Ps_smoothed_qrs{1};
    
    A_new_qrs = A_qrs;
    H_new_qrs = H_qrs*A_new_qrs-Psi_qrs*H_qrs;
    R_new_qrs = R_qrs;
    Q_new_s_qrs = Q_s_qrs;
    R_new_s_qrs = R_s_qrs;
    
    x_smoothed_qrs_old = x_smoothed_qrs;
    x_smoothed_pt_old = x_smoothed_pt;
    


    
    
    
    
        %% Rewarping signals to original signals
        x_updat_total =  zeros(size(x));
        x_updat_total(1,ind_qrs) = x_updat_qrs(1,:);
        x_updat_total(1,ind_P_T) = x_updat_pt(1,:);
        figure(2),
        subplot(2,1,1)
        plot(1:length(x),x_updat_total(1,:),'r')
                legend({'KF EM'})
                        axis('tight')
        title(['EM Iteration number = ' num2str(iteration) '    Max number of Iterations = ' num2str(max_number_of_iterations)])

        subplot(2,1,2),plot(1:length(x),x_noisy,'b')
        legend({'Noisy'})
        axis('tight')

        pause(.01)
    
end
    [warnMsg, warnId] = lastwarn;
    if ~isempty(warnMsg)
     x_smoothed_qrs = x_smoothed_qrs_old;
     x_smoothed_pt = x_smoothed_pt_old;
     
    end

x_smoothed_P_T = zeros(size(x));
x_smoothed_P_T(1,ind_P_T) = x_smoothed_pt(1,:);


x_smoothed_total =  zeros(size(x));
x_smoothed_total(1,ind_qrs) = x_smoothed_qrs(1,:);
x_smoothed_total(1,ind_P_T) = x_smoothed_pt(1,:);


x_updat_P_T = zeros(size(x));
x_updat_P_T(1,ind_P_T) = x_updat_pt(1,:);


x_updat_total =  zeros(size(x));
x_updat_total(1,ind_qrs) = x_updat_qrs(1,:);
x_updat_total(1,ind_P_T) = x_updat_pt(1,:);




figure(3),
subplot(3,1,1)
        plot(1:length(x),x(1,:),'k')
                legend({'Original'})
                        axis('tight')
        title([file])
        subplot(3,1,2)
        plot(1:length(x),x_smoothed_total(1,:),'r')
                legend({'KS EM'})
                        axis('tight')

        subplot(3,1,3),plot(1:length(x),x_noisy,'b')
        legend({'Noisy'})
        axis('tight')



% % % % % % KF_SNR = 10*log10(mean((x(1,:)-y(1,:)).^2)/mean((x(1,:)-x_updat_total(1,:)).^2));
% % % % % % 
% % % % % % KS_SNR = 10*log10(mean((x(1,:)-y(1,:)).^2)/mean((x(1,:)-x_smoothed_total(1,:)).^2));
% % % % % % 
% % % % % % msepwrd_kf = MSEPWRD_error_v2( x,x_updat_total(1,:) );
% % % % % % msepwrd_ks = MSEPWRD_error_v2( x,x_smoothed_total(1,:) );
% % % % % % 
% % % % % % KF_PRD = Calculate_PRD(x(1,:),x_updat_total(1,:));
% % % % % % KS_PRD = Calculate_PRD(x(1,:),x_smoothed_total(1,:));
% % % % % % KF_corr = Calculate_Corr_coef(x(1,:),x_updat_total(1,:));
% % % % % % KS_corr = Calculate_Corr_coef(x(1,:),x_smoothed_total(1,:));
