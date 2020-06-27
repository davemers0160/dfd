format long g
format compact
clc
close all
clearvars

% get the location of the script file to save figures
full_path = mfilename('fullpath');
[startpath,  filename, ext] = fileparts(full_path);
plot_num = 1;

%% first time run the program should initiate the initial gerenation of the population
% and then create the headers files.  Next it should save the population
% and velocity numbers for each population member

% after the creation of the header files the networks need to be trainined
% and evaluated.

% after the evaluation the numbers need to be read back into the program.
% The new velocity needs to be computed and checked.  Next the population 
% members need to be generated and checked.  Then the new header files need
% to be generated.


%% these are the constants that will be used
header_save_path = 'D:\IUPUI\DfD\common\pso';
results_save_path = 'D:\IUPUI\DfD\dfd_dnn_pso\results\';
% -----------------------------------------------------------------------
% Network Description
% -----------------------------------------------------------------------
% net_description_v6 = struct();
% net_description_v6.net_structure = {'con'; 'act'; 'bn';  'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'act'; 'bn'; ...
%                 'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'act'; 'bn'; 'con'; 'cond'; ...
%                 'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'act'; 'bn'; 'con'; 'cond'; ...
%                 'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'act'; 'bn'; 'con'; 'cond'; ...
%               	'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'act'; 'bn'; ...
%                 'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'act'; 'bn'; 'con'; ...
%                 'input'};
% net_description_v6.activations = {'dlib::relu<'; 'dlib::prelu<'; 'dlib::sig<'; 'dlib::htan<'; 'dlib::elu<'; 'dlib::srelu<'};
% net_description_v6.tags = {'DTI_0<'; 'DTO_0<'};
% net_description_v6.concats = {'dlib::concat2<DTO_0, DTI_0,'; 'dlib::concat2<DTO_1, DTI_1,'};
% net_description_v6.cond = {'con21d21<'; 'con32d32<'; 'con33d33<'};
% net_description_v6.con_num = nnz(strcmp(net_description_v6.net_structure,'con')) + nnz(strcmp(net_description_v6.net_structure,'cond'));
% net_description_v6.bn_num = nnz(strcmp(net_description_v6.net_structure,'bn'));
% net_description_v6.act_num = nnz(strcmp(net_description_v6.net_structure,'act'));

net_description_v14 = struct();
net_description_v14.net_structure = {'con'; 'act'; 'bn';  'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'act'; 'bn'; 'con'; ...
                'concat'; 'tags'; 'cont2u'; 'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'act'; 'bn'; 'con'; ...               
                'concat'; 'tags'; 'cont2u'; 'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'cond'; ...
                'tags'; 'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'cond'; ...
                'tags'; 'add_prev1'; 'con'; 'act'; 'bn'; 'con'; 'tag1'; 'act'; 'bn'; 'con'; 'input'};
net_description_v14.activations = {'dlib::relu<'; 'dlib::prelu<'; 'dlib::sig<'; 'dlib::htan<'; 'dlib::elu<'; 'dlib::srelu<'};
net_description_v14.tags = {'DTI_1<'; 'DTI_2<'; 'DTO_2<'; 'DTO_1<'};            
net_description_v14.concats = {'dlib::concat2<DTO_1, DTI_1,'; 'dlib::concat2<DTO_2, DTI_2,'};            
net_description_v14.cond = {'con2d<'; 'con2d<'};
%net_description_v14.con_num = nnz(strcmp(net_description_v14.net_structure,'con')) + nnz(strcmp(net_description_v14.net_structure,'cond')) + nnz(strcmp(net_description_v14.net_structure,'cont2u'));
net_description_v14.con_num = 13;
net_description_v14.bn_num = nnz(strcmp(net_description_v14.net_structure,'bn'));
net_description_v14.act_num = nnz(strcmp(net_description_v14.net_structure,'act'));
net_description_v14.con_map = [1,2,3,2,4,5,6,5,7,8,9,8,10,11,10,12,13,12];


% copy what net is used into the generic container
net_description = net_description_v14;



%% PSO member setup
%g = [0.4, 10, 10, 0.05];        % gamma value for combining NRMSE and SSIM metrics
g = [0.3, 10, 10, 0.05];        % gamma value for combining NRMSE and SSIM metrics
vg = [0.7, 0.4];

c1 = 2.4;
c2 = 2.1;

phi = c1 + c2;
kap = 2/(abs(2 - phi - sqrt(phi^2 - 4*phi)));

N = 20;                         % population size
itr_max = 20;                   % number of iterations

num_classes = 256;

% -----------------------------------------------------------------------
% PSO Population Structure
% -----------------------------------------------------------------------
pso_member = struct();
pso_member.iteration = 0;
pso_member.number = 0;
pso_member.con = 0;
pso_member.bn = 0;
pso_member.act = 0;
pso_member.input = 0;
pso_member.input_file = 0;
pso_member.crop_size = 0;
pso_member.cond = 0;
pso_member.con_map = net_description_v14.con_map;

pso_member_limits = struct();
% pso_member_limits.con = [[256,256,3,3,1,1];[8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];...
%                          [8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];...
%                          [8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];...
%                          [8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4];[8,512,0,4,0,4]];
                     
pso_member_limits.con = repmat([8,512,0,4,0,4], [net_description_v14.con_num, 1]);
pso_member_limits.con(1,:) = [256,256,1,1,1,1];

% pso_member_limits.bn = [[0,1];[0,1];[0,1];[0,1];[0,1];[0,1];[0,1];...
%                         [0,1];[0,1];[0,1];[0,1];[0,1];[0,1]];
pso_member_limits.bn = repmat([0,1], [net_description_v14.bn_num, 1]);

% pso_member_limits.act = [[1,6];[1,6];[1,6];[1,6];[1,6];[1,6];[1,6];...
%                          [1,6];[1,6];[1,6];[1,6];[1,6];[1,6]];
pso_member_limits.act = repmat([1,6], [net_description_v14.act_num, 1]);
pso_member_limits.cond = [1,1];                     
pso_member_limits.input = [5,5];
pso_member_limits.input_file = [1,1];
pso_member_limits.crop_size = [3,16];

% numvars = size(pso_member_limits.con,1)*3 + size(pso_member_limits.bn,1) + size(pso_member_limits.act,1) + ...
%             size(pso_member_limits.cond,1) + size(pso_member_limits.input,1) + size(pso_member_limits.input_file,1) + size(pso_member_limits.crop_size,1);

% -----------------------------------------------------------------------
% PSO Velocity Structure
% -----------------------------------------------------------------------
pso_velocity = struct();
pso_velocity.con = 0;
pso_velocity.bn = 0;
pso_velocity.act = 0;
pso_velocity.input = 0;
pso_velocity.input_file = 0;
pso_velocity.crop_size = 0;
pso_velocity.cond = 0;


pso_velocity_limits = struct();
% pso_velocity_limits.con = [[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];...
%                            [-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];...
%                            [-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];...
%                            [-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1];[-8,8,-1,1,-1,1]];
pso_velocity_limits.con = repmat([-16,16,-1,1,-1,1], [net_description_v14.con_num, 1]);
pso_velocity_limits.con(1,:) = [0,0,0,0,0,0];
% pso_velocity_limits.bn = [[-1,1];[-1,1];[-1,1];[-1,1];[-1,1];[-1,1];[-1,1];...
%                           [-1,1];[-1,1];[-1,1];[-1,1];[-1,1];[-1,1]];  
pso_velocity_limits.bn = repmat([-1,1], [net_description_v14.bn_num, 1]);

% pso_velocity_limits.act = [[-1,1];[-1,1];[-1,1];[-1,1];[-1,1];[-1,1];[-1,1];...
%                            [-1,1];[-1,1];[-1,1];[-1,1];[-1,1];[-1,1]];
pso_velocity_limits.act = repmat([-1,1], [net_description_v14.bn_num, 1]);

pso_velocity_limits.cond = [-1,1];
pso_velocity_limits.input = [0,0];
pso_velocity_limits.input_file = [0,0];
pso_velocity_limits.crop_size = [-1,1];       



%% Ask the user for the needed files

% Ask for the current results file.  If no file is selected assume that we
% are starting at iteration 1 and nothing has been run yet
file_filter = {'*.txt','Text Files';'*.*','All Files' };
mat_file_filter = {'*.mat','Mat Files';'*.*','All Files' };

startpath = 'D:\IUPUI\DfD\dfd_dnn_pso\pso_results\';
[results_file, results_file_path] = uigetfile(file_filter, 'Select Results File', startpath);
if(results_file_path == 0)
    itr = 0;
else
    results = parse_input_parameters(fullfile(results_file_path, results_file));
    %itr = str2double(results{1}{1});
    [mat_save_file, mat_save_path] = uigetfile(mat_file_filter, 'Select Mat File', results_file_path);
    if(mat_save_path == 0)
        return;
    end
    load(fullfile(mat_save_path,mat_save_file));
end

commandwindow;

%% this is the pseudo code for the PSO

% Step 1: Initialize the population and the initial velocity vectors or
% load in values from a previous iteration
% first iteration
%itr = 1;
if (itr == 0)
    
    % F represents the results of evaluating the trainined architecture
    F = zeros(N, itr_max+1);
    
    % X represents the population to be tested for a given iteration
    X = repmat(pso_member,N,itr_max+1);
    
    % V represents the velocity vectors for a given iteration
    V = repmat(pso_velocity,N,itr_max+1);

    % P represents the best poplation member for that iteration
    P = repmat(pso_member,N,itr_max+1);
    
    % G represents the best population member for all iterations
    G = repmat(pso_member,1,itr_max+1);
    
    % p_best contains the info about the best results for each iteration
    p_best = zeros(3, itr_max+1);
    
    % g_best contains the info about the best results overall
    g_best = zeros(1, itr_max+1);
    
    % fill in the first one with the best to date
    idx = 1;
    X(idx,1).iteration = 1;
    X(idx,1).number = idx;
    X(idx,1).con = cat(2,[256;128;128;128;256;256;256;512;512;256;256;128;128], ...
                         repmat([1,1],[net_description_v14.con_num,1]));
    X(idx,1).con(1,1) = num_classes;
    X(idx,1).bn = ones(net_description_v14.bn_num,1);
    X(idx,1).act = 2*ones(net_description_v14.act_num,1);
    X(idx,1).cond = 1;
    X(idx,1).input = 5;
    X(idx,1).input_file = 1;
    X(idx,1).crop_size = 8;

    X(idx,1) = validate_pso_member(X(idx,1), pso_member_limits); 
    
    V(idx,1).con = cat(2,(-8 + (8+8)*rand(net_description_v14.con_num,1)),(-1+(1+1)*rand(net_description_v14.con_num,2)));
    V(idx,1).bn = -1 + (1+1)*rand(net_description_v14.bn_num,1);
    V(idx,1).act = -1 + (1+1)*rand(net_description_v14.act_num,1);
    V(idx,1).cond = -1 + (1+1)*rand(1,1);
    V(idx,1).input = -1 + (1+1)*rand(1,1);
    V(idx,1).input_file = -1 + (1+1)*rand(1,1);
    V(idx,1).crop_size = -1 + (1+1)*randi(1,1);
    
    for idx=2:N
        X(idx,1).iteration = 1;
        X(idx,1).number = idx;
        X(idx,1).con = cat(2,randi([8,512],net_description_v14.con_num,1),randi([0,4],net_description_v14.con_num,2));
        X(idx,1).con(1,1) = num_classes;
        X(idx,1).bn = randi(pso_member_limits.bn(1,:),net_description_v14.bn_num,1);
        X(idx,1).act = randi(pso_member_limits.act(1,:),net_description_v14.act_num,1);
        X(idx,1).cond = randi(pso_member_limits.cond,1,1);
        X(idx,1).input = randi(pso_member_limits.input,1,1);
        X(idx,1).input_file = randi(pso_member_limits.input_file,1,1);
        X(idx,1).crop_size = randi(pso_member_limits.crop_size,1,1);

        X(idx,1) = validate_pso_member(X(idx,1), pso_member_limits);

        % -5 + (5+5)*rand(10,1)
        V(idx,1).con = cat(2,(-8 + (8+8)*rand(net_description_v14.con_num,1)),(-1+(1+1)*rand(net_description_v14.con_num,2)));
        V(idx,1).bn = -1 + (1+1)*rand(net_description_v14.bn_num,1);
        V(idx,1).act = -1 + (1+1)*rand(net_description_v14.act_num,1);
        V(idx,1).cond = -1 + (1+1)*rand(1,1);
        V(idx,1).input = -1 + (1+1)*rand(1,1);
        V(idx,1).input_file = -1 + (1+1)*rand(1,1);
        V(idx,1).crop_size = -1 + (1+1)*randi(1,1);    
    end

%% Step 2: 
% - Set the population best for each member to the initial population values
%   - P(0) = x(0)
%    P(:,1) = X(:,1);

%% Step 3:
% - Generate the header file that contains the following info:
%   - input type
%   - crop size
%   - net parameters
%   - maybe look at using the input file as a parameter
    fprintf('Building Headers...\n');
    for idx=1:N
        header_name = strcat('dfd_net_v14_pso_',num2str(idx,'%02d.h'));
        pso_header_filename = fullfile(header_save_path, header_name);
        fprintf('Building: %s\n', header_name);
        build_dfd_net_rw(pso_header_filename, X(idx,1), net_description_v14);
    end

    
%% save out the results of everything that has been generated to date

% save the P, X, V, G, p_best, g_best, iteration number in a mat file to be read on future
% iterations

    [mat_save_file, mat_save_path] = uiputfile(mat_file_filter, 'Select Save Mat File', fullfile(startpath,'PSO_Test_Results.mat'));
    itr = 1;

    save(fullfile(mat_save_path, mat_save_file), 'P', 'X', 'V', 'G', 'F', 'p_best', 'g_best', 'itr', 'g', 'vg', 'net_description');

elseif(itr == 1)
    
%% Step 4:
% - Evaluate each of the networks
%   - Get the NRMSE and SSIM for each population member
%   - f(x(k)) = g*NRMSE + (1-g)*(1-SSIM)
    fprintf('Iteration #: %d\n', itr);

    for idx=1:size(results,1)
        pop_num = str2double(results{idx}{2});
        
        nmae_tr = str2double(results{idx}{3});
        nrmse_tr = str2double(results{idx}{4});
        ssim_tr = str2double(results{idx}{5});
        var_gt_tr = str2double(results{idx}{6});
        var_dm_tr = str2double(results{idx}{7});
        
        nmae_te = str2double(results{idx}{8});
        nrmse_te = str2double(results{idx}{9});
        ssim_te = str2double(results{idx}{10});
        var_gt_te = str2double(results{idx}{11});
        var_dm_te = str2double(results{idx}{12});
               
        tr = (nmae_tr + nrmse_tr + (1-ssim_tr));
        te = (nmae_te + nrmse_te + (1-ssim_te));
        
        p_nrmse = g(2) * abs(min(nrmse_te-nrmse_tr,0));
        p_ssim = g(3) * abs(min(ssim_tr-ssim_te,0));
        
        p_var_tr = g(4)*max(abs(var_gt_tr - var_dm_tr) - (vg(1)*var_gt_tr), 0);
        p_var_te = g(4)*max(abs(var_gt_te - var_dm_te) - (vg(2)*var_gt_te), 0);
        
        %F(pop_num,itr) =  g(1)*(nrmse_te+nrmse_tr) + (1-g(1))*(2-ssim_te-ssim_tr) + p_nrmse + p_ssim; %+ p_var_tr %+ p_var_te;
        F(pop_num,itr) = g(1)*tr + (1-g(1))*te;
    end
    
    % get the minimum and min index - F(f_min_idx,1) tells us which
    % population member performed the best
    [p_best(1, itr), f_min_idx] = min(F(:,itr));
    p_best(2, itr) = mean(F(:,itr));
    [p_best(3, itr), f_max_idx] = max(F(:,itr));
    g_best(itr) = p_best(1, itr);
    
    % Update P and G
    P(:,itr) = X(:,itr);    
    G(itr) = X(f_min_idx, itr);


%% Step 5:
% - Generate r(k) and s(k) for each of the population members
% - Calculate the updated velocity vector and then check the velocity limits
%   - v(k+1) = kap*(v(k) + c1*r(k) o (P(k) - x(k)) + c2*s(k) o (G - x(k)))
% - Calculate the new population members and check the population limits
%   - x(k+1) = x(k) + v(k+1)

    V = update_velocity(X, V, P, G, c1, c2, kap, N, itr, pso_velocity_limits);

    X = update_population(X, V, N, itr, pso_member_limits);

%% Step 6:
% - Generate the network header file
    fprintf('Building Headers...\n');
    for idx=1:N
        header_name = strcat('dfd_net_v14_pso_',num2str(idx,'%02d.h'));
        pso_header_filename = fullfile(header_save_path, header_name);
        fprintf('Building: %s\n', header_name);
        build_dfd_net_rw(pso_header_filename, X(idx,itr+1), net_description_v14);
    end

    
%% save out the results of everything that has been generated to date

% save the P, X, V, G, p_best, g_best, iteration number in a mat file to be read on future
% iterations
    itr = itr + 1;

    save(fullfile(mat_save_path, mat_save_file), 'P', 'X', 'V', 'G', 'F', 'p_best', 'g_best', 'itr', 'g', 'vg', 'net_description');
   
else
    
    input(strcat('Ready to process PSO Iteration:', 32, num2str(itr,'%03d')));
    
%% Step 7:    
% - Evaluate the new networks (f(x(k+1))
%   - Get the NRMSE and SSIM for each population member
%   - f(x(k+1)) = g*NRMSE + (1-g)*(1-SSIM)
    for idx=1:size(results,1)
        pop_num = str2double(results{idx}{2});
        
        nmae_tr = str2double(results{idx}{3});
        nrmse_tr = str2double(results{idx}{4});
        ssim_tr = str2double(results{idx}{5});
        var_gt_tr = str2double(results{idx}{6});
        var_dm_tr = str2double(results{idx}{7});
        
        nmae_te = str2double(results{idx}{8});
        nrmse_te = str2double(results{idx}{9});
        ssim_te = str2double(results{idx}{10});
        var_gt_te = str2double(results{idx}{11});
        var_dm_te = str2double(results{idx}{12});
               
        tr = (nmae_tr + nrmse_tr + (1-ssim_tr));
        te = (nmae_te + nrmse_te + (1-ssim_te));
        
        p_nrmse = g(2) * abs(min(nrmse_te-nrmse_tr,0));
        p_ssim = g(3) * abs(min(ssim_tr-ssim_te,0));
        
        p_var_tr = g(4)*max(abs(var_gt_tr - var_dm_tr) - (vg(1)*var_gt_tr), 0);
        p_var_te = g(4)*max(abs(var_gt_te - var_dm_te) - (vg(2)*var_gt_te), 0);
        
        %F(pop_num,itr) =  g(1)*(nrmse_te+nrmse_tr) + (1-g(1))*(2-ssim_te-ssim_tr) + p_nrmse + p_ssim; %+ p_var_tr %+ p_var_te;
        F(pop_num,itr) = g(1)*tr + (1-g(1))*te;
    end
  
%% Step 8:
% - Update the population best for the next iteration
%   - if f(x(k+1)) < f(P(k)) -> P(k+1) = x(k+1) else P(k+1) = P(k)
%   - if f(x(k+1)) < f(G(k)) -> G(k+1) = x(k+1) else G(k+1) = G(k)    
    
    % get the minimum and min index - F(f_min_idx,1) tells us which
    % population member performed the best
    [p_best(1, itr), f_min_idx] = min(F(:,itr));
    p_best(2, itr) = mean(F(:,itr));
    [p_best(3, itr), f_max_idx] = max(F(:,itr));
    
    % Update P and G
    for idx=1:N
        if(F(idx,itr) < F(idx,itr-1))
            P(idx,itr) = X(idx,itr);
        else
            P(idx,itr) = P(idx,itr-1);
        end
    end
    
    if(p_best(1, itr) < g_best(itr-1))
        G(itr) = X(f_min_idx, itr);
        g_best(itr) = p_best(1, itr);
    else
        G(itr) = G(itr-1);
        g_best(itr) = g_best(itr-1);
    end 


%% Step 9:
% - Generate r(k) and s(k) for each of the population members
% - Calculate the updated velocity vector and then check the velocity limits
%   - v(k+1) = kap*(v(k) + c1*r(k) o (P(k) - x(k)) + c2*s(k) o (G - x(k)))
% - Calculate the new population members and check the population limits
%   - x(k+1) = x(k) + v(k+1)

    V = update_velocity(X, V, P, G, c1, c2, kap, N, itr, pso_velocity_limits);

    X = update_population(X, V, N, itr, pso_member_limits);

%% Step 10:
% - Generate the network header file

    for idx=1:N
        header_name = strcat('dfd_net_v14_pso_',num2str(idx,'%02d.h'));
        pso_header_filename = fullfile(header_save_path, header_name);
        fprintf('Building: %s\n', header_name);
        build_dfd_net_rw(pso_header_filename, X(idx,itr+1), net_description_v14);
    end

    
%% save out the results of everything that has been generated to date

% save the P, X, V, G, p_best, g_best, iteration number in a mat file to be read on future
% iterations
    itr = itr + 1;

    save(fullfile(mat_save_path, mat_save_file), 'P', 'X', 'V', 'G', 'F', 'p_best', 'g_best', 'itr', 'g', 'vg', 'net_description');

end

%% display the results

if(itr>1)
    fprintf('\ng_best: %2.4f\n', g_best(itr-1));
    fprintf('G Best:\n');
    disp(G(itr-1))

    fprintf('\np_best: %2.4f\n', p_best(1,itr-1));
    fprintf('P Best:\n');
    disp(P(f_min_idx,itr-1));
end

fprintf('\nComplete!\n');





