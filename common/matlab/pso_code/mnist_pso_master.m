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
header_save_path = 'D:\Projects\mnist_pso\include';
results_save_path = 'D:\Projects\mnist_pso\pso_results';
% -----------------------------------------------------------------------
% Network Description
% -----------------------------------------------------------------------
mnist_description = struct();
mnist_description.net_structure = {'fc'; ...
                                   'act'; 'bn'; 'fc'; ...
                                   'act'; 'bn'; 'fc'; ...
                                   'mp'; 'act'; 'bn'; 'con'; ...
                                   'mp'; 'act'; 'bn'; 'con'; ...
                                   'input'};
mnist_description.activations = {'dlib::relu<'; 'dlib::prelu<'; 'dlib::sig<'; 'dlib::htan<'; 'dlib::elu<'; 'dlib::srelu<'};
mnist_description.con_num = 2;
mnist_description.fc_num = 3;
mnist_description.mp = {'con2d<'; 'con2d<'};
mnist_description.bn_num = nnz(strcmp(mnist_description.net_structure,'bn'));
mnist_description.act_num = nnz(strcmp(mnist_description.net_structure,'act'));
mnist_description.con_map = [1,2];
mnist_description.fc_map = [1,2,3];

net_description = mnist_description;

%% PSO member setup
g = 0.3;

c1 = 2.4;
c2 = 2.1;

phi = c1 + c2;
kap = 2/(abs(2 - phi - sqrt(phi^2 - 4*phi)));

N = 20;                         % population size
itr_max = 40;                   % number of iterations

num_classes = 10;

% -----------------------------------------------------------------------
% PSO Population Structure
% -----------------------------------------------------------------------
pso_member = struct();
pso_member.iteration = 0;
pso_member.number = 0;
pso_member.con = 0;
pso_member.fc = 0;
pso_member.bn = 0;
pso_member.act = 0;
% pso_member.input = 0;
% pso_member.input_file = 0;
% pso_member.crop_size = 0;
% pso_member.cond = 0;
pso_member.con_map = mnist_description.con_map;
pso_member.fc_map = mnist_description.fc_map;

pso_member_limits = struct();
                   
pso_member_limits.con = repmat([1,256,0,4,0,4], [mnist_description.con_num, 1]);

pso_member_limits.fc = repmat([1,1000], [mnist_description.fc_num, 1]);
pso_member_limits.fc(1,:) = [10,10];

pso_member_limits.bn = repmat([0,1], [mnist_description.bn_num, 1]);

pso_member_limits.act = repmat([1,6], [mnist_description.act_num, 1]);
% pso_member_limits.cond = [1,1];                     
% pso_member_limits.input = [5,5];
% pso_member_limits.input_file = [1,1];
% pso_member_limits.crop_size = [3,16];

% numvars = size(pso_member_limits.con,1)*3 + size(pso_member_limits.bn,1) + size(pso_member_limits.act,1) + ...
%             size(pso_member_limits.cond,1) + size(pso_member_limits.input,1) + size(pso_member_limits.input_file,1) + size(pso_member_limits.crop_size,1);

% -----------------------------------------------------------------------
% PSO Velocity Structure
% -----------------------------------------------------------------------
pso_velocity = struct();
pso_velocity.con = 0;
pso_velocity.fc = 0;
pso_velocity.bn = 0;
pso_velocity.act = 0;
% pso_velocity.input = 0;
% pso_velocity.input_file = 0;
% pso_velocity.crop_size = 0;
% pso_velocity.cond = 0;


pso_velocity_limits = struct();
pso_velocity_limits.con = repmat([-8,8,-1,1,-1,1], [mnist_description.con_num, 1]);

pso_velocity_limits.fc = repmat([-5,5], [mnist_description.fc_num, 1]);
pso_velocity_limits.fc(1,:) = [0,0];
 
pso_velocity_limits.bn = repmat([-1,1], [mnist_description.bn_num, 1]);

pso_velocity_limits.act = repmat([-1,1], [mnist_description.bn_num, 1]);

% pso_velocity_limits.cond = [-1,1];
% pso_velocity_limits.input = [0,0];
% pso_velocity_limits.input_file = [0,0];
% pso_velocity_limits.crop_size = [-1,1];       



%% Ask the user for the needed files

% Ask for the current results file.  If no file is selected assume that we
% are starting at iteration 1 and nothing has been run yet
file_filter = {'*.txt','Text Files';'*.*','All Files' };
mat_file_filter = {'*.mat','Mat Files';'*.*','All Files' };

startpath = results_save_path;
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
    input(strcat('Ready to process PSO Iteration:', 32, num2str(itr,'%03d')));
    
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
    
    for idx=1:N
        X(idx,1).iteration = 1;
        X(idx,1).number = idx;
        X(idx,1).con = cat(2,randi([1,256],mnist_description.con_num,1),randi([0,4],mnist_description.con_num,2));
        X(idx,1).fc = randi([1,1000],mnist_description.fc_num,1);
        X(idx,1).fc(1,1) = num_classes;
        
        X(idx,1).bn = randi(pso_member_limits.bn(1,:),mnist_description.bn_num,1);
        X(idx,1).act = randi(pso_member_limits.act(1,:),mnist_description.act_num,1);
%         X(idx,1).cond = randi(pso_member_limits.cond,1,1);
%         X(idx,1).input = randi(pso_member_limits.input,1,1);
%         X(idx,1).input_file = randi(pso_member_limits.input_file,1,1);
%         X(idx,1).crop_size = randi(pso_member_limits.crop_size,1,1);

        X(idx,1) = validate_mnist_pso_member(X(idx,1), pso_member_limits);

        % -5 + (5+5)*rand(10,1)
        V(idx,1).con = cat(2,(-8 + (16)*rand(mnist_description.con_num,1)),(-1+(2)*rand(mnist_description.con_num,2)));
        V(idx,1).fc = (-10 + (20)*rand(mnist_description.fc_num,1));
        V(idx,1).bn = -1 + (1+1)*rand(mnist_description.bn_num,1);
        V(idx,1).act = -1 + (1+1)*rand(mnist_description.act_num,1);
%         V(idx,1).cond = -1 + (1+1)*rand(1,1);
%         V(idx,1).input = -1 + (1+1)*rand(1,1);
%         V(idx,1).input_file = -1 + (1+1)*rand(1,1);
%         V(idx,1).crop_size = -1 + (1+1)*randi(1,1);    
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
        header_name = strcat('mnist_net_pso_',num2str(idx,'%02d.h'));
        pso_header_filename = fullfile(header_save_path, header_name);
        fprintf('Building: %s\n', header_name);
        build_mnist_net(pso_header_filename, X(idx,1), mnist_description);
    end

    
%% save out the results of everything that has been generated to date

% save the P, X, V, G, p_best, g_best, iteration number in a mat file to be read on future
% iterations

    [mat_save_file, mat_save_path] = uiputfile(mat_file_filter, 'Select Save Mat File', fullfile(startpath,'MNIST_PSO_Test_Results.mat'));
    itr = 1;

    save(fullfile(mat_save_path,mat_save_file),'P','X','V','G','F','p_best','g_best','itr');

elseif(itr == 1)
    
%% Step 4:
% - Evaluate each of the networks
%   - Get the NRMSE and SSIM for each population member
%   - f(x(k)) = g*NRMSE + (1-g)*(1-SSIM)
    input(strcat('Ready to process PSO Iteration:', 32, num2str(itr,'%03d')));

    for idx=1:size(results,1)
        pop_num = str2double(results{idx}{2});
        
        acc_tr = str2double(results{idx}{5});
        acc_te = str2double(results{idx}{8});
        
        F(pop_num,itr) = 1.0 - acc_te;
        
        
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

    V = update_mnist_velocity(X, V, P, G, c1, c2, kap, N, itr, pso_velocity_limits);

    X = update_mnist_population(X, V, N, itr, pso_member_limits);

%% Step 6:
% - Generate the network header file
    fprintf('Building Headers...\n');
    for idx=1:N
        header_name = strcat('mnist_net_pso_',num2str(idx,'%02d.h'));
        pso_header_filename = fullfile(header_save_path, header_name);
        fprintf('Building: %s\n', header_name);
        build_mnist_net(pso_header_filename, X(idx,itr+1), mnist_description);
    end

    
%% save out the results of everything that has been generated to date

% save the P, X, V, G, p_best, g_best, iteration number in a mat file to be read on future
% iterations
    itr = itr + 1;

    save(fullfile(mat_save_path, mat_save_file), 'P', 'X', 'V', 'G', 'F', 'p_best', 'g_best', 'itr', 'net_description');
    
else
    
    input(strcat('Ready to process PSO Iteration:', 32, num2str(itr,'%03d')));
    
%% Step 7:    
% - Evaluate the new networks (f(x(k+1))
%   - Get the NRMSE and SSIM for each population member
%   - f(x(k+1)) = g*NRMSE + (1-g)*(1-SSIM)
    for idx=1:size(results,1)
        pop_num = str2double(results{idx}{2});
        
        acc_tr = str2double(results{idx}{5});
        acc_te = str2double(results{idx}{8});
        
        F(pop_num,itr) = 1.0 - acc_te;
        
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
            
        elseif(F(idx,itr) == F(idx,itr-1))
            x_sum_con = g*sum(X(idx, itr).con(:,1));
            x_sum_fc = (1-g)*sum(X(idx, itr).fc(:,1));
            p_sum_con = g*sum(P(idx, itr-1).con(:,1));
            p_sum_fc = (1-g)*sum(P(idx, itr-1).fc(:,1));
            
            if((x_sum_con + x_sum_fc) < (p_sum_con + p_sum_fc))
                P(idx,itr) = X(idx,itr);
            else
                P(idx,itr) = P(idx,itr-1);
                P(idx,itr).iteration = itr;
            end            
        else
            P(idx,itr) = P(idx,itr-1);
            P(idx,itr).iteration = itr;
        end
    end
    
    if(p_best(1, itr) < g_best(itr-1))
        G(itr) = X(f_min_idx, itr);
        g_best(itr) = p_best(1, itr);
        
    % deal with the tie breaker
    elseif(p_best(1, itr) == g_best(itr-1))
        x_sum_con = g*sum(X(f_min_idx, itr).con(:,1));
        x_sum_fc = (1-g)*sum(X(f_min_idx, itr).fc(:,1));
        g_sum_con = g*sum(G(itr-1).con(:,1));
        g_sum_fc = (1-g)*sum(G(itr-1).fc(:,1));
        
        if((x_sum_con+x_sum_fc) < (g_sum_con + g_sum_fc))
            G(itr) = X(f_min_idx, itr);
            g_best(itr) = p_best(1, itr);
        else
            G(itr) = G(itr-1);
            g_best(itr) = g_best(itr-1);                        
        end
        
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

    V = update_mnist_velocity(X, V, P, G, c1, c2, kap, N, itr, pso_velocity_limits);

    X = update_mnist_population(X, V, N, itr, pso_member_limits);

%% Step 10:
% - Generate the network header file

    for idx=1:N
        header_name = strcat('mnist_net_pso_',num2str(idx,'%02d.h'));
        pso_header_filename = fullfile(header_save_path, header_name);
        fprintf('Building: %s\n', header_name);
        build_mnist_net(pso_header_filename, X(idx,itr+1), mnist_description);
    end

    
%% save out the results of everything that has been generated to date

% save the P, X, V, G, p_best, g_best, iteration number in a mat file to be read on future
% iterations
    itr = itr + 1;

    save(fullfile(mat_save_path, mat_save_file), 'P', 'X', 'V', 'G', 'F', 'p_best', 'g_best', 'itr', 'net_description');

end

%% end

if(itr>1)
    fprintf('\ng_best: %2.4f\n', g_best(itr-1));
    fprintf('G Best:\n');
    disp(G(itr-1))

    fprintf('\np_best: %2.4f\n', p_best(1,itr-1));
    fprintf('P Best:\n');
    disp(P(f_min_idx,itr-1));
end

fprintf('\nComplete!\n');

return;

%% save a specific header file

idx_itr = 31;
idx_N = 13;

header_name = strcat('mnist_net_pso_',num2str(idx_N,'%02d.h'));
pso_header_filename = fullfile(header_save_path, header_name);
fprintf('Building: %s\n', header_name);
build_mnist_net(pso_header_filename, X(idx_N,idx_itr), mnist_description);


