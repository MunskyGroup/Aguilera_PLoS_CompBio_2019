function [X_output,Ncol] =  SSA_longNames(k,t_array,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP)
%% Defining parameters
exclusion = 9;                      % exclusion volume
maximum_Number_Ribosomes = 1000;       % maximum number of ribosomes
k_elongation = k(2:end-1);           % elongation constants
k_bind = k(1);                       % Initiation rate
k_termination = k(end);              % Termination rate.
gene_Length = length(k_elongation);  % number of codons.
X_State = zeros(maximum_Number_Ribosomes,1);                       % Initial state (length of the longest protein = 1)
t = t_array(1);                      % time array
%% Preallocating memory
number_TimePoints = length(t_array);
t_final = t_array(number_TimePoints);
X_output = zeros(maximum_Number_Ribosomes,number_TimePoints);
% collisions =0;
Ncol = zeros(1,0);
col = zeros(1,maximum_Number_Ribosomes);

%% Run SSA
iteration = 1;
while t < t_final
    
    if t >= timePerturbationApplication && evaluatingInhibitor ==1
        Inhibitor_Condition = 0; % Set to 0 during its application
    else
        Inhibitor_Condition = 1; % Set to 1 before its application
    end
    if evaluatingFRAP ==1 && t >= timePerturbationApplication && t <= timePerturbationApplication+10
        X_State = zeros(maximum_Number_Ribosomes,1);
    end
    %% Compute propensity functions
    Nribosomes = sum(X_State > 0);
    X_State(Nribosomes + 1,1) = 0;
    temp_Stoichiometry = eye(Nribosomes + 1);
    temp_Propensities = zeros(Nribosomes + 1,1);
    temp_Propensities(X_State > 0,1) = k_elongation(X_State(X_State > 0));
    % termination
    if X_State(1) >= gene_Length -2
        temp_Stoichiometry(:,1) = [X_State(2:Nribosomes+1);0] - X_State(1:Nribosomes+1);
        temp_Propensities(1,1) = k_termination;
    end
    % initiation
    if Nribosomes == 0 || X_State(Nribosomes) > exclusion
        temp_Propensities(Nribosomes+1,1) = k_bind * Inhibitor_Condition;
    end
    % elongation
    elongation_Condition = (X_State(2:Nribosomes) + exclusion) < X_State(1:Nribosomes-1);
    temp_Propensities(2:Nribosomes,1) = temp_Propensities(2:Nribosomes,1) .* elongation_Condition;
    
    % deffining prop and stoich
    propensities = temp_Propensities;
    stoichiometry = temp_Stoichiometry;
    %% Update time
    sum_Propensities = sum(propensities);
    t = t - log(rand) / sum_Propensities;
    %% Generate output
    while iteration <= number_TimePoints && t > t_array(iteration)
        X_output(1:size(X_State,1),iteration) = X_State;
        iteration = iteration + 1;
    end
    
    %% Update state
    
    if t < t_final
        %% Select reaction
        r2 = sum_Propensities * rand;
        i = 1;
        tmp = propensities(i);
        while tmp < r2
            i = i + 1;
            tmp = tmp + propensities(i);
        end
        
        %% update collisions
        if sum(stoichiometry(:,i)) > 0
          %  j = find(stoichiometry(:,i));
            if i ~= 1
                                 j = find(stoichiometry(:,i));
                if X_State(j-1) == X_State(j)+exclusion+1
                    col(j) = col(j) + 1;
                end
            end
        else
            Ncol = [Ncol,col(1)];
            col = [col(2:end),0];
        end
        
        %% update state
        X_State(1:Nribosomes+1) = X_State(1:Nribosomes+1) + stoichiometry(:,i);
    end
    
end




% function [X_output,collisions] =  SSA_longNames(k,t_array,timePerturbationApplication,evaluatingInhibitor,evaluatingFRAP)
% %% Defining parameters
% exclusion = 9;                      % exclusion volume
% maximum_Number_Ribosomes = 200;       % maximum number of ribosomes
% k_elongation = k(2:end-1);           % elongation constants
% k_bind = k(1);                       % Initiation rate
% k_termination = k(end);              % Termination rate.
% gene_Length = length(k_elongation);  % number of codons.
% X_State = zeros(maximum_Number_Ribosomes,1);                       % Initial state (length of the longest protein = 1)
% t = t_array(1);                      % time array
% %% Preallocating memory
% number_TimePoints = length(t_array);
% t_final = t_array(number_TimePoints);
% X_output = zeros(maximum_Number_Ribosomes,number_TimePoints);
% counter_No_Ribosomes = 0;
% collisions =0;
% %% Run SSA
% iteration = 1;
% while t < t_final
%
%     if t >= timePerturbationApplication && evaluatingInhibitor ==1
%         Inhibitor_Condition = 0; % Set to 0 during its application
%     else
%         Inhibitor_Condition = 1; % Set to 1 before its application
%     end
%     if evaluatingFRAP ==1 && t >= timePerturbationApplication && t <= timePerturbationApplication+10
%         X_State = zeros(maximum_Number_Ribosomes,1);
%     end
%     %% Compute propensity functions
%     Nribosomes = sum(X_State > 0);
%     X_State(Nribosomes + 1,1) = 0;
%     %X_State = X_State(1:current_NO_Ribosomes + 1);
%     temp_Stoichiometry = eye(Nribosomes + 1);
%     temp_Propensities = zeros(Nribosomes + 1,1);
%     temp_Propensities(X_State > 0,1) = k_elongation(X_State(X_State > 0));
%     % termination
%     if X_State(1) >= gene_Length -2
%         temp_Stoichiometry(:,1) = [X_State(2:Nribosomes+1);0] - X_State(1:Nribosomes+1);
%         temp_Propensities(1,1) = k_termination;
%     end
%     % initiation
%     if Nribosomes == 0 || X_State(Nribosomes) > exclusion
%         temp_Propensities(Nribosomes+1,1) = k_bind * Inhibitor_Condition;
%     end
%     % elongation
%     elongation_Condition = (X_State(2:Nribosomes) + exclusion) < X_State(1:Nribosomes-1);
%     temp_Propensities(2:Nribosomes,1) = temp_Propensities(2:Nribosomes,1) .* elongation_Condition;
%
%
%
%     % deffining prop and stoich
%     propensities = temp_Propensities;
%     stoichiometry = temp_Stoichiometry;
%     %% Update time
%     sum_Propensities = sum(propensities);
%     t = t - log(rand) / sum_Propensities;
%     %% Generate output
%     while iteration <= number_TimePoints && t > t_array(iteration)
%         X_output(1:size(X_State,1),iteration) = X_State;
%         iteration = iteration + 1;
%     end
%
%
%     if t> 5000
%         counter_No_Ribosomes = counter_No_Ribosomes + sum (X_State==1);
%     end
%     %% Update state
%             old_X_State = X_State; % temporal value of x state to be compared with new value and determine if collisions are considered.
%
%     if t < t_final
%         %% Select reaction
%         r2 = sum_Propensities * rand;
%         i = 1;
%         tmp = propensities(i);
%         while tmp < r2
%             i = i + 1;
%             tmp = tmp + propensities(i);
%         end
%         %% update state
%         X_State(1:Nribosomes+1) = X_State(1:Nribosomes+1) + stoichiometry(:,i);
%     end
%
%
%     if t>5000 && (old_X_State ~= X_State)
%     collisions = collisions + nnz(~ elongation_Condition);
%     end
% end

