function [solutions, complete_times, ode_warn_flags, slope_norms] = ... 
    perturb_Ksets(Network_Data, All_K, All_fractions,t_interval, ...
    uptake_rxns, uptake_values, perturbed_rxn, expression_level)

    n_Ksets = size(All_K,2);
    n_rxns = length(Network_Data.rxns);
    n_metabs = size(Network_Data.S,1);
    enz_enzComplex = Network_Data.enz_enzComplex;
    enzymatic_rxns = find(Network_Data.rxn_type ~= 2);
    
    solutions = zeros(n_rxns,n_Ksets);
    complete_times = zeros(n_Ksets,1);
    ode_warn_flags = zeros(n_Ksets,1);
    slope_norms = zeros(n_Ksets,1);
    
    enzyme = find(enzymatic_rxns == perturbed_rxn);
    
    rxn_ind_tmp = find(enz_enzComplex(:,enzyme) ~=0);
    rxn_ind_tmp2 = [enzyme; rxn_ind_tmp] + n_metabs;
    
    %for x = 1:n_Ksets
    parfor x = 1:n_Ksets
        initial_conc_perturbed = [ones(n_metabs,1); All_fractions(:,x)];
        initial_conc_perturbed(rxn_ind_tmp2) = ...
            initial_conc_perturbed(rxn_ind_tmp2)*expression_level;
               
        [solutions(:,x), complete_times(x), ode_warn_flags(x), ...
            slope_norms(x)] = calculate_rates(Network_Data, t_interval, ...
            All_K(:,x), initial_conc_perturbed, uptake_rxns, uptake_values);
    end
         
end
