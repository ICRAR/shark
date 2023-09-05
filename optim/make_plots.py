import numpy as np
import os
import argparse
import matplotlib.pyplot as plt 
import pandas as pd

import plots
 
hobs = 0.7
h0 = 0.6751 # TODO: better not hard code these

def parse_args():
    parser = argparse.ArgumentParser(description='Proctor in prep. plots')
    parser.add_argument("-f", default='/scratch/pawsey0119/clagos/SHARK_Out/PSO83_Proctor/',
                         type=str, help="Optimisation folder name")
    args = parser.parse_args()
    f = args.f

    if not f.endswith("/"):
        f += "/"
    return f


def get_pso_info(pso_folder):
    '''Reads shark_pso.log to retrieve relevant information on PSO runs for plotting'''

    params_substr = "Search space parameters: "
    constr_substr_start = "Constraints:"
    constr_substr_end= "HPC mode:"
    n_substr = "Maximum iterations: "

    constr_list = []

    with open (pso_folder + "shark_pso.log", 'rt') as logs:

        const = False

        for line in logs:

            # extract parameters optimised 
            params_match_ind = line.find(params_substr)
            if params_match_ind!= -1: 
                params = np.array(line.split(params_substr)[1].strip('\n').split(" "))

            # extract number iterations performed 
            n_match_ind = line.find(n_substr)
            if n_match_ind!= -1: 
                n = int(line.split(n_substr)[1])    

            # store list of constraints
            if line.find(constr_substr_start) != -1: 
                const = True
                continue
            elif line.find(constr_substr_end) != -1:
                const = False
                continue
            elif const:
                constr_list.append(line.split(" ")[-1].strip('\n'))
    
    constr_list = np.array(constr_list)
    
    return params, n, constr_list


def flatten_list(l):
    return [item for sublist in l for item in sublist]                


def read_model_outputs(pso_folder, n_pso, plot_dir):
    '''Reads PSO info and model output for all PSO iterations,
      removes instances where Shark errored from data'''
    modelvals = []
    modelerrorvals = []
    pos = []
    fx = []

    for i in range(n_pso):
        fn_mv = pso_folder + "PSOSMF_" + str(i) + "/modelvals.npy"
        fn_mev = pso_folder + "PSOSMF_" + str(i) + "/modelerrorvals.npy"
        fn_pos = pso_folder + "tracks/track_0" + str(i).zfill(2)  +"_pos.npy"
        fn_fx = pso_folder + "tracks/track_0" + str(i).zfill(2) +"_fx.npy"

        modelvals.append(np.load(fn_mv, allow_pickle = True))
        modelerrorvals.append(np.load(fn_mev, allow_pickle = True))
        pos.append(np.load(fn_pos, allow_pickle = True))
        fx.append(np.load(fn_fx, allow_pickle = True))

    # plot log-likelihood vs. iteration
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    plots.plot_ll_iteration(fx, plot_dir)

    # flatten data so that length of lists = [n_pso x n_swarm]
    # (we no longer care about data on the level of individual particles)
    modelvals = flatten_list(modelvals)
    modelerrorvals = flatten_list(modelerrorvals)
    pos = np.array(flatten_list(pos))
    fx = np.array(flatten_list(fx))
    len_init = len(fx)

    # modelvals/modelerrorvals is nested list w/ dimensions: [n_pso_iterations x n_swarm][nconstr][len_constr]
    # pos/fx is list w/ dimensions: [n_pso_iterations x n_swarm] (independent of constraint number or length)

    # drop instances where parameter combination results in a Shark error
    # >= necessary for when there is more than one constraint (combined LL)
    fail_mask = fx >= 1e20
    pos = pos[~fail_mask]
    fx = fx[~fail_mask]
    modelvals = [i for (i, fail) in zip(modelvals, fail_mask) if not fail]
    modelerrorvals = [i for (i, fail) in zip(modelerrorvals, fail_mask) if not fail]
    # print("Dropped ", len_init - len(fx), "error instances")

    return modelvals, pos, fx


def get_best_fits(models, pos, fx, param_list, plot_dir):
    
    '''Returns best-fitting models for each observational constraint,
    and outputs values of free parameters that produces the best-fit.'''

    # get index of best fit
    # in practice we work with negative LL, so want the min value
    bf_ind = np.argmin(fx)
    bf_params = pos[bf_ind]
    bf_ll = -1*fx[bf_ind]

    # append loglikelihood for output
    param_list = np.append(param_list, "LOG LIKELIHOOD")
    bf_params = np.append(bf_params, bf_ll)

    df = pd.DataFrame(np.array([param_list, bf_params]).T, columns = ["param", "value"])
    with open(plot_dir + "best_fit_params.txt", 'w') as f:
        dfAsString = df.to_string(header=True, index=False)
        f.write(dfAsString)

    bf_mask = fx == np.min(fx)
    bf_models = [model for (model, bf) in zip(models, bf_mask) if bf][0]

    return bf_models


def load_obs(constr_list, obs_dir):
    '''Will load the relevant obs data for the list of constraints implemented'''
    #TODO: would be better not to have mass limits hard coded in

    smfz0_check = [True for constr in constr_list if "SMF_z0" in constr]
    rm_check = [True for constr in constr_list if "RM" in constr]

    obs_constr_dict = {}

    if smfz0_check:

        smfz0_tmp = np.loadtxt(obs_dir + "mf/SMF/SMF_Li2009.dat", unpack = True)

        smfz0 = []
        for i in range(len(smfz0_tmp)):
            tmp_var = smfz0_tmp[i]

            # correct for little h and mask to relevant mass range
            if i == 0: 
                tmp_var = tmp_var - 2.0 * np.log10(hobs) + 2.0 * np.log10(hobs/h0)
                mask = (tmp_var > 8.3) & (tmp_var < 12) # TODO:temp fix, this should be checked - using 8 doesn't work...
            if i ==1:
                tmp_var = tmp_var + 3.0 * np.log10(hobs) - 3.0 * np.log10(hobs/h0)

            smfz0.append(tmp_var[mask])

        obs_constr_dict['SMF_z0'] = smfz0

    if rm_check:

        rm_tmp = np.loadtxt(obs_dir + "SizeMass/GAMA_H-band_dlogM_0.25_reff.txt", unpack = True)
        mask = (rm_tmp[0] > 8) & (rm_tmp[0] < 12)


        rm = []
        for i in range(len(rm_tmp)):
            tmp_var = rm_tmp[i]
            rm.append(tmp_var[mask])

        obs_constr_dict['RM'] = rm

    return obs_constr_dict


def get_model_dict(models, constr_list):
    '''Takes cleaned modelvals.npy data as input (nested list).
    Returns dict of model vals that can be accessed by models["constraint"]'''
    # reshape models so that models[n] gives list of nth obs constraint 
    models = list(map(list, zip(*models)))

    models_dict = {}

    for i, constr in enumerate(constr_list):
        constr = constr[:constr.index("(")] 
        models_dict[constr] = models[i]

    return models_dict


def get_parameter_lookup():
    '''Creates dictionary that returns mathmode plot labels for each parameter
    TODO: log_10 labels need to be added where appropriate.
    Ideally shark_pos.log would contain list of whether a variable is logged or not to automate this.å'''

    param_lookup = {'stellar_feedback.beta_disk': r'$\beta$',
                    'stellar_feedback.min_beta': r'$\beta_{ {\rm min} }$',
                    'stellar_feedback.redshift_power': r'$z_{{\rm p}}$',
                    'stellar_feedback.v_sn': r'$v_{{\rm hot}}$',
                    'star_formation.nu_sf': r'$\nu_{{\rm SF}}$',
                    'star_formation.boost_starburst': r'$\nu_{{\rm burst}}$',
                    'agn_feedback.kappa_agn': r'$\kappa$',
                    'agn_feedback.kappa_radio': r'$\kappa_{{\rm radio}}$',
                    'agn_feedback.alpha_cool': r'$\alpha_{{\rm cool}}$',
                    'agn_feedback.mseed': r'$m_{{\rm seed}}$',
                    'agn_feedback.mhalo_seed': r'$m_{{\rm halo,\ seed}}$',
                    'agn_feedback.hot_halo_threshold': r'$\Gamma_{{\rm thresh}}$',
                    'reincorporation.tau_reinc': r'$\tau_{{\rm reinc}}$',
                    'reincorporation.mhalo_norm': r'$M_{{\rm norm}}$',
                    'reincorporation.halo_mass_power': r'$\gamma$',
                    'disk_instability.stable': r'$\epsilon_{{\rm disk}}$'}
    
    return param_lookup


def main():
    # set directories etc.
    pso_folder = parse_args()
    param_list, n, constr_list = get_pso_info(pso_folder)
    param_dict = get_parameter_lookup()
    plot_dir = pso_folder + "plots/"
    obs_dir = "../data/"

    print("Creating plots in ", plot_dir)
    print("These parameters: ", param_list)
    print("Were constrained against observations of: ", constr_list)
    print("For", n, "PSO iterations")

    # Load relevant observations
    obs_constr_dict = load_obs(constr_list, obs_dir)

    # Process model data and get best fit models
    models, pos, fx = read_model_outputs(pso_folder, n, plot_dir)
    bf_models = get_best_fits(models, pos, fx, param_list, plot_dir)
    models_dict = get_model_dict(models, constr_list)

    # useful diagnostic plots
    plots.plot_best_fits(models_dict, bf_models, obs_constr_dict, plot_dir)
    plots.plot_param_corrs(param_list, pos, fx, plot_dir)

    # make these plots for each observational constraint
    for constr in constr_list:
        constr = constr[:constr.index("(")] # removes mass limits
        print("Making plots for: ", constr)
        plots.plot_model_val_by_param(param_list, obs_constr_dict, pos, fx, models_dict,
                                      param_dict, plot_dir,
                                      param="stellar_feedback.v_sn", constr=constr)
        
        plots.plot_parameter_summaries(param_list, obs_constr_dict, pos, 
                                    models_dict, param_dict,
                                    plot_dir, constr=constr)
    

if __name__ == '__main__':
    main()