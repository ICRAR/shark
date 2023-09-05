import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as cols

import pandas as pd
import scipy.stats as stats
import seaborn as sns

plt.style.use("mnras.mplstyle")

fig_w, fig_h = 3.321, 3
label_size = 14
tick_size = 12
legend_size = 12


def flatten_list(l):
    return [item for sublist in l for item in sublist] 


def plot_ll_iteration(fx, plot_dir):
    '''Plots the log-likelihood of each particle in the swarm 
    as a function of the iteration number'''

    # transpose so that each array is an individual particle
    fx = np.array(fx).T
    fig, ax = plt.subplots(1,1, figsize = (fig_w, fig_h))

    for particle in fx:
        particle[particle >= 1e20] = np.nan 
        # plot negative of fx, corresponds to log-likelihood value
        ax.plot(np.arange(len(particle))+1, -1*particle, alpha = .9)

    ax.set_ylim(-1*np.nanmax(flatten_list(fx)), 10+-1*np.nanmin(flatten_list(fx)))
    ax.set_xlabel("Iteration number")
    ax.set_ylabel("Log-likelihood")
    plt.savefig(plot_dir + "ll_iteration.pdf")
    plt.close()


def plot_best_fits(models, bf_models, obs_constrs, plot_dir):
    '''For each observtaional constraint, plots all models and the best-fit'''

    n_constr = len(models)
    fig, ax = plt.subplots(1, n_constr, figsize=(n_constr*fig_w, fig_h))

    # plot each constraint in separate panel
    for i, constr in enumerate(obs_constrs):
        model_tmp = models[constr]
        bf_tmp = bf_models[i]

        # get obs
        constr = obs_constrs[constr]
        xplot = constr[0] 
        yobs = constr[1]

        # for each PSO evaluation
        for j in range(len(model_tmp)):
            yplot = model_tmp[j]

            if n_constr > 1:
                ax[i].plot(xplot, yplot, color = "grey", lw = 1., alpha =0.2)
            else:    
                ax.plot(xplot, yplot, color = "grey", lw = 1., alpha =0.2)

        if n_constr > 1:        
            ax[i].plot(xplot, bf_tmp, color = "orange")
            ax[i].scatter(xplot, yobs, marker = 'o', s = 5, color = "k", zorder = 5)
            ax[i].set_xlabel(r"$\log\,M_\star / {\rm M_\odot}$")
        else:   
            ax.plot(xplot, bf_tmp, color = "orange") 
            ax.scatter(xplot, yobs, marker = 'o', s = 5, color = "k", zorder = 5)
            ax.set_xlabel(r"$\log\,M_\star / {\rm M_\odot}$")

    plt.savefig(plot_dir + "model_evaluations.pdf")
    plt.close()


def get_corr_dat(df):
    '''Get dataframe containing Spearman correlation coefficients between free parameters. 
    Takes pandas df of param values as input.'''

    cors = df.corr(method='spearman')
    cors = cors.where(np.triu(np.ones(cors.shape)).astype(bool))
    cors = cors.stack().reset_index()
    cors.columns = ['par1','par2','rho']

    return cors


def plot_param_corrs(param_list, pos, fx, plot_dir):

    # drop first 100 iterations - these usually confuse the correlation
    # discussed in Sec 4.2 of draft
    ndrop = 100
    fx = fx[ndrop:]
    pos = [params for i, params in enumerate(pos) if i >= ndrop]

    df = pd.DataFrame(pos, columns = param_list)
    corrs = get_corr_dat(df)

    # only plot significant correlations (abs(rho) > 0.35)
    sig_lim = 0.35
    plot_corrs = corrs[(abs(corrs['rho']) >= sig_lim) & (corrs['rho'] != 1)]
    plot_corrs.reset_index(inplace=True)

    # color points by LL value
    df['ll'] = fx*-1

    if(len(plot_corrs) == 0):
        print("Skipping sig_corrs.pdf - no signifcant correlations between parameters.")
        return

    if len(plot_corrs) % 2 == 0:
        nrows = int(len(plot_corrs)/2)
    else:
        nrows = int(len(plot_corrs)/2 + 0.5)

    fig, axs = plt.subplots(nrows, 2, figsize=(fig_w*2, 2*nrows))

    cmap = plt.cm.plasma_r
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=plt.Normalize(vmin=min(df['ll'] ),
                                                    vmax=max(df['ll'] )))
    
    for i, ax in enumerate(axs.flat): 

        # if odd number of vars, last plot should be empty
        if (i == 2*nrows-1) & (len(plot_corrs) % 2 != 0):
            ax.set_visible(False)
            continue

        par1 = plot_corrs['par1'][i]
        par2 = plot_corrs['par2'][i]
        
        ax.scatter(x=df[str(par1)], y=df[str(par2)], 
                c = df['ll'], cmap = 'plasma_r', s=35, ec = 'w', linewidths=0.5,
                    alpha = 0.75)
        ax.tick_params(direction='in', axis='both',
                    which='both', bottom=True,top=True,
                    left=True, right=True, labelsize = 16)
        ax.set_xlabel(str(par1), fontsize = label_size)
        ax.set_ylabel(str(par2), fontsize = 12)
        
        # spearman coeff
        label = str(r'$\rho = $' + str(round(plot_corrs['rho'][i],2)))
        ax.annotate(label,  xy=(0.7, 0.9), xycoords="axes fraction", 
                            size = label_size)
        ax.minorticks_on()
        
    plt.savefig(plot_dir + "sig_corrs.pdf")
    plt.close()
    

def label_corrs(x, y, **kwargs):
    coef = stats.spearmanr(x, y)[0]
    label = r'$\rho$ = ' + str(round(coef, 2))
    return label

def get_constr_ylab(constr):
    label_lookup = {'SMF_z0': r"$\Phi(M_{\star})/\rm{dlog M_{\star}/Mpc^{-3}}$",
                    'RM': r"$\log\,(r_{50}\,/\,{\rm{kpc}})$" }
    
    ylab = label_lookup[constr]
    return ylab


def plot_model_val_by_param(param_list, obs_dict, pos, fx, model_dict, param_dict, plot_dir,
                            param="agn_feedback.kappa_agn", constr="SMF_z0"):  
    '''Plots the values of a given free parameter against model values for 4 bins of stellar mass'''

    # select parameter and constraints
    param_ind = np.where(param_list == param)
    params = np.array(pos).T[param_ind]
    params = params[0]

    # relevant data for obs constraint
    models_plot = model_dict[constr]
    obs_plot = obs_dict[constr]
    ylab = get_constr_ylab(constr)  

    # only plot 4 example mass bins
    # equally spaced across the range, excluding extremes
    df = np.zeros(shape = (4, len(models_plot)))
    nobs = len(obs_plot[0])
    mass_bin_inds = np.linspace(1, nobs-1, 4).astype(int)
    constr_masses = obs_plot[0][mass_bin_inds]

    for i in range(len(models_plot)):
        model_massbins = models_plot[i]
        df[:, i] = model_massbins[mass_bin_inds]
    
    fig, axs = plt.subplots(2,2, figsize = (fig_w*2, fig_h*2), sharex=True)
    fig.subplots_adjust(hspace = .3, wspace=.25)
    cmap = plt.cm.plasma_r
    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=plt.Normalize(vmin=min(fx*-1),
                                                    vmax=max(fx*-1)))

    for i, ax in enumerate(axs.flat): 
        ax.scatter(params, df[i],
                       c = fx*-1,cmap =cmap, alpha =0.55, s=7)
        ax.set_title(r'$\log M_\star / {\rm M_{\odot}} = $'
                            +str(round(constr_masses[i], 1)),pad = 10, 
                            fontsize=label_size)
        ax.tick_params(direction='in', axis='both',
                    which='both', bottom=True,top=True,
                    left=True, right=True, labelsize = tick_size)
    
        # Add the spearman correlation to the plot
        ax.annotate(label_corrs(params, df[i]), size = label_size,
                            xy=(0.03, 0.05), xycoords='axes fraction')
    
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    fig.colorbar(sm, cax=cbar_ax)
    fig.subplots_adjust(right=0.8)
    cbar_ax.set_ylabel(ylabel = "Log-likelihood",size=12)
    fig.text(0.5, 0.04, param_dict[param], ha='center', fontsize = label_size)
    fig.text(0.04, 0.5, ylab, va='center', rotation='vertical', fontsize = label_size)

    plt.savefig(plot_dir + param.split('.')[1] + "_" + constr + '_massbins.pdf');    


def plot_parameter_summaries(param_list, obs_dict, pos, model_dict, param_dict,
                  plot_dir, constr="SMF_z0"):
    
    model = model_dict[constr]
    obs_plot = obs_dict[constr]
    constr_lab = get_constr_ylab(constr)  

    # reshape pos so len(pos) = n_params
    pos = np.array(pos).T    

    # get stellar mass values
    lm = obs_plot[0] 
    lm = np.round(lm-0.05, 1)

    # rho and fitted gradient dfs
    corr_data = np.zeros((len(pos), len(lm)))
    grad_data = np.zeros((len(pos), len(lm)))

    # labels for heatmap plot
    columns = lm
    rows = [param_dict[param] for param in param_list]
    
    # for each free parameter
    for param in range(len(pos)):
        paramvals = pos[param]
        
        # for each bin of log galaxy mass
        for massval in range(len(lm)):

            # get model values at mass bin
            y = list(zip(*model))[massval]

            # calc spearman coeff between free param and model vals
            corr = stats.spearmanr(paramvals, y)[0]
            corr_data[param, massval] = round(corr,2)

            # calc gradient of fitted line between free param and model vals
            slope = stats.linregress(paramvals, y)[0]
            grad_data[param, massval] = round(slope, 2)
                      
    corr_data = pd.DataFrame(corr_data,columns=columns,index=rows)
    grad_data = pd.DataFrame(grad_data,columns=columns,index=rows)

    ## FIG 1: rho heatmap
    grid_spec = {"width_ratios": (.9, .05)} # necessary for cbar formatting
    fig, (ax, cbar_ax) = plt.subplots(1,2, figsize=(2*fig_w, 2*fig_h), gridspec_kw=grid_spec) 
    ax = sns.heatmap(corr_data,
                      ax=ax, cbar_ax=cbar_ax,
                     annot_kws={"size": label_size},
                     linewidth=1, 
                     vmin = -1,
                     vmax = 1,
                     xticklabels = 3, # plot every 3rd tick label for readability
                     norm = cols.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1),
                     cmap='coolwarm',
                     cbar_kws={'label': r'$\rho (X, $' + constr_lab +'$)$'})
    
    ax.tick_params(labelsize=tick_size)  
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    cbar_ax.tick_params(labelsize = tick_size)  
    cbar_ax.yaxis.label.set_size(tick_size)
    ax.set_xlabel(r"$\log\,(M_{\star} / {\rm M}_{\odot})$",
               fontsize = label_size)
    plt.savefig(plot_dir + "corr_heatmap_" + constr + ".png")


    ## FIG 2: rho vs. mstar for each param
    indexNamesArr=corr_data.index.values

    # significant rho cut-off for plot mask
    rho_cut = 0.4
    sig_vars = corr_data[(abs(corr_data).max(axis=1) >= rho_cut)].index.values
    cmapvals = ['#219ebc', '#e56b6f', "#fb8500",
                '#007f5f', '#E9C46A', '#8338EC', 
                 '#00b4d8','#A4243B', '#B5838D']
    
    # plot param in grey if corr not significant, otherwise use cmapvals
    final_cols = []
    n_insig = 0
    for i,param in enumerate(indexNamesArr):
        if param in sig_vars:
            col_ind = i - n_insig
            final_cols.append(cmapvals[col_ind])
        else:
            n_insig += 1
            final_cols.append('grey')
    coldict = dict(zip(indexNamesArr, final_cols))    
    
    fig, ax = plt.subplots(1,1,figsize=(2*fig_w, 2*fig_h))
    
    for param in range(len(indexNamesArr)):
        
        dat = corr_data.iloc[param]
                
        # significance threshold
        below = abs(dat.T) < rho_cut
        insig = dat.T
        sig = np.ma.masked_where(below, dat.T) # masks insignificant values
        
        ax.plot(lm, sig, 
                    label = indexNamesArr[param], 
                    color = coldict[str(indexNamesArr[param])],
                   linewidth = 2.5)
        
        ax.plot(lm, insig, 
                   alpha=0.4,
                    color = coldict[str(indexNamesArr[param])],
                     linestyle='--',
                   linewidth = 2.)
        
    ax.set_ylim(-1,1)
    ax.legend(loc='best', prop={'size':legend_size},ncol=3,frameon=False)

    ax.set_xlabel(r"$\log\,(M_{\star}/ {\rm M}_{\odot})$",
          fontsize = label_size)
    ax.set_ylabel(r'$\rho (X, $' + constr_lab +'$)$',
          fontsize = label_size)
    ax.tick_params(direction='in', axis='both',
                 which='both', bottom=True,top=True,
                left=True, right=True, labelsize = tick_size)
    ax.axhline(0, color = 'grey', linestyle = '--')
    plt.savefig(plot_dir + "corr_param_mstar_" + constr + ".png")  


    ## FIG 3: fitted gradient vs. mstar for significant params only
    indexNamesArr= grad_data.index.values  
    
    fig, ax = plt.subplots(1,1,figsize=(2*fig_w, 2*fig_h))
    
    ylim = 0

    for param in sig_vars:
        
        grad = grad_data.loc[param]
        sig = corr_data.loc[param]
                
        # significance threshold - still mask by corr_data
        below = abs(sig.T) < rho_cut
        insig = grad.T
        sig = np.ma.masked_where(below, grad.T) # masks insignificant grad values

        ylim = np.max([ylim, np.max(np.abs(np.array(grad)))])

        ax.plot(lm, sig, 
                    label = param, 
                    color = coldict[param],
                   linewidth = 2.5)
        
        ax.plot(lm, insig, 
                   alpha=0.4,
                    color = coldict[param],
                     linestyle='--',
                   linewidth = 2.)
        
    ax.set_ylim(-ylim-0.1, ylim+0.1)
    ax.legend(loc='best', prop={'size':legend_size},
              ncol=3,frameon=False)

    ax.set_xlabel(r"$\log\,(M_{\star}/ {\rm M}_{\odot})$",
          fontsize = label_size)
    ax.set_ylabel(r'$\alpha (X, $' + constr_lab +'$)$',
          fontsize = label_size)
    ax.tick_params(direction='in', axis='both',
                 which='both', bottom=True,top=True,
                left=True, right=True, labelsize = tick_size)
    ax.axhline(0, color = 'grey', linestyle = '--')
    plt.savefig(plot_dir + "fitted_grad_mstar_" + constr + ".png")  