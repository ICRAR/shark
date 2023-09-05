#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2019
# Copyright by UWA (in the framework of the ICRAR)
#
# Originally contributed by Mawson Sammons
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""
Constraints for optimizers to evaluate shark models against observations
"""

import logging
import os
import re

import analysis
import common
import numpy as np
import smf
import sizes
import sizes_with_icl
import global_quantities


logger = logging.getLogger(__name__)

GyrToYr = 1e9

#######################
# Binning configuration
mlow = 5
mupp = 14
dm = 0.125
mbins = np.arange(mlow,mupp,dm)
xmf = mbins + dm/2.0

mlow2 = 5
mupp2 = 14
dm2 = 0.3
mbins2 = np.arange(mlow2,mupp2,dm2)
xmf2 = mbins2 + dm2/2.0

ssfrlow = -6
ssfrupp = 4
dssfr = 0.2
ssfrbins = np.arange(ssfrlow,ssfrupp,dssfr)


sfrlow = -3
sfrupp = 1.5
dsfr = 0.2
sfrbins = np.arange(sfrlow,sfrupp,dsfr)
xsfr    = sfrbins + dsfr/2.0


# These are two easily create variables of these different shapes without
# actually storing a reference ourselves; we don't need it
zeros1 = lambda: np.zeros(shape=(1, 3, len(xmf)))
zeros2 = lambda: np.zeros(shape=(1, 3, len(xmf2)))
zeros3 = lambda: np.zeros(shape=(1, len(mbins)))
zeros4 = lambda: np.empty(shape=(1), dtype=np.bool_)
zeros5 = lambda: np.zeros(shape=(1, 4, len(ssfrbins)))
zeros6 = lambda: np.zeros(shape = (1, 3, len(xsfr)))

class Constraint(object):
    """Base classes for constraint objects"""

    convert_to_multiple_batches = True

    def __init__(self):
        self.redshift_table = None
        self.weight = 1
        self.rel_weight = 1

    def _load_model_data(self, modeldir, subvols):

        if  len(subvols) > 1 and self.convert_to_multiple_batches:
            subvols = ["multiple_batches"]

        # Histograms we are interested in
        hist_smf = zeros3()
        hist_HImf = zeros3()
        hist_smf_err = zeros3()
        hist_HImf_err = zeros3()
        hist_smf_comp = np.zeros(shape = (1, 4, len(mbins)))

        fields = {
            'galaxies': (
                'sfr_disk', 'sfr_burst', 'mstars_disk', 'mstars_bulge',
                'rstar_disk', 'm_bh', 'matom_disk', 'mmol_disk', 'mgas_disk',
                'matom_bulge', 'mmol_bulge', 'mgas_bulge',
                'mgas_metals_disk', 'mgas_metals_bulge',
                'mstars_metals_disk', 'mstars_metals_bulge', 'type',
                'mvir_hosthalo', 'rstar_bulge', 'mstars_burst_mergers', 
                'mstars_burst_diskinstabilities', 'mstars_bulge_mergers_assembly', 
                'mstars_bulge_diskins_assembly', 'specific_angular_momentum_disk_gas_atom',
                'vmax_subhalo', 'rgas_disk')

        }

        for index, z in enumerate(self.z):
            hdf5_data = common.read_data(modeldir, self.redshift_table[z], fields, subvols)
            h0 = hdf5_data[0]
            smf.prepare_data(hdf5_data, index, hist_smf, zeros3(), zeros3(),
                             zeros3(), zeros3(), hist_HImf, zeros3(), zeros3(),
                             zeros3(), zeros3(), zeros3(), zeros3(), zeros1(), zeros1(),
                             zeros1(), zeros1(), zeros1(), zeros1(), zeros1(),
                             zeros1(), zeros1(), zeros1(), zeros1(), zeros1(),
                             zeros1(), zeros4(), zeros4(), zeros2(), zeros5(),
                             zeros1(), zeros1(), zeros1(), zeros1(), zeros1(), 
                             zeros1(), hist_smf_err, hist_HImf_err, hist_smf_comp, 
                             zeros6())

        #########################
        # take logs
        ind = np.where(hist_smf > 0.)
        hist_smf_err[ind] = abs(np.log10(hist_smf[ind]) - np.log10(hist_smf_err[ind]))
        hist_smf[ind] = np.log10(hist_smf[ind])

        ind = np.where(hist_HImf > 0.)
        hist_HImf_err[ind] = abs(np.log10(hist_HImf[ind]) - np.log10(hist_HImf_err[ind]))
        hist_HImf[ind] = np.log10(hist_HImf[ind])

   	    #### Read CSFR model data ####
        fields = {'global': ('redshifts', 'm_hi', 'm_h2', 'mcold', 'mcold_metals',
                         'mhot_halo', 'mejected_halo', 'mbar_lost', 'mbar_created', 'mstars', 
                         'mstars_bursts_mergers', 'mstars_bursts_diskinstabilities',
                         'm_bh', 'sfr_quiescent', 'sfr_burst', 'm_dm', 'mcold_halo', 'number_major_mergers', 
                         'number_minor_mergers', 'number_disk_instabilities', 'smbh_maximum')}

    	# Read data from each subvolume at a time and add it up
    	# rather than appending it all together
        for idx, subvol in enumerate(subvols):
            subvol_data = common.read_data(modeldir, self.redshift_table[0], fields, [subvol])
            max_bhs_subvol = subvol_data[20].copy()
            if idx == 0:
                hdf5_data_sfr        = subvol_data
                max_smbh         = max_bhs_subvol
            else:
                max_smbh = np.maximum(max_smbh, max_bhs_subvol)
                for subvol_datum, hdf5_datum in zip(subvol_data[3:], hdf5_data_sfr[3:]):
                    hdf5_datum += subvol_datum
                    #select the most massive black hole from the last list item

        # Also make sure that the total volume takes into account the number of subvolumes read
        hdf5_data_sfr[1] = hdf5_data_sfr[1] * len(subvols)

        h0_sfr, redshifts = hdf5_data_sfr[0], hdf5_data_sfr[2]

        (mstar_plot, mcold_plot, mhot_plot, meje_plot,
        mstar_dm_plot, mcold_dm_plot, mhot_dm_plot, meje_dm_plot, mbar_dm_plot,
        sfr, sfrd, sfrb, mstarden, mstarbden_mergers, mstarbden_diskins, sfre, sfreH2, mhrat,
        mHI_plot, mH2_plot, mH2den, mdustden, omegaHI, mdustden_mol, mcoldden, mhotden,
        mejeden, history_interactions, mDMden, mlost_dm_plot, mcreated_dm_plot) = global_quantities.prepare_data(hdf5_data_sfr, redshifts)

	#### Size-mass relation ####
        mlow3 = 6.5
        mupp3 = 12.5
        dm3 = 0.2
        mbins3 = np.arange(mlow3,mupp3,dm3)
        xmf_rm = mbins3 + dm3/2.0
  
        dmobs = 0.4
        mbins_obs = np.arange(mlow,mupp,dmobs)
        xmf_obs = mbins_obs + dmobs/2.0

        vlow = 1.0
        vupp = 3.0
        dv   = 0.1
        vbins = np.arange(vlow, vupp, dv)
        xv_rm = vbins + dv/2.0
 	
        fields_icl = {'galaxies': ('mstars_disk', 'mstars_bulge', 'rstar_disk', 'rstar_bulge', 'type',
                    'mstellar_halo', 'cnfw_subhalo', 'vvir_hosthalo', 'mvir_hosthalo')}
	
        # Loop over redshift and subvolumes
        rcomb = np.zeros(shape = (len(self.z), 3, len(xmf_rm)))
        rcomb_icl = np.zeros(shape = (len(self.z), 3, len(xmf_rm)))
        bs_error = np.zeros(shape = (len(self.z), len(xmf_rm)))

        for index, z in enumerate(self.z):
            hdf5_data = common.read_data(modeldir, self.redshift_table[z], fields_icl, subvols)
            # constrain to sizes without icl component (rcomb) 
            sizes_with_icl.prepare_data(hdf5_data, index, rcomb, rcomb_icl, bs_error)
        # change rcomb to rcomb_icl below to constrain to sizes+icl component (10*r50)
        return h0, hist_smf, hist_HImf, hist_smf_err, hist_HImf_err, redshifts, sfr, xmf_rm, rcomb, bs_error

    def load_observation(self, *args, **kwargs):
        obsdir = os.path.normpath(os.path.abspath(os.path.join(__file__, '..', '..', 'data')))
        return common.load_observation(obsdir, *args, **kwargs)

    def _get_raw_data(self, modeldir, subvols):
        """Gets the model and observational data for further analysis.
        The model data is interpolated to match the observation's X values."""

        h0, hist_smf, hist_HImf, hist_smf_err, hist_HImf_err, redshifts, sfr, xmf_rm, rcomb, bs_error = self._load_model_data(modeldir, subvols)
        x_obs, y_obs, y_dn, y_up = self.get_obs_x_y_err(h0)
        x_mod, y_mod, y_mod_err = self.get_model_x_y(h0, hist_smf, hist_smf_err, hist_HImf, hist_HImf_err, redshifts, sfr, xmf_rm, rcomb, bs_error)
        return x_obs, y_obs, y_dn, y_up, x_mod, y_mod, y_mod_err

    def get_data(self, modeldir, subvols, plot_outputdir=None):

        x_obs, y_obs, y_dn, y_up, x_mod, y_mod, y_mod_err = self._get_raw_data(modeldir, subvols)

        # Both observations and model values don't come necessarily in order,
        # but if at the end of the day we want to perform array-wise operations
        # over them (to calculate chi2 or student-t) then they should be both in
        # ascending order
        sorted_obs = np.argsort(x_obs)
        x_obs = x_obs[sorted_obs]
        y_obs = y_obs[sorted_obs]
        y_dn = y_dn[sorted_obs]
        y_up = y_up[sorted_obs]
        sorted_mod = np.argsort(x_mod)
        x_mod = x_mod[sorted_mod]
        y_mod = y_mod[sorted_mod]
        y_mod_err = y_mod_err[sorted_mod]

        # Linearly interpolate model Y values respect to the observations'
        # X values, and only take those within the domain. We do the same 
        # for the errors of the model.
        # We also consider the biggest relative error as "the" error, in case
        # they are different and add it in quadrature to the Poisson error
        # of the model.
        y_mod_interp = np.interp(x_obs, x_mod, y_mod)
        y_mod_err_interp = np.interp(x_obs, x_mod, y_mod_err)
        sel = np.where((x_obs >= self.domain[0]) & (x_obs <= self.domain[1]))
        x_obs_sel = x_obs[sel]
        y_obs_sel = y_obs[sel]
        y_mod_sel = y_mod_interp[sel]
        y_mod_err_sel = y_mod_err_interp[sel]
        y_obs_err_sel = np.maximum(np.abs(y_dn[sel]), np.abs(y_up[sel]))
        err = np.sqrt( y_obs_err_sel ** 2.0 + y_mod_err_sel ** 2.0)

        if plot_outputdir:
            self.plot(plot_outputdir,
                      x_obs, y_obs, y_dn, y_up,
                      x_mod, y_mod, y_mod_interp,
                      x_obs_sel, y_obs_sel, y_mod_sel, err)

        return y_obs_sel, y_mod_sel, err

    def plot(self, plot_outputdir, x_obs, y_obs, obs_err_dn, obs_err_up,
             x_mod, y_mod, y_mod_interp, x_obs_sel, y_obs_sel, y_mod_sel, err):
        fig = common.load_matplotlib().figure(figsize=(4.5,4.5))
        ax = fig.add_subplot(111)
        ax.axvline(self.domain[0], ls='dotted', c='red')
        ax.axvline(self.domain[1], ls='dotted', c='red')
        ax.plot(x_obs_sel, y_obs_sel, marker='v', ls='None', c='blue', label="Selected observations")
        ax.plot(x_mod, y_mod, marker='^', ls='solid', c='orange', label="Model")
        ax.plot(x_obs, y_mod_interp, ls='solid', c='green', label="Interpolated model")
        ax.plot(x_obs_sel, y_mod_sel, ls='solid', c='brown', label="Selected model")
        common.errorbars(ax, x_obs, y_obs, obs_err_dn, obs_err_up, 'black', '+',
                         err_absolute=False, label="Observations")

        common.prepare_legend(ax, ['blue', 'orange', 'green', 'brown', 'black'])

        chi2 = analysis.chi2(y_obs_sel, y_mod_sel, err)
        st = analysis.studentT(y_obs_sel, y_mod_sel, err)
        ax.set_title('%s\nchi2 = %g, student-t = %g' % (str(self), chi2, st))

        common.savefig(plot_outputdir, fig, str(self))

    def __str__(self):
        s = '%s(%.1f-%.1f)'
        args = self.__class__.__name__, self.domain[0], self.domain[1]
        if self.weight != 1:
            s += ', weight=%.2f, rel_weight=%.2f'
            args += self.weight, self.rel_weight
        return s % args


class HIMF(Constraint):
    """The HI Mass Function constraint"""

    domain = (7, 12)
    z = [0]

    def get_obs_x_y_err(self, h0):
        # Load Zwaan05 data and correct data for their choice of cosmology
        lmHI, pHI, dpHIdn, dpHIup = self.load_observation('mf/GasMF/HIMF_Zwaan2005.dat', [0,1,2,3])

        #correct data for their choice of cosmology
        hobs = 0.75
        x_obs = lmHI + np.log10(pow(hobs,2)/pow(h0,2))
        y_obs = pHI + np.log10(pow(h0,3)/pow(hobs,3))
        y_dn = dpHIdn
        y_up = dpHIup
        #lmHI, pHI, pdnHI, pduHI = self.load_observation('mf/GasMF/HIMF_Jones18.dat', cols=[0,1,2,3])
        #dpdnHI = pHI - pdnHI
        #dpupHI = pduHI - pHI
        #hobs = 0.7
        #x_obs = lmHI + np.log10(pow(hobs, 2) / pow(h0, 2))
        #y_obs = pHI + np.log10(pow(h0, 3) / pow(hobs, 3))

        return x_obs, y_obs, y_dn, y_up

    def get_model_x_y(self, _, __, ___, hist_HImf, hist_HImf_err, ____, _____, ______, _______, ________):
        y = hist_HImf[0]
        yerr = hist_HImf_err[0]
        ind = np.where(y < 0.)
        return xmf[ind], y[ind], yerr[ind]

class SMF(Constraint):
    """Common logic for SMF constraints"""

    domain = (8, 13)

    def get_model_x_y(self, _, hist_smf, hist_smf_err, __, ___, ____, _____, ______, _______, ________):
        y = hist_smf[0,:]
        yerr = hist_smf_err[0,:]
        ind = np.where(y < 0.)
        return xmf[ind], y[ind], yerr[ind]

class SMF_z0(SMF):
    """The SMF constraint at z=0"""

    z = [0]

    def get_obs_x_y_err(self, h0):

        #lm, p, dp = self.load_observation('mf/SMF/GAMAIV_Driver22.dat', cols=[0,1,2])
        #hobs = 0.7
        #x_obs = lm + 2.0 * np.log10(hobs/h0)
        #y_obs = p - 3.0 * np.log10(hobs/h0)
        #y_dn = dp
        #y_up = dp

        lm, p, dpdn, dpup = self.load_observation('mf/SMF/SMF_Li2009.dat', cols=[0,1,2,3])
        hobs = 0.7
        x_obs = lm - 2.0 * np.log10(hobs) + 2.0 * np.log10(hobs/h0)
        y_obs = p + 3.0 * np.log10(hobs) - 3.0 * np.log10(hobs/h0)
        y_dn = dpdn
        y_up = dpup

        return x_obs, y_obs, y_dn, y_up

class SMF_z0p5(SMF):
    """The SMF constraint at z=0.5"""

    z = [0.5]

    def get_obs_x_y_err(self, h0):

        # Wright et al. (2018, several reshifts). Assumes Chabrier IMF.
        #zD17, lmD17, pD17, dp_dn_D17, dp_up_D17, dp_cv = self.load_observation('mf/SMF/Wright18_CombinedSMF.dat', cols=[0,1,2,3,4,5])
        #hobs = 0.7
        #binobs = 0.25
        #pD17 = pD17 + 2.0 * np.log10(hobs/h0) - np.log10(binobs)
        #lmD17 = lmD17 - 3.0 * np.log10(hobs/h0) - np.log10(binobs)
        #in_redshift = np.where(zD17 == 0.5)

        #x_obs = lmD17[in_redshift]
        #y_obs = pD17[in_redshift]
        #y_dn = dp_dn_D17[in_redshift]
        #y_up = dp_up_D17[in_redshift]
        #cv = np.log10(1 + dp_cv[in_redshift]) 

        # combine cosmic variance and model variance in quadrature
        #y_dn = np.sqrt(y_dn**2 + cv**2)
        #y_up = np.sqrt(y_up**2 + cv**2)

        #SMF from Weaver et al. (2022)
        lm, pD, dn, du = self.load_observation('mf/SMF/COSMOS2020/SMF_Farmer_v2.1_0.5z0.8_total.txt', cols = [0,2,3,4])
        hobs = 0.7
        y_obs = np.log10(pD) +  3.0 * np.log10(hobs/h0)
        y_dn = np.log10(pD) - np.log10(dn)
        y_up = np.log10(du) - np.log10(pD)
        x_obs = lm -  2.0 * np.log10(hobs/h0)

        return x_obs, y_obs, y_dn, y_up

class SMF_z1(SMF):
    """The SMF constraint at z=1"""

    z = [1]

    def get_obs_x_y_err(self, h0):

        # Wright et al. (2018, several reshifts). Assumes Chabrier IMF.

        #zD17, lmD17, pD17, dp_dn_D17, dp_up_D17, dp_cv = self.load_observation('mf/SMF/Wright18_CombinedSMF.dat', cols=[0,1,2,3,4,5])
        #hobs = 0.7
        #binobs = 0.25
        #pD17 = pD17 + 2.0 * np.log10(hobs/h0) - np.log10(binobs)
        #lmD17 = lmD17 - 3.0 * np.log10(hobs/h0) - np.log10(binobs)
        #in_redshift = np.where(zD17 == 1)
        #x_obs = lmD17[in_redshift]
        #y_obs = pD17[in_redshift]
        #y_dn = dp_dn_D17[in_redshift]
        #y_up = dp_up_D17[in_redshift]
        #cv = np.log10(1 + dp_cv[in_redshift])

        ## combine cosmic variance and model variance in quadrature
        #y_dn = np.sqrt(y_dn**2 + cv**2)
        #y_up = np.sqrt(y_up**2 + cv**2)

        ## take first set of values
        #ind = np.arange(0,12)
        #x_obs = x_obs[ind]
        #y_obs = y_obs[ind]
        #y_dn = y_dn[ind]
        #y_up = y_up[ind]

        #SMF from Weaver et al. (2022)
        lm, pD, dn, du = self.load_observation('mf/SMF/COSMOS2020/SMF_Farmer_v2.1_0.8z1.1_total.txt', cols = [0,2,3,4])
        hobs = 0.7
        y_obs = np.log10(pD) +  3.0 * np.log10(hobs/h0)
        y_dn = np.log10(pD) - np.log10(dn)
        y_up = np.log10(du) - np.log10(pD)
        x_obs = lm -  2.0 * np.log10(hobs/h0)

        return x_obs, y_obs, y_dn, y_up
 
class SMF_z2(SMF):
    """The SMF constraint at z=2"""

    z = [2]

    def get_obs_x_y_err(self, h0):

        #SMF from Weaver et al. (2022)
        lm, pD, dn, du = self.load_observation('mf/SMF/COSMOS2020/SMF_Farmer_v2.1_2.0z2.5_total.txt', cols = [0,2,3,4])
        hobs = 0.7
        y_obs = np.log10(pD) +  3.0 * np.log10(hobs/h0)
        y_dn = np.log10(pD) - np.log10(dn)
        y_up = np.log10(du) - np.log10(pD)
        x_obs = lm -  2.0 * np.log10(hobs/h0)

        return x_obs, y_obs, y_dn, y_up


class CSFR(Constraint):
    """The Cosmic Star Formation Rate constraint"""

    domain = (0, 5)
    z = [0]

    def get_obs_x_y_err(self, h0):
	#Driver (Chabrier IMF), ['Baldry+2012, z<0.06']
        redD17d, redD17u, sfrD17, err1, err2, err3, err4 = self.load_observation('Global/Driver18_sfr.dat', [0,1,2,3,4,5,6])

        hobs = 0.7
        xobsD17 = (redD17d+redD17u)/2.0
        yobsD17 = sfrD17 + np.log10(hobs/h0)
        errD17 = yobsD17*0. - 999.
        errD17 = np.sqrt(pow(err1,2.0)+pow(err2,2.0)+pow(err3,2.0)+pow(err4,2.0))
        y_dn, y_up = errD17, errD17 # symmetric error
        x_obs = xobsD17
        y_obs = yobsD17

        return x_obs, y_obs, y_dn, y_up

    def get_model_x_y(self, h0,_, __, ___, ____, redshifts, sfr, _____, ______, _______):
	#note that only h^2 is needed because the volume provides h^3, and the SFR h^-1.
        ind = np.where(sfr > 0)
        y = np.log10(sfr[ind]*pow(h0,2.0))
        xz = redshifts[ind]

        # y error - zero for now to avoid breaking things
        yerr = y*0.0
        
        return xz, y, yerr

class RM(Constraint):
    """The combined size-mass relation"""
    domain = (7,12)
    z = [0]

    def get_obs_x_y_err(self, h0):
        # use the semi-major axis sizes
        x_obs, y_obs, count, y_err = self.load_observation('SizeMass/GAMA_H-band_dlogM_0.25_reff.txt', [0,1,2,3])
        y_dn = y_err
        y_up = y_err

        return x_obs, y_obs, y_dn, y_up

    def get_model_x_y(self, _, __, ___, ____, _____, ______, _______, xmf_rm, rcomb,bs_error):

        ind = np.where(rcomb[0,0,:] != 0)
        x_mod = xmf_rm[ind]
        y_mod = rcomb[0,0,ind][0]
        y_err = bs_error[0][ind]

        return x_mod, y_mod, y_err

def _evaluate(constraint, stat_test, modeldir, subvols, plot_outputdir):
    try:
        y_obs, y_mod, err = constraint.get_data(modeldir, subvols,
                                                plot_outputdir=plot_outputdir)
        return stat_test(y_obs, y_mod, err) * constraint.weight
    except:
        logger.exception('Error while evaluating constraint, returning 1e20')
        return 1e20


def evaluate(constraints, stat_test, modeldir, subvols, plot_outputdir=None):
    """Returns the evaluation of all constraints, as a total number (default)
    or as individual numbers for each constraint"""
    return [_evaluate(c, stat_test, modeldir, subvols, plot_outputdir)
            for c in constraints]

def _get_y_mod(constraint, modeldir, subvols, plot_outputdir):
    try:
        y_obs, y_mod, err = constraint.get_data(modeldir, subvols,
						plot_outputdir=plot_outputdir)
        return y_mod
    except:
        logger.exception('Error')
        return [1e20]

def get_y_mod(constraints, modeldir, subvols, plot_outputdir=None):
    return[_get_y_mod(c, modeldir, subvols, plot_outputdir)
            for c in constraints]


def _get_y_err(constraint, modeldir, subvols, plot_outputdir):
    try:
        y_obs, y_mod, err = constraint.get_data(modeldir, subvols,
                                                plot_outputdir=plot_outputdir)
        return err
    except:
        logger.exception('Error')
        return [1e20]

def get_y_err(constraints, modeldir, subvols, plot_outputdir=None):
    return[_get_y_err(c, modeldir, subvols, plot_outputdir)
                for c in constraints]


def log_results(constraints, results):
    """Emits a log message showing the function evaluation for `constraints`"""

    sums = [sum(result) for result in results]
    min_idx = min(enumerate(sums), key=lambda enumerated_sum: enumerated_sum[1])[0]
    min_flags = [idx == min_idx for idx in range(len(sums))]

    n_cols = len(constraints) + 1
    msg = 'Particle evaluation results per-particle, per-constraint:\n'
    args = ()
    msg += ' ' * 3 + ' '.join(["%20.20s"] * n_cols) + ' Min\n'
    args += tuple(constraints)
    args += 'Total',
    msg += ' ' * 3 + ' '.join(["=" * 20] * n_cols) + ' ===\n'
    for particle_num, (result, min_flag) in enumerate(zip(results, min_flags)):
        msg += '%2d' + ' ' + ' '.join(['%20e'] * len(constraints)) + ' %20e %2s\n'
        args += particle_num,
        args += tuple(result)
        args += sum(result),
        args += '*' if min_flag else '',
    logger.info(msg, *args)


_constraint_re = re.compile((r'([0-9_a-zA-Z]+)' # name
                              '(?:\(([0-9\.]+)-([0-9\.]+)\))?' # domain boundaries
                              '(?:\*([0-9\.]+))?')) # weight
def parse(spec):
    """Parses a comma-separated string of constraint names into a list of
    Constraint objects. Specific domain values can be specified in `spec`"""

    _constraints = {
        'HIMF': HIMF,
        'SMF_z0': SMF_z0,
	'SMF_z0p5': SMF_z0p5,
        'SMF_z1': SMF_z1,
	'CSFR': CSFR,
	'RM': RM
    }

    def _parse(s):
        m = _constraint_re.match(s)
        if not m or m.group(1) not in _constraints:
            raise ValueError('Constraint does not specify a valid constraint: %s' % s)
        c = _constraints[m.group(1)]()
        if m.group(2):
            dn, up = float(m.group(2)), float(m.group(3))
            if dn < c.domain[0]:
                raise ValueError('Constraint low boundary is lower than lowest value possible (%f < %f)' % (dn, c.domain[0]))
            if up > c.domain[1]:
                raise ValueError('Constraint up boundary is higher than lowest value possible (%f > %f)' % (up, c.domain[1]))
            c.domain = (dn, up)
        if m.group(4):
            c.weight = float(m.group(4))
        return c

    constraints = [_parse(s) for s in spec.split(',')]
    total_weight = sum([c.weight for c in constraints])
    for c in constraints:
        c.rel_weight = c.weight / total_weight
    return constraints
