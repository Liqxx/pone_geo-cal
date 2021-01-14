# usual
import gc
import os
import sys
import copy
import numpy as np
import scipy.optimize as opt
from scipy.special import factorial
from scipy.interpolate import interp1d, interp2d

def llh(t, d,
        scale_pocam=1.0, delta_pocam=0.041,
        wt=1.0,
        chi2unc=False, chi2mod=False,
        chi2modunc=False, chi2modpen=False):
    '''
    Global main function to calculate LLH profiles
    based on POCAM ppc simulation data.
    
    Author: F. Henningsen
    Date: 2020-10-26
    '''
    # case handling
    chi2 = True if chi2unc or chi2mod or chi2modunc or chi2modpen else False

    # use poisson-binnd MLE
    if not chi2:
        
        # handle for loop usage
        if type(d) == np.int or type(d) == np.float:
            
            # use poisson if d < 20 else gaussian
            if d < 20.:
                
                l = -2*( d*np.log(t/wt)-np.log(np.math.factorial(d))-t/wt)
        
            else:
                
                l = (((d - t/wt)**2)/(t/wt))+ 2*np.log((np.sqrt(2*np.pi*t/wt)))
        
        # handle array usage
        else:
            m = (d < 20)
            l = copy.deepcopy(d)
            l[m]  = -2*( d*np.log(t/wt)-np.log(factorial(d))-t/wt)
            l[~m] = (((d - t/wt)**2)/(t/wt))+ 2*np.log((np.sqrt(2*np.pi*t/wt)))
    
    # use chi2 approximation, valid for bin counts n > 10
    else:
        
        # standard pearson chi2
        l = (t/wt - d)**2 / (t/wt)
        
        # standard pearson chi2 including POCAM intensity systematics
        if chi2unc:
            
            l = (t/wt - d)**2 / (t/wt + (delta_pocam*t)**2)
            
        # modified pearson chi2
        if chi2mod:
            
            l = (t/wt - d)**2 / (t/wt + np.sqrt(t/wt)**2)
            
        # modfied pearson chi2 including POCAM intensity systematics
        if chi2modunc:
            
            l = (t/wt - d)**2 / (t/wt + np.sqrt(t/wt)**2 + (delta_pocam*t/wt)**2)
            
        # modified pearson chi2 including penalty for POCAM systematics
        if chi2modpen:
            l  = (t/wt * scale_pocam - d)**2 / (t/wt * scale_pocam + t * (scale_pocam/wt)**2)
            #l += (scale_pocam - 1.00)**2 / (delta_pocam**2)
        
    # return llh for grid point
    return l


def llh_pocam_tot(scale,
                  truth, sim,
                  pocam, dom_keys,
                  wt=1,
                  delta_pocam=0.041,
                  min_q_tot=0, min_q_bin=10,
                  df=False,
                  smart_binning=True, n_data=30,
                  pocam_minimize=False,
                  chi2unc=False, chi2mod=False,
                  chi2modunc=False, chi2modpen=False):
    '''
    Function to calculate total LLH for one POCAM.
    '''
    
    # list for individual dom llhs
    llh_doms = np.zeros((len(dom_keys)))
    
    hit_doms = 0
    for k, key in enumerate(dom_keys):
        
        # list of llh for datapoints
        llh_datapoints = []
        
        # get om and string
        string, om = [float(i) for i in key.split('-')]
        
        # get om data
        if df:
            m_truth = (truth['om'] == om) & (truth['string'] == string)
            m_sim   = (sim['om'] == om) & (sim['string'] == string)
            truth_i = truth.loc[m_truth].values.squeeze()[:n_data]
            sim_i   = sim.loc[m_sim].values.squeeze()[:n_data]
        else:
            m_truth = (truth[:,n_data + 3] == string) & (truth[:,n_data + 4] == om)
            m_sim   = (sim[:,n_data + 3] == string) & (sim[:,n_data + 4] == om)
            truth_i = truth[m_truth].squeeze()[:n_data]
            sim_i   = sim[m_sim].squeeze()[:n_data]
                                
        if sum(truth_i) > min_q_tot and sum(sim_i) > min_q_tot: 

            hit_doms += 1

            if not smart_binning:

                index    = next((i for i, x in enumerate(truth) if x), None) # first non-zero
                truth    = truth_i[index:index+n_data]
                data_sim = sim_i[index:index+n_data]

            for i in range(n_data):
                
                t = float(truth_i[i])
                d = float(sim_i[i])

                # remove non-hit doms (and optionally apply some minimum bin value)
                if t > 0 and t > min_q_bin:
                    
                    # get global llh definition
                    llhi = llh(t, d,
                               scale_pocam=scale, delta_pocam=delta_pocam, wt=wt,
                               chi2unc=chi2unc, chi2mod=chi2mod, 
                               chi2modunc=chi2modunc, chi2modpen=chi2modpen)
                    
                    # append it 
                    llh_datapoints.append(llhi)
                    
            #print(key, truth_i, sim_i, llh_datapoints)
            #print()
            
            # append llh sum from DOM:
            llh_doms[k] = np.sum(llh_datapoints)
    
    # include penalty for chi2modpen
    penalty = (scale - 1.00)**2 / (delta_pocam**2) if chi2modpen else 0
    
    if pocam_minimize:
        return np.sum(llh_doms) + penalty
    
    else:
        return llh_doms, hit_doms

def llh_grid5d(df_t, df_d,
               t, d, 
               params,
               pocam_keys, dom_keys,
               scale_pocam=1.0, delta_pocam=0.041,
               min_q_tot=0, min_q_bin=10,
               df=False,
               store_indiv=False,
               smart_binning=True, n_data=30,
               pocam_minimize=False,
               chi2unc=False, chi2mod=False,
               chi2modunc=False, chi2modpen=True):
    '''
    Return 5-dimensional likelihood grid based on <data_binned> of 
    all five simulated parameters for POCAM PPC simulations.
    '''

    ### get parameters
    grid = params['grid']

    # data / scan
    scan_dict  = params['scan_dict']
    scan_n     = params['n_data_read']
    scan_N     = scan_dict['n']
    scan_Nph   = scan_dict['nph']
    scan_abs   = scan_dict['abs']
    scan_sca   = scan_dict['sca']
    scan_dome  = scan_dict['domeff']
    scan_p0    = scan_dict['p0']
    scan_p1    = scan_dict['p1']
    
    # truth
    truth_dict = params['truth_dict']
    t_N        = truth_dict['n']
    t_Nph      = truth_dict['nph']
    t_abs      = truth_dict['abs']
    t_sca      = truth_dict['sca']
    t_dome     = truth_dict['domeff']
    t_p0       = truth_dict['p0']
    t_p1       = truth_dict['p1']
        
    # combinations
    combinations = int(len(grid))
    
    # get weights for truth if statistics higher 
    wt = np.sqrt(t_N / scan_n)
    
    # Define empty llh-landscape for sum:
    llh_array = np.zeros((combinations,)) # grid [abs, sca, domeff, p0, p1]
    
    # Define empty llh-landscape for individual DOMs:
    string_doms = np.zeros((len(dom_keys), 2))
    for i, dk in enumerate(dom_keys):
        string_doms[i] = [int(dk.split('-')[0]), int(dk.split('-')[1])]
    strings      = string_doms[:,0]
    doms         = string_doms[:,1]
    llh_array_iv = np.zeros((len(pocam_keys), len(string_doms), combinations))

    combo_cntr = 0
    for gp in grid:
        
        # define parameter for data dictionary
        abso, sca, dome, p0, p1 = gp
        
        # print progress
        print('%i / %i (parameters: %.3f %.3f %.3f %.3f %.3f)' %(combo_cntr+1, combinations, abso, sca, dome, p0, p1))
                        
        # list of llh for pocams
        llh_pocams = []
                                    
        # whether to use dataframe or data_binned
        if df:
            # apply parameter mask
            m_truth = ()
            m_sim   = (df_d['abs'] == abso) & (df_d['sca'] == sca) & (df_d['dom_eff'] == dome) & (df_d['p0'] == p0) & (df_d['p1'] == p1)
            truth   = df_t.loc[m_truth]
            sim     = df_d.loc[m_sim]
        else:
            m_truth = ()
            m_sim   = (d[:,-5] == abso) & (d[:,-4] == sca) & (d[:,-3] == dome) & (d[:,-2] == p0) & (d[:,-1] == p1)
            truth   = t[m_truth]
            sim     = d[m_sim]
                                
        # format '<string>-<dom>'
        for k, pocam in enumerate(pocam_keys):
            
            # apply pocam mask
            pocam  = pocam.decode('utf-8')
            ps, po = [int(i) for i in pocam.split('-')]
            
            # whether to use dataframe or data_binned
            if df:
                m_truth = (truth['pocam'] == pocam)
                m_sim   = (sim['pocam'] == pocam)
                truth_i = truth.loc[m_truth]
                sim_i   = sim.loc[m_sim]
            else:
                m_truth = (truth[:,n_data+1] == ps) & (truth[:,n_data+2] == po)
                m_sim   = (sim[:,n_data+1] == ps) & (sim[:,n_data+2] == po)
                truth_i = truth[m_truth]
                sim_i   = sim[m_sim]
                                        
            # minimize systematic of POCAM scaling
            # define parameter dependent pocam function to minimize
            def llh_pocam_minimize(scale):
                llhsum = llh_pocam_tot(scale, 
                                       truth_i, sim_i,
                                       pocam, dom_keys,
                                       wt=wt,
                                       delta_pocam=delta_pocam,
                                       min_q_tot=min_q_tot, min_q_bin=min_q_bin,
                                       df=df,
                                       smart_binning=smart_binning, n_data=n_data,
                                       pocam_minimize=True,
                                       chi2unc=chi2unc, chi2mod=chi2mod, 
                                       chi2modunc=chi2modunc, chi2modpen=chi2modpen)
                return llhsum
            
            # and get minimal scale if penalty llh is use
            if chi2modpen:

                # minimize
                res = opt.minimize(llh_pocam_minimize, [1.0],
                                #method='Nelder-Mead', options={'fatol':0.01})
                                method='L-BFGS-B', bounds=[(1-5*delta_pocam, 1+5*delta_pocam)], options={'ftol':1e-3})
                scale = res['x'][0]
                
            else:
                scale = 1

            # and get final llh
            llh_doms, hit_doms= llh_pocam_tot(scale, 
                                              truth_i, sim_i,
                                              pocam, dom_keys,
                                              wt=wt,
                                              delta_pocam=delta_pocam,
                                              min_q_tot=min_q_tot, min_q_bin=min_q_bin,
                                              df=df,
                                              smart_binning=smart_binning, n_data=n_data,
                                              pocam_minimize=False,
                                              chi2unc=chi2unc, chi2mod=chi2mod, 
                                              chi2modunc=chi2modunc, chi2modpen=chi2modpen)

            # include penalty for chi2modpen
            penalty = (scale - 1.00)**2 / (delta_pocam**2) if chi2modpen else 0
            
            # append llh sum from POCAM
            llh_pocam_value = sum(llh_doms) / max([float(hit_doms), 1.]) + penalty
            llh_pocams.append(llh_pocam_value) 
            
            # print progress
            print('\tPOCAM %s --> LLH %.2f / scale %.5f' %(pocam, llh_pocam_value, scale))
            
            # optionally append individual llhs
            # index: pocam_idx, [strings x doms], [abs_idx x sca_idx x dome_idx x p0_idx x p1_idx]
            if store_indiv:
                llh_array_iv[k, :, combo_cntr] = llh_doms 

        # putting the final llh into the landscape array
        llh_array[combo_cntr] = np.sum(llh_pocams) / float(len(pocam_keys))
        
        # increment combinations counter
        combo_cntr += 1
        
        # collect garbage
        gc.collect()
    
    # return either llh sum or individual llh
    if not store_indiv:
        return llh_array
    
    else:
        return llh_array, llh_array_iv












































