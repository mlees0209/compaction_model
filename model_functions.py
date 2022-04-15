#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 16:42:49 2020

Contains functions required for the compaction modelling. 

@author: mlees
"""
print('Importing required modules and functions; please wait.')
import sys
import re
import distutils.util as ut
import os
from time import process_time 
from matplotlib.dates import date2num
from matplotlib.dates import num2date
from datetime import datetime as dt
import numpy as np
import platform
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import scipy
import glob
import math

### Define parameter default values. We will then read and overwrite any from the parameter file.

Defaults_Dictionary={'internal_time_delay':True,'overwrite':True,'run_name':False,'output_folder':False,'no_layers':True,'layer_names':True,'layer_types':True,'layer_thicknesses':True,'layer_compaction_switch':True,'interbeds_switch':True,'interbeds_type':False,'clay_Ssk_type':False,'clay_Ssk':False,'sand_Ssk':True,'compressibility_of_water':True,'dz_clays':True,'dt_gwaterflow':True,'create_output_head_video':True,'overburden_stress_gwflow':True,'overburden_stress_compaction':True,'rho_w':True,'g':True,'specific_yield':True,'save_effective_stress':True,'time_unit':True,'compaction_solver_debug_include_endnodes':True,'save_internal_compaction':True,'mode':True,'resume_directory':False,'resume_date':False,'layer_thickness_types':True,'compaction_solver_compressibility_type':False,'save_s':True,'initial_stress_type':False,'initial_stress_offset':False,'initial_stress_offset_unit':True,'n_z':True,'CTL_value':True} # Define which variables have pre-defined defaults

Default_Values={'internal_time_delay':0.5,'overwrite':False,'no_layers':2,'layer_names':['Upper Aquifer', 'Lower Aquifer'],'layer_types':{'Upper Aquifer': 'Aquifer', 'Lower Aquifer': 'Aquifer'},'layer_thicknesses':{'Upper Aquifer': 100.0,'Lower Aquifer': 100.0},'layer_compaction_switch':{'Upper Aquifer': True, 'Lower Aquifer': True},'interbeds_switch':{'Upper Aquifer': False, 'Lower Aquifer': False},'sand_Ssk':1,'compressibility_of_water':4.4e-10,'dz_clays':0.3,'dt_gwaterflow':1,'create_output_head_video':False,'overburden_stress_gwflow':False,'overburden_stress_compaction':False,'rho_w':1000,'g':9.81,'specific_yield':0.2,'save_effective_stress':False,'time_unit':'days','compaction_solver_debug_include_endnodes':False,'save_internal_compaction':False,'mode':'Normal','layer_thickness_types':'constant','save_s':False,'initial_stress_offset_unit':'stress','n_z':12,'CTL_value':0.5}


class Logger(object):
    def __init__(self,output_folder,run_name):
        self.terminal = sys.stdout
        self.log = open("%s/%s/logfile.log" % (output_folder,run_name), "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

def round_to_multiple(number, multiple, direction='nearest'):
    if direction == 'nearest':
        return multiple * round(number / multiple)
    elif direction == 'up':
        return multiple * math.ceil(number / multiple)
    elif direction == 'down':
        return multiple * math.floor(number / multiple)
    else:
        return multiple * round(number / multiple)


def solve_head_equation_singlevalue(dt,t,dx,x,bc,ic,k):        
    print ( '' )
    print ( '\t\t\tFD1D_HEAT_EXPLICIT_SINGLEVALUE_Ssk:' )
    print ( '\t\t\t  Compute an approximate solution to the time-dependent' )
    print ( '\t\t\t  one dimensional heat equation:' )
    print ( '' )
    print ( '\t\t\t    dH/dt - K * d2H/dx2 = f(x,t)' )
    print ( '' )
    
    cfl = k * dt / dx / dx # This is the coefficient that determines convergence 

    if ( 0.5 <= cfl ):
        print ( '\t\t\tFD1D_HEAT_EXPLICIT_CFL - Fatal error!' )
        print ( '\t\t\t  CFL condition failed.' )
        print ( '\t\t\t  0.5 <= K * dT / dX / dX = %f' % ( cfl ) )
        sys.exit(1)

    
    print ( '\t\t\t  Number of X nodes = %d' % ( len(x) ) )
    print ( '\t\t\t  X interval is [%f,%f]' % ( min(x), max(x) ) )
    print ( '\t\t\t  X spacing is %f' % ( dx ) )
    print ( '\t\t\t  Number of T values = %d' % ( len(t) ) )
    print ( '\t\t\t  T interval is [%f,%f]' % ( min(t), max(t) ) )
    print ( '\t\t\t  T spacing is %f' % ( dt ) )
    print ( '\t\t\t  Constant K = %g' % ( k ) )
    print ( '\t\t\t  CFL coefficient = %g' % ( cfl ) )

    hmat = np.zeros ( ( len(x), len(t) ) ) # running the code creates an output matrix of heads at each position and time

    for j in range ( 0, len(t) ):
        if j % (int(len(t)/20)) == 0:
            printProgressBar(j,len(t))
        if ( j == 0 ):
            h = ic*np.ones(len(x))
            h[0] = bc[0,0]
            h[-1] = bc[1,0]
        else:
            h_new = np.zeros ( len(x) )

            for c in range ( 1, len(x) - 1 ):
                l = c - 1
                r = c + 1
                h_new[c] = h[c] + cfl * ( h[l] - 2.0 * h[c] + h[r] )
            h_new[0] = bc[0,j]
            h_new[-1] = bc[1,j]
            h = h_new
        for i in range ( 0, len(x) ):
            hmat[i,j] = h[i]
    print(' ')
    print('\t\t\t SOLVER COMPLETE')
    return hmat


def solve_head_equation_elasticinelastic(dt,t,dx,x,bc,ic,k_elastic,k_inelastic,overburdenstress=False,overburden_data=[],preset_initial_maxstress=False,initial_maxstress=[]):        
    '''This is the main solver of the diffusion equation for effective stress, with an elastic/inelastic switch for K. Description of (some-incomplete) options:
        - t is the array of times over which to solve (t = np.array([t0, t1, t2, t3, ..., tn]). Entries to t should be equally spaced with spacing dt. 
        - preset_initial_maxstress: BOOL flag to indicate whether the initial stress condition and the preconsolidation stress are equivalent.
        - initial_condition_maxstress: FLOAT ARRAY: Only used if preset_initial_maxstress is false, in which case it gives the initial stress value above which inelastic behaviour will occur. Useful if you are resuming a model run.
    '''
    
    print ( '' )
    print ( '\t\t\tFD1D_HEAT_EXPLICIT_ELASTICINELASTIC_Ssk:' )
    print ( '\t\t\t  Python version: %s' % ( platform.python_version ( ) ) )
    print ( '\t\t\t  Compute an approximate solution to the time-dependent' )
    print ( '\t\t\t  one dimensional heat equation:' )
    print ( '' )
    print ( '\t\t\t    dH/dt - K * d2H/dx2 = f(x,t)' )
    print ( '' )
    if overburdenstress:
        print(' \t\t\tSOLVING WITH OVERBURDEN STRESS INCLUDED;  ')
        print('\t\t\tOverburden data read in as ')
        print(overburden_data)
    else:
        overburden_data = np.zeros_like(t)
    
    h_precons = np.zeros ( ( len(x), len(t) ) )
    if not preset_initial_maxstress:
        h_precons[:,0] = ic
    else:
        print('\t\t\t preset max stress found to be',initial_maxstress)
        h_precons[:,0]=initial_maxstress
    inelastic_flag = np.zeros ( ( len(x), len(t) ) ) 

    cfl_elastic = k_elastic * dt / dx / dx # This is the coefficient that determines convergence 
    cfl_inelastic = k_inelastic * dt / dx / dx # This is the coefficient that determines convergence 



    
    print ( '\t\t\t  Number of X nodes = %d' % ( len(x) ) )
    print ( '\t\t\t  X interval is [%f,%f]' % ( min(x), max(x) ) )
    print ( '\t\t\t  X spacing is %f' % ( dx ) )
    print ( '\t\t\t  Number of T values = %d' % ( len(t) ) )
    print ( '\t\t\t  T interval is [%f,%f]' % ( min(t), max(t) ) )
    print ( '\t\t\t  T spacing is %f' % ( dt ) )
    print ( '\t\t\t  Elastic K = %g' % ( k_elastic ) )
    print ( '\t\t\t  Inelastic K = %g' % ( k_inelastic ) )
    print ( '\t\t\t  CFL elastic coefficient = %g' % ( cfl_elastic ) )
    print ( '\t\t\t  CFL inelastic coefficient = %g' % ( cfl_inelastic ) )
    
    if ( 0.5 <= max(cfl_elastic,cfl_inelastic) ):
        print ( '\t\t\tFD1D_HEAT_EXPLICIT_CFL - Potentially Fatal Warning!' )
        print ( '\t\t\t  CFL condition failed.' )
        print ( '\t\t\t  0.5 <= K * dT / dX / dX = %f' % max(cfl_elastic,cfl_inelastic))
        print ( "\t\t\t THIS MAY MEAN THAT YOUR SOLUTION HAS LARGE ERRORS!")


    hmat = np.zeros ( ( len(x), len(t) ) ) # running the code creates an output matrix of heads at each position and time

    for j in range ( 0, len(t) ):
        if j % (int(len(t)/20)) == 0:
            printProgressBar(j,len(t))
        if ( j == 0 ):
            h = ic*np.ones(len(x))
            h[0] = bc[0,0]
            h[-1] = bc[1,0]
        else:
            h_new = np.zeros ( len(x) )

            for c in range ( 1, len(x) - 1 ):
                l = c - 1
                r = c + 1
                if inelastic_flag[c,j-1]:
                    h_new[c] = h[c] + cfl_inelastic * ( h[l] - 2.0 * h[c] + h[r] ) + (overburden_data[j] - overburden_data[j-1])
                    #print(dt * overburden_data[j])
#                    h_new[c] = h[c] + cfl * ( h[l] - 2.0 * h[c] + h[r] ) + dt * f[c] This is the forcing version

                elif not inelastic_flag[c,j-1]:
                    h_new[c] = h[c] + cfl_elastic * ( h[l] - 2.0 * h[c] + h[r] ) + (overburden_data[j] - overburden_data[j-1])
                    #print(dt * overburden_data[j])

                else:
                    print('Uh oh! Neither inelastic or elastic..something went wrong.')

            h_new[0] = bc[0,j]
            h_new[-1] = bc[1,j]
            h = h_new
        for i in range ( 0, len(x) ):
            hmat[i,j] = h[i]
            if h[i] - overburden_data[j] < h_precons[i,j]:
                if j <= len(t)-3:
                    h_precons[i,j+1]=h[i] - overburden_data[j]
                    inelastic_flag[i,j] =1
            else:
                if j <= len(t)-3:
                    h_precons[i,j+1] = h_precons[i,j]
            if j == len(t)-3:
                if h[i] - overburden_data[j] < h_precons[i,j]:
                    inelastic_flag[i,j] =1

#        if j <= len(t)-3:
#            h_precons[:,j+1] = np.min(hmat[:,:j+1],axis=1)

    print(' ')
    print('\t\t\t SOLVER COMPLETE')
    return hmat, inelastic_flag



def read_parameter(name,typ,length,paramfilelines,printlots=True):
    '''Reads in a parameter from paramfilelines. 
    name=parameter name
    typ=parameter type (str,int,float). if length>1 this is the type that each entry in the list/dictionary will be. Note that special value 'dict' is for clay interbeds only.
    length=how many comma separated entries are there for the parameter'''
        
    par_paramfile=[x for x in paramfilelines if x.replace(' ','').startswith('%s=' % name)]
    if len(par_paramfile)==0:
        if Defaults_Dictionary[name]==False:
            print("\tReading parameters: NOTICE. No '%s' found in parameter file and no default set." % name)
            return

        else:
            print("\t\tReading parameters error: WARNING. No '%s' found in parameter file. Using default value." % name)
            par=Default_Values[name]
            if printlots:
                print('\t%s=%s' % (name,par))
            return par
        
    par=par_paramfile[0].split('#')[0].split('=')[1]
                     
    if length==1:
        par = par.strip()
    else:
        if ':' not in par:
            par=par.split(',')
            if len(par) != length:
                print('\t\tReading parameters error: terminal. %s should have %i entries but only has %i.' % (name,length,len(par)))
                sys.exit(1)
            par=[x.strip() for x in par]
    
    if ':' in par:
        if typ!=dict:
            par=re.split(':|,',par)
            if len(par) != length*2:
                print('\tERROR: %s=%s' % (name,par))          
                print('\t\tReading parameters error: terminal. %s should have %i entries but has %i.' % (name,length,0.5*len(par)))
                sys.exit(1)
        else:
            par = re.split('},',par)
            b = [re.split(':{',p.replace('}','').strip()) for p in par]
            return b
        if typ!=bool:
            par=dict([(par[2*i].strip(),typ(par[2*i+1].strip())) for i in range(length)])
        else:
            par=dict([(par[2*i].strip(),bool(ut.strtobool(par[2*i+1].strip()))) for i in range(length)])

    if length==1:
        if type(par) != dict:
            if typ == bool:
                par=bool(ut.strtobool(par)) 
            else: 
                par=typ(par)

    if printlots:
        print('\t%s=%s' % (name,par))
    
    if length==1:
        if type(par) != dict:
            if type(par) != typ:
                print("\t\tReading parameters error: WARNING. %s is of class %s but may need to be %s. This may lead to errors later." % (name, type(par),typ))
        else:
            if type(par)!=dict:
                iscorrect=[type(x)==typ for x in par]
                if False in iscorrect:
                    print("\t\tReading parameters error: WARNING. elements of %s should be %s but 1 or more is not. This may lead to errors later." % (name,typ))
    
    if len(par_paramfile)>1:
        print("\t\tReading parameters error: WARNING. Multiple '%s's found in parameter file. using the first." % name)
        
    if type(par)==dict:
        iscorrect=[type(x) ==typ for x in par.values()]
        if False in iscorrect:
            print("\t\tReading parameters error: WARNING. values of %s should be %s but 1 or more is not. This may lead to errors later." % (name,typ))

    return par

def read_parameter_layerthickness_multitype(name,paramfilelines,printlots=True):
      par_out = {}
      par_paramfile=[x for x in paramfilelines if x.replace(' ','').startswith('%s=' % name)]
      par=par_paramfile[0].split('#')[0].split('=')[1]
      par=re.split(':|,',par)
      # Find if there are any dictionaries
      contains_dict = ['{' in s or '}' in s for s in par]
      if sum(contains_dict) % 2 != 0:
          print('\tERROR: odd number of { or } brackets found.' % (name,par))          
          print('\t\tReading parameters error: terminal. %s should have even number of { } brackets.' % name)
      elif sum(contains_dict) > 0:
          for i_tmp in range(int(len(np.where(contains_dict)[0])/2)):
              layername_tmp=par[np.where(contains_dict)[0][2*i_tmp]-1]
              par_out[layername_tmp] = {}
              dic_length_tmp = int((np.where(contains_dict)[0][2*i_tmp + 1] - np.where(contains_dict)[0][2*i_tmp] + 1)/2)
              for j_tmp in range(dic_length_tmp):
                  keytmp = par[np.where(contains_dict)[0][2*i_tmp]+ j_tmp*2].split('{')[-1]
                  val_tmp = par[np.where(contains_dict)[0][2*i_tmp]+ 1+ j_tmp*2].split('}')[0]
                  par_out[layername_tmp][keytmp] = float(val_tmp)
      dict_idxs_full=[]
      for i_tmp in range(int(len(np.where(contains_dict)[0])/2)):
          dict_idxs_full= np.append(dict_idxs_full,np.arange(np.where(contains_dict)[0][2*i_tmp]-1,np.where(contains_dict)[0][2*i_tmp+1]+1))
      all_idxs = np.arange(len(par))
      nondict_idxs = [idx for idx in all_idxs if idx not in dict_idxs_full]
      nondict_par = np.array(par)[nondict_idxs]
      for i_tmp in range(int(len(nondict_par)/2)):
          par_out[nondict_par[2*i_tmp]]=float(nondict_par[2*i_tmp+1])
      if printlots:
          print('\t%s=%s' % (name,par_out))
      return par_out

def subsidence_solver_aquitard_elasticinelastic(hmat,Sske,Sskv,b0,n_z,TESTn=1,overburden=False,unconfined=False,overburden_data=0,debuglevel=0,endnodes=False,preset_initial_maxstress=False,ic_maxstress=[]):
    ''' TESTn is a temporary variable, referring to the number of midpoints done. If you start with 20 nodes and TESTn=1, you integrate over 20 nodes. If TESTn=2 you intergrate over 40 nodes, and so on. It can be used to reduce error from using the Riemann sum. NOTE AS OF Apr 2022, this no longer works.'''
    print('Running subsidence solver. Overburden=%s, unconfined=%s.' % (overburden,unconfined))
    if overburden:
        if not unconfined:
            print(' \t\t\tSOLVING WITH OVERBURDEN STRESS INCLUDED;  ')
            print('\t\t\tOverburden data read in as ')
            print(overburden_data)
        else:
            print(' \t\t\tSOLVING WITH OVERBURDEN STRESS INCLUDED;  ')
            print('\t\t\tThis aquifer is unconfined. ')
            print(overburden_data)

    print('Aquitard solver is done at midpoints. Applying linear interpolation to hmat.')
    if not endnodes:
        hmat_interp = np.zeros((np.shape(hmat)[0]*(2*TESTn)-(2*TESTn-1),np.shape(hmat)[1]))
        for i in range(np.shape(hmat)[1]):
            if i % (int(np.shape(hmat)[1]/20)) == 0:
                printProgressBar(i,np.shape(hmat)[1])
            # if len(hmat[:,i]) != len( 0.000000001*np.arange(0,1000000000*np.shape(hmat)[0]*dz,1000000000*dz)):
            #     print('ERROR: hmat is not the same length as 0.001*np.arange(0,1000*np.shape(hmat)[0]*dz,1000*dz). If dz_clays is not a multiple of the layer thickness, you may need to give it to more significant figures for this to work.')
            #     print(0.000000001*np.arange(0,1000000000*np.shape(hmat)[0]*dz,1000000000*dz))
            #     sys.exit(1)
            a = scipy.interpolate.interp1d(np.linspace(0,b0,n_z),hmat[:,i],kind='linear')
#            print(np.arange(0,np.shape(hmat)[0]*dz,dz))
##            print(np.shape(hmat_interp)[0])
##            print(dz)
##            print(np.shape(a))
#            print(np.arange(0,np.shape(hmat_interp)[0]*(dz/2),(dz/2)))
            hmat_interp[:,i] = a(np.linspace(0,b0,2*n_z-1))
        if TESTn != 1:
            hmat_midpoints = hmat_interp[1:-1,:]
        else:
            hmat_midpoints = hmat_interp[1::2,:]
    else:
        print('DEBUG: not using midpoints for head solver.')
        hmat_midpoints = hmat
    #hmat_midpoints_precons = np.array([np.min(hmat_midpoints[:,:i+1],axis=1) for i in range(np.shape(hmat_midpoints)[1])]).T
    
    if overburden:
        overburden_data_midpoints = np.tile(overburden_data, (np.shape(hmat_midpoints)[0],1))
    else:
        overburden_data_midpoints = np.zeros_like(hmat_midpoints)
    
    stress_midpoints = overburden_data_midpoints - hmat_midpoints
    
    stress_midpoints_precons = np.zeros_like(hmat_midpoints)
    inelastic_flag_midpoints = np.zeros_like(hmat_midpoints)
    if preset_initial_maxstress:
        print('preset precons found. interpolating to midpoints.')
        print('starting with',ic_maxstress)
        a = scipy.interpolate.interp1d(np.linspace(0,b0,n_z),ic_maxstress)
        ic_precons_interp = a(np.linspace(0,b0,n_z*2-1))
        ic_precons_initial = ic_precons_interp[1::2]
        print('interpoolated to',ic_precons_initial)
        stress_midpoints_precons[:,0] = ic_precons_initial
    else:
        stress_midpoints_precons[:,0] = overburden_data_midpoints[:,0] - hmat_midpoints[:,0]
    
    for i in range(np.shape(stress_midpoints)[1]-1):
        if i % (int((np.shape(stress_midpoints)[1]-1)/20)) == 0:
            printProgressBar(i,np.shape(stress_midpoints)[1])        
        for j in range(np.shape(stress_midpoints)[0]):
            if stress_midpoints[j,i] > stress_midpoints_precons[j,i]:
                stress_midpoints_precons[j,i+1]=stress_midpoints[j,i]
                inelastic_flag_midpoints[j,i]=1
            else:
                stress_midpoints_precons[j,i+1]=stress_midpoints_precons[j,i]
                inelastic_flag_midpoints[j,i]=0

    if debuglevel==1:
        plt.figure()
        plt.imshow(stress_midpoints,aspect='auto')
        plt.colorbar()
        plt.title('Stress at midpoints')

        plt.figure()
        plt.imshow(inelastic_flag_midpoints,aspect='auto')
        plt.colorbar()
        plt.title('Inelastic flag at midpoints')
        
        plt.figure()
        plt.plot(stress_midpoints[10,:],'-',label='stress')
        plt.plot(hmat_midpoints[10,:],'-',label='head')
        plt.legend()
        ax2 = plt.twinx(plt.gca())
        ax2.plot(inelastic_flag_midpoints[10,:],'k--',label='inelastic flag')
        plt.title('stress head and inelastic flag at node 10')
        plt.show()
    
    inelastic_flag_midpoints= np.array(inelastic_flag_midpoints,dtype=bool)    

    # print('doing db')
    # db = [dz*( (inelastic_flag_midpoints[:,i] * Sskv * (stress_midpoints[:,i+1] - stress_midpoints[:,i])) + ~inelastic_flag_midpoints[:,i] * Sske * (stress_midpoints[:,i+1] - stress_midpoints[:,i])) for i in range(np.shape(stress_midpoints)[1]-1)]
    # print('doing ds')
    # ds = [dz*( np.dot(inelastic_flag_midpoints[:,i] * Sskv, stress_midpoints[:,i+1] - stress_midpoints[:,i]) + np.dot(~inelastic_flag_midpoints[:,i] * Sske, stress_midpoints[:,i+1] - stress_midpoints[:,i])) for i in range(np.shape(stress_midpoints)[1]-1)]
    # print('doing ds elastic')
    # ds_elastic = [dz*(np.dot(~inelastic_flag_midpoints[:,i] * Sske, stress_midpoints[:,i+1] - stress_midpoints[:,i])) for i in range(np.shape(stress_midpoints)[1]-1)]
    # print('doing ds inelastic')
    # ds_inelastic = [dz*( np.dot(inelastic_flag_midpoints[:,i] * Sskv, stress_midpoints[:,i+1] - stress_midpoints[:,i]))  for i in range(np.shape(stress_midpoints)[1]-1)]

    # s = np.zeros(np.shape(hmat)[1])
    # s_elastic = np.zeros(np.shape(hmat)[1])
    # s_inelastic = np.zeros(np.shape(hmat)[1])

    # print('\tIntegrating deformation over time.')
    # for i in range(1,np.shape(hmat)[1]):
    #     if i % (int(np.shape(hmat)[1]/20)) == 0:
    #         printProgressBar(i,np.shape(hmat)[1]-1)
    #     s[i] = s[i-1]-ds[i-1]
    #     s_elastic[i] = s_elastic[i-1]-ds_elastic[i-1]
    #     s_inelastic[i] = s_inelastic[i-1]-ds_inelastic[i-1]

    b = np.zeros(np.shape(hmat)[1])
    for ti in range(1,np.shape(hmat)[1]):
        if ti % (int(np.shape(hmat)[1]/20)) == 0:
            printProgressBar(ti,np.shape(hmat)[1]-1)
        b[ti] = b0/n_z * ( Sskv * np.sum(stress_midpoints_precons[:,ti] - stress_midpoints_precons[:,0]) - Sske * np.sum(stress_midpoints_precons[:,ti] - stress_midpoints[:,ti]))
    
    b = -1 * np.array(b)

#    return db,s,s_elastic,s_inelastic,inelastic_flag_midpoints
    return b,inelastic_flag_midpoints

def create_head_video_elasticinelastic(hmat,z,inelastic_flag,dates_str,outputfolder,layer,delt=100,startyear=None,endyear=None,datelabels='year'):
    # I think delt is in units of days; see what happens with the variable t_jumps below. startyear and endyear should be YYYY. datelabels can be 'year' or 'monthyear' and changes the title of the plots.
        if not os.path.isdir('%s/headvideo_%s' % (outputfolder, layer.replace(' ','_'))):
            os.mkdir('%s/headvideo_%s' % (outputfolder, layer.replace(' ','_')))
        
        t1_start = process_time()  
 
        # Make the video frames; use ffmpeg -f image2 -i %*.png vid.mp4 to make the vid itself
        plt.ioff()
        sns.set_context('talk')
        print('\t\tCreating frames.')
        plt.figure(figsize=(6,6))
        plt.xlabel('Head in the clay (m)')
        plt.ylabel('Z (m)')
        plt.xlim([np.max(hmat)+1,np.min(hmat)-1])
        plt.ylim([np.min(z)-0.1*(z[-1]-z[1]),np.max(z)+0.1*(z[-1]-z[1])])
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()
        t_tmp = date2num([dt.strptime(date, '%d-%b-%Y').date() for date in dates_str])
        
        t_jumps = [t_tmp[i] - t_tmp[0] for i in range(len(t_tmp))] # It looks like we get t in days.
        if not delt in t_jumps:
            print("\tERROR MAKING VIDEO! Selected dt not compatible with dates given.")
            sys.exit(1)
        
        if startyear:
            starting_t=date2num(dt.strptime('01-09-%s' % startyear, '%d-%m-%Y'))      
            print('Movie starting at date %s' % num2date(starting_t).strftime('%d-%b-%Y'))
        else:
            starting_t = t_tmp[0]
            print('Movie starting at date %s' % num2date(starting_t).strftime('%d-%b-%Y'))


        if endyear:
            ending_t=date2num(dt.strptime('01-09-%s' % endyear, '%d-%m-%Y'))      
            print('Movie ending at date %s' % num2date(ending_t).strftime('%d-%b-%Y'))
        else:
            ending_t = t_tmp[-1]
            print('Movie ending at date %s' % num2date(ending_t).strftime('%d-%b-%Y'))

        
        ts_to_do = np.arange(starting_t,ending_t+0.0001,delt)
        
        firsttime=0
        
        for t in ts_to_do:
            i=np.argwhere(t_tmp==t)[0][0] # If there are multiple frames on the same day, this line means we take the first of them
            if firsttime==2: # This whole "firsttime" bit is just to help print every 5%, it's really not important.
                if i % (int(len(ts_to_do)/20)) <= di: # This should make it so only 20 progress bars are printed
                    printProgressBar(np.argwhere(ts_to_do==t)[0][0],len(ts_to_do))
            else:
                printProgressBar(np.argwhere(ts_to_do==t)[0][0],len(ts_to_do))

            #plt.plot(hmat[:,0],np.linspace(0,40,20),'k--',label='t=0 position')
            if firsttime==1:
                firsttime=2
                di = np.argwhere(ts_to_do==t)[0][0]
            if t==starting_t:
                l1, = plt.plot(np.min(hmat[:,:i+1],axis=1),z,'b--',label='Preconsolidation head')
                l2,= plt.plot(hmat[:,i],z,'r.')
                l3, = plt.plot(hmat[:,i][~inelastic_flag[:,i]],z[~inelastic_flag[:,i]],'g.')
                plt.legend()
                firsttime=1
            else:
                l1.set_xdata(np.min(hmat[:,:i+1],axis=1))
                l2.set_xdata(hmat[:,i])
                l3.remove()
                l3, = plt.plot(hmat[:,i][~inelastic_flag[:,i]],z[~inelastic_flag[:,i]],'g.')
            if datelabels=='year':
                plt.title('t=%s' % num2date(t_tmp[i]).strftime('%Y'))
            if datelabels=='monthyear':
                plt.title('t=%s' % num2date(t_tmp[i]).strftime('%b-%Y'))

#            set_size(plt.gcf(), (12, 12))
            plt.savefig('%s/headvideo_%s/frame%06d.jpg' % (outputfolder, layer.replace(' ','_'),i),dpi=60,bbox_inches='tight')
        t1_stop = process_time() 
        print("")     
        print("\t\tElapsed time in seconds:",  t1_stop-t1_start)  
        print('')
        print('\t\tStitching frames together using ffmpeg.')
        #cmd='ffmpeg -hide_banner -loglevel warning -r 10 -f image2 -i %*.jpg vid.mp4'
        cmd='ffmpeg -hide_banner -loglevel warning -r 10 -f image2 -i %%*.jpg -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" vid_years%sto%s.mp4' % (num2date(starting_t).strftime('%Y'),num2date(ending_t).strftime('%Y'))

        print('\t\t\t%s.' % cmd)
        cd=os.getcwd()
        os.chdir('%s/headvideo_%s' % (outputfolder, layer.replace(' ','_')))
        subprocess.call(cmd,shell=True)
        
        if os.path.isfile('vid_years%sto%s.mp4' % (num2date(starting_t).strftime('%Y'),num2date(ending_t).strftime('%Y'))):
              print('\tVideo seems to have been a success; deleting excess .jpg files.')
              jpg_files_tmp = glob.glob('*.jpg')
              for file in jpg_files_tmp:
                  os.remove(file)
              print('\tDone.')

        os.chdir(cd)
        
        
def create_compaction_video(outputfolder,layer,db_plot,z,time_sim,Inelastic_Flag,delt=20,startyear=None,endyear=None,datelabels='year'):
    
    t1_start = process_time()  
    
    if not os.path.isdir('%s/compactionvideo_%s' % (outputfolder, layer.replace(' ','_'))):
        os.mkdir('%s/compactionvideo_%s' % (outputfolder, layer.replace(' ','_')))


    plt.ioff()
    sns.set_context('talk')
    print('\t\tCreating frames.')
    plt.figure(figsize=(8,8))
    plt.xlabel('node compaction rate (cm yr-1)')
    plt.ylabel('Z (m)')
    plt.xlim([np.max(db_plot),np.min(db_plot)])
    plt.ylim([np.min(z)-0.1*(z[-1]-z[1]),np.max(z)+0.1*(z[-1]-z[1])])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.xscale('symlog')
    t_tmp = np.round(time_sim,5)
    
    t_jumps = np.round([t_tmp[i] - t_tmp[0] for i in range(len(t_tmp))],4) # It looks like we get t in days.
    if not delt in t_jumps:
        print("\tERROR MAKING VIDEO! Selected dt not compatible with dates given.")
        sys.exit(1)
    
    if startyear:
        starting_t=date2num(dt.strptime('01-09-%s' % startyear, '%d-%m-%Y'))      
        print('Movie starting at date %s' % num2date(starting_t).strftime('%d-%b-%Y'))
    else:
        starting_t = t_tmp[0]
        print('Movie starting at date %s' % num2date(starting_t).strftime('%d-%b-%Y'))
    
    
    if endyear:
        ending_t=date2num(dt.strptime('01-09-%s' % endyear, '%d-%m-%Y'))      
        print('Movie ending at date %s' % num2date(ending_t).strftime('%d-%b-%Y'))
    else:
        ending_t = t_tmp[-1]
        print('Movie ending at date %s' % num2date(ending_t).strftime('%d-%b-%Y'))


    ts_to_do = 1/1000* np.arange(1000*starting_t,1000*ending_t+0.0001,1000*delt) # avoid weird rounding issues
    
    firsttime=0
    
    for t in ts_to_do:
        i=np.argwhere(t_tmp==t)[0][0] # If there are multiple frames on the same day, this line means we take the first of them
        if firsttime==2: # This whole "firsttime" bit is just to help print every 5%, it's really not important.
            if i % (int(len(ts_to_do)/20)) <= di: # This should make it so only 20 progress bars are printed
                printProgressBar(np.argwhere(ts_to_do==t)[0][0],len(ts_to_do))
        else:
            printProgressBar(np.argwhere(ts_to_do==t)[0][0],len(ts_to_do))
    
        #plt.plot(hmat[:,0],np.linspace(0,40,20),'k--',label='t=0 position')
        if firsttime==1:
            firsttime=2
            di = np.argwhere(ts_to_do==t)[0][0]
        if t==starting_t:
            l2,= plt.plot(db_plot[:,i],z,'r.',label='inelastic node')
            l3, = plt.plot(db_plot[:,i][~Inelastic_Flag[:,i]],z[~Inelastic_Flag[:,i]],'g.',label='elastic node')
            plt.legend()
            firsttime=1
        else:
            l2.set_xdata(db_plot[:,i])
            l3.remove()
            l3, = plt.plot(db_plot[:,i][~Inelastic_Flag[:,i]],z[~Inelastic_Flag[:,i]],'g.')
        if datelabels=='year':
            plt.title('t=%s' % num2date(t_tmp[i]).strftime('%Y'))
        if datelabels=='monthyear':
            plt.title('t=%s' % num2date(t_tmp[i]).strftime('%b-%Y'))
    
    #            set_size(plt.gcf(), (12, 12))
        plt.savefig('%s/compactionvideo_%s/frame%06d.jpg' % (outputfolder, layer.replace(' ','_'),i),dpi=60,bbox_inches='tight')
    print("")     
    print('')
    print('\t\tStitching frames together using ffmpeg.')
    #cmd='ffmpeg -hide_banner -loglevel warning -r 10 -f image2 -i %*.jpg vid.mp4'
    cmd='ffmpeg -hide_banner -loglevel warning -r 5 -f image2 -i %%*.jpg -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" vid_b_inst_years%sto%s.mp4' % (num2date(starting_t).strftime('%Y'),num2date(ending_t).strftime('%Y'))
    
    print('\t\t\t%s.' % cmd)
    cd=os.getcwd()
    os.chdir('%s/compactionvideo_%s' % (outputfolder, layer.replace(' ','_')))
    subprocess.call(cmd,shell=True)
    
    if os.path.isfile('vid_b_inst_years%sto%s.mp4' % (num2date(starting_t).strftime('%Y'),num2date(ending_t).strftime('%Y'))):
          print('\tVideo seems to have been a success; deleting excess .jpg files.')
          jpg_files_tmp = glob.glob('*.jpg')
          for file in jpg_files_tmp:
              os.remove(file)
          print('\tDone.')
    
    os.chdir(cd)



## Misc functions 
        
from matplotlib.image import imread
from tempfile import NamedTemporaryFile

def get_size(fig, dpi=100):
    with NamedTemporaryFile(suffix='.png') as f:
        fig.savefig(f.name, bbox_inches='tight', dpi=dpi)
        height, width, _channels = imread(f.name).shape
        return width / dpi, height / dpi

def set_size(fig, size, dpi=100, eps=1e-2, give_up=2, min_size_px=10):
    target_width, target_height = size
    set_width, set_height = target_width, target_height # reasonable starting point
    deltas = [] # how far we have
    while True:
        fig.set_size_inches([set_width, set_height])
        actual_width, actual_height = get_size(fig, dpi=dpi)
        set_width *= target_width / actual_width
        set_height *= target_height / actual_height
        deltas.append(abs(actual_width - target_width) + abs(actual_height - target_height))
        if deltas[-1] < eps:
            return True
        if len(deltas) > give_up and sorted(deltas[-give_up:]) == deltas[-give_up:]:
            return False
        if set_width * dpi < min_size_px or set_height * dpi < min_size_px:
            return False

