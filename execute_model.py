#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### This is the main model script. This is where all the "under the hood" operations occur. Do not edit unless you know what you are doing!! The correct call should be: python execute_model.py parameter_file.par [logfileon=True]

#from __future__ import print_function
from model_functions import *
import datetime
import csv
import distutils.util as ut
from datetime import date
from netCDF4 import Dataset
import time
import shutil
from distutils.dir_util import copy_tree
from shutil import copy2
import pandas as pd
import copy
import pygmt

t_total_start = process_time()

gmt=True # If this is true, head outputs can be saved as .netCDF grid files which are tiny. It uses gmt xyz2grd command, hence requires gmt to be installed.

sns.set_context('talk')

if len(sys.argv) != 2:
    if len(sys.argv)==3:
        logfileon= sys.argv[2]
        if '=' in logfileon:
            if logfileon.split('=')[0]=='logfileon':
                if logfileon.split('=')[1]=='True':
                    logfileon=True
                elif logfileon.split('=')[1]=='False':
                    logfileon=False
                else:
                    print('Execute model error; terminal. 3rd input argument is not being understood. Check if it looks correct; it is %s.' % sys.argv[2])
                    sys.exit(1)
            else:
                print('Execute model error; terminal. 3rd input argument is not being understood. Check if it looks correct; it is %s.' % sys.argv[2])
                sys.exit(1)
        else:
            print('Execute model error; terminal. 3rd input argument is not being understood. Check if it looks correct; it is %s.' % sys.argv[2])
            sys.exit(1)


    else:
        print('Execute model error; terminal. Incorrect number of input arguments. Correct usage: python execute_model.py parameter_file.par [logfileon=True]')
        sys.exit(1)
else:
    logfileon=True

if logfileon:
    output_folder='.'
    run_name='.'
    sys.stdout = Logger(output_folder,run_name)
else:
    print('\tOption %s detected. No log file being produced.' % sys.argv[2])

param_filename=sys.argv[1]


print()
print()
print(''.center(80, '*'))
print('  READING PARAMETERS  '.center(80, '*'))
print(''.center(80, '*'))
print()
param_read_start = process_time()
time.sleep(0.5)


print('Reading parameters from file: %s' % param_filename)

# Open the paramfile and parse it
f = open("%s" % param_filename, "r")
paramfilelines = f.readlines()
f.close()
paramfilelines = [line.strip() for line in paramfilelines]
paramfilelines = [x for x in paramfilelines if not x.startswith('#')]
paramfilelines[:] = [x for x in paramfilelines if x]

#print(paramfilelines) # bugfix line; remove later


# Read in parameters

def make_output_folder(outdestination,overwrite):
    # Function which makes the output folder, and copies across the paramfile; I just defined this as a function for neatness.
    if not os.path.isdir(outdestination):
        print("\t\tMaking output directory at %s. This stdout will now be directed to a log in that folder as well as displayed here." % (outdestination))
        os.mkdir(outdestination)
        os.mkdir('%s/figures' % outdestination)
        os.mkdir('%s/s_outputs' % outdestination)
        os.mkdir('%s/head_outputs' % outdestination)
    else:
        if overwrite==False:
            print('\t\tAdmin error: terminal. Output folder %s already exists and overwrite flag not specified.' % outdestination)
            sys.exit(1)
        if overwrite==True:
            print('\t\tOutput folder %s already exists and overwrite flag specified. Do you want to overwrite? (Y/N)' % outdestination)
            check= input()
            check = ut.strtobool(check)
            if check:
                if os.path.isdir('%s_old' % outdestination):
                    shutil.rmtree('%s_old' % outdestination)
                copy_tree(outdestination,'%s_old' % outdestination)
                shutil.rmtree(outdestination)
                os.mkdir(outdestination)
                os.mkdir('%s/figures' % outdestination)
                os.mkdir('%s/s_outputs' % outdestination)
                os.mkdir('%s/head_outputs' % outdestination)
            else:
                print('\t\tNot overwriting. Aborting.')
                sys.exit(1)
    #    OVERWRITE = input("\t\tOutput directory %s already exists. Do you want to overwrite this directory? WARNING: may delete existing data." % (outdestination))
    
    if logfileon:
        shutil.move('logfile.log','%s/logfile.log' % outdestination)
    if os.path.exists(param_filename):
        copy2(param_filename,"%s/paramfile.par" % outdestination)
    else:
        if param_filename.split('/')[-1]=='paramfile.par':
            print('\tAssuming this is a rerun of old run; copying paramfile over.')
            copy2('%s_old/paramfile.par' % outdestination, "%s/paramfile.par" % outdestination)
        else:
            print('\tSomething has gone wrong setting up output directories, ABORT.')
            sys.exit(1)

MODE = read_parameter('mode',str,1,paramfilelines)

def read_parameters_admin(paramfilelines):
    internal_time_delay = read_parameter('internal_time_delay',float,1,paramfilelines)
    overwrite=read_parameter('overwrite',bool,1,paramfilelines)
    run_name = read_parameter('run_name',str,1,paramfilelines)
    output_folder = read_parameter('output_folder',str,1,paramfilelines)
    outdestination="%s/%s" % (output_folder,run_name)
    make_output_folder(outdestination,overwrite)    
    return internal_time_delay,overwrite,run_name,output_folder,outdestination

def read_parameters_noadmin(paramfilelines):
    # Function to read all parameters; defined as a function again for neatness.
          
    save_output_head_timeseries = read_parameter('save_output_head_timeseries',bool,1,paramfilelines)
    save_effective_stress = read_parameter('save_effective_stress',bool,1,paramfilelines)
    save_internal_compaction = read_parameter('save_internal_compaction',bool,1,paramfilelines)
    no_layers = read_parameter('no_layers',int,1,paramfilelines)
    layer_names=read_parameter('layer_names',str,no_layers,paramfilelines)
    if no_layers==1:
        layer_names = np.array([layer_names])
    layer_types=read_parameter('layer_types',str,no_layers,paramfilelines)
    no_aquifers = list(layer_types.values()).count('Aquifer')
    no_aquitards = list(layer_types.values()).count('Aquitard')
    print('\t\tNumber of aquifer layers calculated to be %i.' % no_aquifers)
    print('\t\tNumber of aquitard layers calculated to be %i.' % no_aquitards)
    layer_thickness_types=read_parameter('layer_thickness_types',str,no_layers,paramfilelines)
    layer_thicknesses=read_parameter_layerthickness_multitype('layer_thicknesses',paramfilelines)
    layer_compaction_switch=read_parameter('layer_compaction_switch',bool,no_layers,paramfilelines)
    interbeds_switch=read_parameter('interbeds_switch',bool,list(layer_types.values()).count('Aquifer'),paramfilelines)
    initial_stress_type=read_parameter('initial_stress_type',str,no_layers,paramfilelines)
    no_layers_with_past_stepdrop_initialstress = np.sum(np.isin(list(initial_stress_type.values()),'past_stepdrop'))
    layers_with_past_stepdrop = np.array(list((initial_stress_type.keys())))[np.where(np.isin(list(initial_stress_type.values()),'past_stepdrop'))]
    if no_layers_with_past_stepdrop_initialstress >= 1:
        print('\t\tFound %i layers with past_stepdrop initial stress. Reading additional related parameters initial_stress_paststepdrop_size and initial_stress_paststepdrop_time.' % no_layers_with_past_stepdrop_initialstress)
        initial_stress_paststepdrop_size = read_parameter('initial_stress_paststepdrop_size',float,no_layers_with_past_stepdrop_initialstress,paramfilelines)
        for layer in layers_with_past_stepdrop:
            if layer not in initial_stress_paststepdrop_size:
                print('READ PARAMETERS ERROR: TERMINAL. Layer %s with past_stepdrop does not have initial_stress_paststepdrop_size. Check the input paramfile.' % layer)
                sys.exit()
        initial_stress_paststepdrop_time = read_parameter('initial_stress_paststepdrop_time',float,no_layers_with_past_stepdrop_initialstress,paramfilelines)

        for layer in layers_with_past_stepdrop:
            if layer not in initial_stress_paststepdrop_time:
                print('READ PARAMETERS ERROR: TERMINAL. Layer %s with past_stepdrop does not have initial_stress_paststepdrop_time. Check the input paramfile.' % layer)
                sys.exit()
    else:
        initial_stress_paststepdrop_time=None
        initial_stress_paststepdrop_size=None

    initial_stress_offset=read_parameter('initial_stress_offset',float,no_aquifers,paramfilelines)        
    initial_stress_unit=read_parameter('initial_stress_offset_unit',str,1,paramfilelines)
    # Import interbeds_distributions -- an awkward parameter as its a dictionary of dictionaries!
    interbeds_distributions1=read_parameter('interbeds_distributions',dict,sum(interbeds_switch.values()),paramfilelines)
    interbeds_distributions1=np.array(interbeds_distributions1)
    if np.shape(interbeds_distributions1)[0]==1:
        interbeds_distributions1=interbeds_distributions1[0]
        minidics = [dict([(float(re.split(',|:',interbeds_distributions1[2*i + 1])[2*j]),float( re.split(',|:',interbeds_distributions1[2*i + 1])[2*j+1])) for j in range(int(len( re.split(',|:',interbeds_distributions1[2*i + 1]))/2))]) for i in range(sum(interbeds_switch.values()))]
        interbeds_distributions = dict([(interbeds_distributions1[2*i],minidics[i]) for i in range(sum(interbeds_switch.values()))])
        print('\tinterbeds_distributions=%s' % interbeds_distributions)
    else:
        interbeds_distributions = {}
        for abc in interbeds_distributions1:
            interbeds_distributions[abc[0]] = dict([(float(re.split(':|,',abc[1])[2*i]),float(re.split(':|,',abc[1])[2*i+1])) for i in range(len(re.split(',',abc[1])))])
        print('\tinterbeds_distributions=%s' % interbeds_distributions)
    
    aquitards = [name for name,value in layer_types.items() if value=='Aquitard']
    interbedded_layers= [name for name,value in interbeds_switch.items() if value==True]
    no_layers_containing_clay = len(aquitards) + len(interbedded_layers)
    #no_layers_containing_clay = list(layer_types.values()).count('Aquitard') + sum(list(interbeds_switch.values()))
    print('\t\tNumber of layers containing clay calculated to be %i.' % no_layers_containing_clay)
    
    layers_requiring_solving = interbedded_layers + aquitards
    create_output_head_video= read_parameter('create_output_head_video',bool,no_layers_containing_clay,paramfilelines)
    groundwater_flow_solver_type=read_parameter('groundwater_flow_solver_type',str,len(layers_requiring_solving),paramfilelines)
    if False in [x == 'singlevalue' or x == 'elastic-inelastic' for x in groundwater_flow_solver_type.values()]:
        print("\t\tReading parameters error: terminal. Only groundwater_flow_solver_type of 'singlevalue' or 'elastic-inelastic' currently supported.")
        sys.exit(1)
    overburden_stress_gwflow = read_parameter('overburden_stress_gwflow',bool,1,paramfilelines)
    compaction_solver_compressibility_type = read_parameter('compaction_solver_compressibility_type',str,1,paramfilelines)
    compaction_solver_debug_include_endnodes = read_parameter('compaction_solver_debug_include_endnodes',bool,1,paramfilelines)
        
    clay_Ssk = read_parameter('clay_Ssk',float,sum(value == 'singlevalue' for value in groundwater_flow_solver_type.values()),paramfilelines)
    clay_Sse = read_parameter('clay_Sse',float,sum(groundwater_flow_solver_type[layer]=='elastic-inelastic' or compaction_solver_compressibility_type[layer]=='elastic-inelastic' for layer in layer_names),paramfilelines)
    clay_Ssv = read_parameter('clay_Ssv',float,sum(groundwater_flow_solver_type[layer]=='elastic-inelastic' or compaction_solver_compressibility_type[layer]=='elastic-inelastic' for layer in layer_names),paramfilelines)
    sand_Sse = read_parameter('sand_Sse',float,no_aquifers,paramfilelines)
    
    time_unit = read_parameter('time_unit',str,1,paramfilelines)
    
    #clay_porosity = read_parameter('clay_porosity',float,no_layers_containing_clay,paramfilelines)
    sand_Ssk = read_parameter('sand_Ssk',float,no_aquifers,paramfilelines)
    compressibility_of_water = read_parameter('compressibility_of_water',float,1,paramfilelines)
    rho_w = read_parameter('rho_w',float,1,paramfilelines)
    g = read_parameter('g',float,1,paramfilelines)
    #dt_master = read_parameter('dt_master',float,no_layers_containing_clay,paramfilelines)
    #dz_clays = read_parameter('dz_clays',float,no_layers_containing_clay,paramfilelines)
    n_z = read_parameter('n_z',int,1,paramfilelines)
    CTL_value = read_parameter('CTL_value',float,1,paramfilelines)
    maxdt_ManualCap = read_parameter('maxdt',int,1,paramfilelines)
    vertical_conductivity = read_parameter('vertical_conductivity',float,len(layers_requiring_solving),paramfilelines)
    overburden_stress_compaction = read_parameter('overburden_stress_compaction',bool,1,paramfilelines)
    #overburden_compaction = read_parameter('overburden_compaction',bool,1,paramfilelines)
    save_s = read_parameter('save_s',bool,1,paramfilelines)
    if overburden_stress_gwflow or overburden_stress_compaction: # Only used if we're doing overburden anywhere
        specific_yield = read_parameter('specific_yield',float,1,paramfilelines)
    else:
        specific_yield=None
    return save_output_head_timeseries,save_effective_stress,save_internal_compaction,no_layers,layer_names,layer_types,no_aquifers,no_aquitards,layer_thickness_types,layer_thicknesses,layer_compaction_switch,interbeds_switch,interbeds_distributions,aquitards,interbedded_layers,no_layers_containing_clay,layers_requiring_solving,create_output_head_video,groundwater_flow_solver_type,overburden_stress_gwflow,compaction_solver_compressibility_type,compaction_solver_debug_include_endnodes,clay_Sse,clay_Ssv,clay_Ssk,sand_Sse,time_unit,sand_Ssk,compressibility_of_water,rho_w,g,n_z,CTL_value,maxdt_ManualCap,vertical_conductivity,overburden_stress_compaction,specific_yield,initial_stress_type,initial_stress_offset,initial_stress_unit,save_s,initial_stress_paststepdrop_time,initial_stress_paststepdrop_size

internal_time_delay,overwrite,run_name,output_folder,outdestination = read_parameters_admin(paramfilelines)


if MODE=='resume':
    resume_directory=read_parameter('resume_directory',str,1,paramfilelines)
    if not resume_directory:
        print('\t\tTerminal error: resume_directory not set.')
        sys.exit(1)
    if resume_directory==run_name:
        print('\t\tTerminal error: resume_directory same as output folder.')
        sys.exit(1)

    resume_date=read_parameter('resume_date',str,1,paramfilelines)
    resume_date=dt.strptime(resume_date,'%b-%d-%Y')
    print('\t\tResume date read in as %s' % resume_date)
    
    no_layers = read_parameter('no_layers',int,1,paramfilelines)
    layer_types=read_parameter('layer_types',str,no_layers,paramfilelines)
    no_aquifers = list(layer_types.values()).count('Aquifer')
    resume_head_value=read_parameter('resume_head_value',str,no_aquifers,paramfilelines)
    resume_layer_thicknesses=read_parameter_layerthickness_multitype('layer_thicknesses',paramfilelines)

    
    print('')
    print('*** MODE is RESUME; reading all non-admin parameters from paramfile %s ***' % (resume_directory+'/paramfile.par'))
    copy2('%s/paramfile.par' % resume_directory,"%s/resume_paramfile.par" % outdestination)

    f = open("%s/paramfile.par" % resume_directory, "r")
    paramfilelines = f.readlines()
    f.close()
    paramfilelines = [line.strip() for line in paramfilelines]
    paramfilelines = [x for x in paramfilelines if not x.startswith('#')]
    paramfilelines[:] = [x for x in paramfilelines if x]
    print('')


save_output_head_timeseries,save_effective_stress,save_internal_compaction,no_layers,layer_names,layer_types,no_aquifers,no_aquitards,layer_thickness_types,layer_thicknesses,layer_compaction_switch,interbeds_switch,interbeds_distributions,aquitards,interbedded_layers,no_layers_containing_clay,layers_requiring_solving,create_output_head_video,groundwater_flow_solver_type,overburden_stress_gwflow,compaction_solver_compressibility_type,compaction_solver_debug_include_endnodes,clay_Sse,clay_Ssv,clay_Ssk,sand_Sse,time_unit,sand_Ssk,compressibility_of_water,rho_w,g,n_z,CTL_value,maxdt_ManualCap,vertical_conductivity,overburden_stress_compaction,specific_yield,initial_stress_type,initial_stress_offset,initial_stress_unit,save_s,initial_stress_paststepdrop_time,initial_stress_paststepdrop_size = read_parameters_noadmin(paramfilelines)

if MODE=='resume':
    print('Mode=RESUME, therefore setting initial_stress_type to preset for ALL LAYERS.')
    for layer in layer_names:
        initial_stress_type[layer] = 'preset'
    print('\tInitial_stress_type=%s' % initial_stress_type)


# Check that the layer thicknesses were correctly imported
print()
print('PARAMETER CHECK: checking layer_thickness_types and layer_thicknesses...')

layers_cst_thickness = [k for k,v in layer_thickness_types.items() if v == 'constant']
layers_var_thickness = [k for k,v in layer_thickness_types.items() if v != 'constant']
for layer in layers_cst_thickness:
    if type(layer_thicknesses[layer])==dict:
        print("ERROR, terminal, layer thicknesses for %s is %s. It's a constant thickness layer so shouldn't be a dictionary." % (layer,layer_thicknesses[layer]))
        sys.exit(1)
    else:
        print('\t%s look good.' % layer)
for layer in layers_var_thickness:
    if type(layer_thicknesses[layer])!=dict:
        print("ERROR, terminal, layer thicknesses for %s is %s. It's a variable thickness layer so should be a dictionary." % (layer,layer_thicknesses[layer]))
        sys.exit(1)
    else:
        print('\t%s look good.' % layer)

for layer in layers_var_thickness:
    if layer_thickness_types[layer]=='step_changes':
        tmp = layer_thicknesses[layer]
        pre = ['pre' in s for s in tmp]
        if sum(pre)==1:
            print("\tExactly 1 'pre' entry for %s, looks good." % layer)
        else:
            print("\tERROR:terminal. %s is variable thickness but doesn't have a pre entry. Needs fixing!." % layer)
            sys.exit(1)

# Check maxdt is a multiple of 2

print('Checking maxdt is set to acceptable value.')
if maxdt_ManualCap > 1:
    remainder_tmp = np.log2(maxdt_ManualCap)
    if int(remainder_tmp)!=remainder_tmp:
        print('\tERROR:terminal. maxdt is %i, which is neither <1 nor a multiple of 2. Change maxdt and try again.' % maxdt_ManualCap)
        sys.exit(1)
    else:
        print('\tMaxdt value is fine.')


if len(layers_var_thickness)>=1:
    initial_thicknesses={}
    for layer in layers_var_thickness:
        prekeyname = np.array(list(layer_thicknesses[layer].keys()))[np.where(['pre' in key for key in list(layer_thicknesses[layer].keys())])[0][0]]
        initial_thicknesses[layer] = layer_thicknesses[layer][prekeyname]
    print('\tInitial thicknesses for varying aquifer thicknesses are %s.' % initial_thicknesses)

if MODE=='resume':
    print('\t MODE is RESUME, therefore overriding original layer thicknesses with resume layer thicknesses. Currently, resume does not allow for varying layer thicknesses.')
    layer_thicknesses = resume_layer_thicknesses
    layers_var_thickness=[]
    for layer in layer_names:
        layer_thickness_types[layer]='constant'
    
param_read_stop = process_time()
param_read_time = param_read_start - param_read_stop


#%% Next section is "READING INPUT DATA MODULE"
print()
print()
print(''.center(80, '*'))
print('  READING INPUT DATA  '.center(80, '*'))
print(''.center(80, '*'))
print()
reading_head_start = process_time()
time.sleep(internal_time_delay)



aquifer_layer_names = [name for name, layer_type in layer_types.items() if layer_type=='Aquifer']
compactable_aquifers_names = [name for name, switch in layer_compaction_switch.items() if name in aquifer_layer_names if switch==True]
aquitard_locations = [layer_names.index(name) for name, layer_type in layer_types.items() if layer_type=='Aquitard']
aquifers_above_aquitards = [layer_names[i-1] for i in aquitard_locations]
aquifers_below_aquitards = [layer_names[i+1] for i in aquitard_locations]
all_aquifers_needing_head_data = list(set(aquifers_below_aquitards +aquifers_above_aquitards + compactable_aquifers_names))
dt_headseries={}
print('Preparing to read in head data. Aquifers for which head data is required are: '+', '.join(map(str,all_aquifers_needing_head_data))+'.')
if len(all_aquifers_needing_head_data) >=0:
    try:
        os.mkdir('%s/input_data' % outdestination)
    except FileExistsError:
        pass
    head_data_files=read_parameter('head_data_files',str,len(all_aquifers_needing_head_data),paramfilelines)
    head_data=copy.deepcopy(head_data_files)
    for aquifer in all_aquifers_needing_head_data:
        fileloc=head_data_files[aquifer]
        print('\tFile for %s specified as %s. Looking for file.' % (aquifer,fileloc))
        if os.path.isfile(fileloc):
            print('\t\tFile %s exists. Storing copy in output folder.' % fileloc)
            print('\t\tReading in head time series.')
            copy2(fileloc,'%s/input_data/%s' % (outdestination,fileloc.split('/')[-1]))
            #dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d') # IF you have trouble reading dates in, use these lines!
            try:
                #data=pd.read_csv(fileloc,parse_dates=[0],date_parser=dateparse) # IF you have trouble reading dates in, use these lines!
                data=pd.read_csv(fileloc,parse_dates=[0])

            except Exception:
                print('\t\tReading head data error: terminal. Input file does not seem to be valid csv format. Format expected is two columns, date then measurement. Date should be "dd-MMM-YYYY".')
                sys.exit(1)
            
            try:
                dates=date2num(data.iloc[:,0].values)
            except Exception:
                print("\t\tPandas couldn't parse the date head. Going to try treating the date as a float. If it's not, things may fail from hereon.")
                dates=np.array([float(ting) for ting in data.iloc[:,0].values])
                if time_unit=='years':
                    dates=365*dates
                    
            data=data.iloc[:,1].values
            head_data[aquifer]=np.array([dates,data]).T
            print('\t\tSuccessfully read in. Head data for %s printing now.' % aquifer)
            print(head_data[aquifer])
            dt_headseries[aquifer] = np.diff(head_data[aquifer][:,0])[0]
        else:
            print('\t\tReading head data error: terminal. File %s not found.' % fileloc)
            sys.exit(1)
        
    print('\tReading head data complete.')
else:
    print('\tNo aquifers requiring head data; skipping reading head data.')
starttimes = [np.min(head_data[aquifer][:,0]) for aquifer in all_aquifers_needing_head_data]
print(starttimes)
starttime = np.max(starttimes)
endtimes = [np.max(head_data[aquifer][:,0]) for aquifer in all_aquifers_needing_head_data]
endtime = np.min(endtimes)




print('')
print('CLIPPING HEAD TIMESERIES TO HAVE CONSISTENT START/END DATES ACROSS AQUIFERS.')
print('\tLatest startdate found is %s and earliest end date is %s. These will be used as model start/end times.' % (num2date(starttime).strftime('%d-%b-%Y'),num2date(endtime).strftime('%d-%b-%Y')))
# Clip to have common starting date
print('Clipping input series to model starttime.')
if MODE=='resume':
    starttime = date2num(resume_date)
    print('\tHead startdate is resume date; clipping series accordingly to start at %s.' % resume_date.strftime('%d-%b-%Y'))


for aquifer in all_aquifers_needing_head_data:
    idx_to_keep = ((head_data[aquifer][:,0]>=starttime) & (head_data[aquifer][:,0]<=endtime)) 
    datesnew = head_data[aquifer][:,0][idx_to_keep]
    datanew = head_data[aquifer][:,1][idx_to_keep]
    head_data[aquifer]=np.array([datesnew,datanew]).T
print('Clipping done.')

if MODE=='resume':
    for aquifer in all_aquifers_needing_head_data:
            if resume_head_value[aquifer]=='cst':
                print('\tResume head value=cst. Setting head in aquifer layers to be constant at the value it was in %s.' % resume_date)
                print('Setting head in %s to be %.2f.' % (aquifer,head_data[aquifer][:,1][0]))
                head_data[aquifer][:,1] = head_data[aquifer][:,1][0]
                print(head_data[aquifer])

for aquifer in all_aquifers_needing_head_data:
    with open('%s/input_data/input_time_series_%s.csv' % (outdestination, aquifer.replace(' ','_')), "w+") as myCsv:
        csvWriter = csv.writer(myCsv, delimiter=',')
        csvWriter.writerows(head_data[aquifer])

if overburden_stress_gwflow or overburden_stress_compaction:
    for aquifer in all_aquifers_needing_head_data:
        if aquifer == layer_names[0]:
            print('\t%s is the uppermost aquifer; therefore it is unconfined. Calculating overburden stress from changing water levels in this aquifer.' % aquifer)
            unconfined_aquifer_name = aquifer
            overburden_dates = head_data[aquifer][:,0]
            overburden_data = [specific_yield * rho_w * g * (head_data[aquifer][i,1] - head_data[aquifer][0,1]) for i in range(len(head_data[aquifer][:,1]))]


if overburden_stress_gwflow or overburden_stress_compaction:
    print('Clipping overburden stress.')
    idx_to_keep = ((overburden_dates>=starttime) & (overburden_dates<=endtime)) 
    overburden_dates = overburden_dates[idx_to_keep]
    overburden_data = np.array(overburden_data)[idx_to_keep]
    print('Clipping done.')
    plt.plot_date(overburden_dates,overburden_data)
    plt.savefig('%s/input_data/overburden_stress_series.png' % outdestination)
    with open('%s/input_data/overburden_data.csv' % outdestination, "w+") as myCsv:
        csvWriter = csv.writer(myCsv, delimiter=',')
        csvWriter.writerows([overburden_dates,overburden_data])
    print('\tOverburden stress calculated and saved in input_data.')                    
effective_stress={}


plt.figure(figsize=(18,12))
sns.set_style('darkgrid')
sns.set_context('notebook')
for aquifer in all_aquifers_needing_head_data:
    plt.plot_date(head_data[aquifer][:,0],head_data[aquifer][:,1],label='%s' % aquifer)
plt.ylabel('Head (masl)')
plt.legend()
plt.savefig('%s/input_data/input_head_timeseries.png' % outdestination)
plt.savefig('%s/input_data/input_head_timeseries.pdf' % outdestination)
#plt.savefig('%s/input_data/input_head_timeseries.svg' % outdestination)
plt.close()
sns.set_style('white')

# Check there are no inconsistencies with any aquifers which change thickness
layers_var_thickness = [k for k,v in layer_thickness_types.items() if v != 'constant']
if len(layers_var_thickness)>=1:
    for layer in layers_var_thickness:
        firstyear_thickness = int(np.array(list(layer_thicknesses[layer].keys()))[np.where(['pre' in key for key in list(layer_thicknesses[layer].keys())])[0][0]].split('-')[1])
        if date2num(datetime.datetime(firstyear_thickness,9,1)) < starttime:
            print('\t\tReading head dates error: TERMINAL. The layer %s changes thickness in %i, but the model starts in %s. Please check the layer_thicknesses variable and/or the input head.' % (layer,firstyear_thickness,num2date(starttime).strftime('%Y')))
            sys.exit(1)

# Convert the stress-head units for offsets if needed
print()
if initial_stress_unit=='stress':
    print('Initial stress unit is STRESS. Converting to head.')
    for layer in aquifer_layer_names:
        initial_stress_offset[layer] = initial_stress_offset[layer] / (rho_w*g)
    print('\tDONE.')
    

print()
if len(layers_requiring_solving)>= 0:
    print('Making input clay distribution plot.')
    thicknesses_tmp = []
    for layer in layers_requiring_solving:
        if layer_types[layer]=='Aquifer':
            for key in interbeds_distributions[layer].keys():
                thicknesses_tmp.append(key)
        if layer_types[layer]=='Aquitard':
            if layer_thickness_types[layer]=='constant':
                thicknesses_tmp.append(layer_thicknesses[layer])
            else:
                print('Error, aquitards with varying thickness not (yet) supported.')
                sys.exit(1)
    # Find the smallest difference between two clay layer thicknesses. If that is greater than 1, set the bar width to be 1. Else, set the bar width to be that difference.
    print(thicknesses_tmp)
    if len(thicknesses_tmp)>1:
        diffs=[np.array(thicknesses_tmp) - t for t in thicknesses_tmp]
        smallest_width = np.min(np.abs(np.array(diffs)[np.where(diffs)]))
    else:
        smallest_width=1
    if smallest_width>0.5:
        smallest_width=0.5
    barwidth=smallest_width/len(layers_requiring_solving)
    
    # Make the clay distribution plot
    sns.set_context('poster')
    plt.figure(figsize=(18,12))
    layeri=0
    for layer in layers_requiring_solving:
        if layer_types[layer]=='Aquifer':
            plt.bar(np.array(list(interbeds_distributions[layer].keys()))+layeri*barwidth,list(interbeds_distributions[layer].values()),width=barwidth,label=layer,color='None',edgecolor=sns.color_palette()[layeri],linewidth=5,alpha=0.6)
            layeri+=1
        if layer_types[layer]=='Aquitard':
            plt.bar(layer_thicknesses[layer]+layeri*barwidth,1,width=barwidth,label=layer,color='None',edgecolor=sns.color_palette()[layeri],linewidth=5,alpha=0.6)
            layeri+=1
    plt.legend()
    plt.xticks(np.arange(0,np.max(thicknesses_tmp)+barwidth+1,np.max([1,barwidth])))
    plt.xlabel('Layer thickness (m)')
    plt.ylabel('Number of layers')
    plt.savefig('%s/input_data/clay_distributions.png' % outdestination,bbox_inches='tight')
    plt.close()

reading_head_stop = process_time()
reading_head_time = reading_head_stop - reading_head_start


#%% New section, head solver.
print()
print()
print(''.center(80, '*'))
print('  SOLVING FOR HEAD TIME SERIES IN CLAY LAYERS  '.center(80, '*'))
print(''.center(80, '*'))
print()
solving_head_start = process_time()
time.sleep(internal_time_delay)

print('Hydrostratigraphy:')
print(''.center(60,'-'))
for layer in layer_names:
    text=layer+' (%s)' % layer_types[layer]
    print(text.center(60,' '))
    print(''.center(60,'-'))
print()
print()

print('Head time series to be solved within the following layers: %s' % layers_requiring_solving)

inelastic_flag = {}
inelastic_flag_compaction = {}
Z={}
t_gwflow={}
max_solver_dt = 0 # This variable is a weird one, but it's the largest dt used at any time. I return to this when summing up the results!
dt_head_solver = {}
dz_head_solver = {}

head_series=copy.deepcopy(head_data)

initial_maxstress={}
#if save_output_head_timeseries:
#    os.mkdir('%s/head_outputs' % outdestination)

if len(layers_requiring_solving)>=0:
    groundwater_solution_dates={}
    for layer in layers_requiring_solving:
        print('')
        print('\tBeginning solving process for layer %s.' % layer)
        if layer_types[layer]=='Aquitard':
            print('\t\t%s is an aquitard.' % layer)
            aquitard_position= layer_names.index(layer)
            top_boundary = layer_names[aquitard_position-1]
            bot_boundary = layer_names[aquitard_position+1]
            initial_maxstress[layer]=np.array([])
            print('\t\tHead time series required for overlying layer %s and lower layer %s.' % (top_boundary,bot_boundary))
            if top_boundary in head_data.keys() and bot_boundary in head_data.keys():
                print('\t\t\tHead time series found.')
            
            t_top=head_data[top_boundary][:,0] # This is raw input timeseries for the boundary condition.
            dt_top_tmp = dt_headseries[top_boundary] 
            t_bot=head_data[bot_boundary][:,0] # This is raw input timeseries for the boundary condition.
            dt_bot_tmp = dt_headseries[bot_boundary] 
            dt_aquitard_tmp = np.max([dt_top_tmp,dt_bot_tmp])
            z_tmp = np.linspace(0,layer_thicknesses[layer],n_z)
            dz_tmp = np.diff(z_tmp)[0]
            Z[layer]=z_tmp

            # At this point, calculate the timestep to give a CTL condition just greater than 0.5
            # Assume it will only be the elastic that matters.
            if clay_Sse[layer] < clay_Ssv[layer]:
                print('\t\t\tFinding a good timestep at which to solve...')
                k_tmp = vertical_conductivity[layer]/clay_Sse[layer] 
                dt_optimal = CTL_value * dz_tmp**2 / k_tmp # I must pick a timestep equal to or smaller than this
                print('\t\t\t\tOptimal (or maximum convergent) timestep is %.4f to four decimal places.' % dt_optimal)
                print('\t\t\t\tRaw timestep given by the max of the over/underlying head data is %.2f.' % dt_aquitard_tmp)
                print('\t\t\t\tSelecting timestep to avoid complicated resampling of input data...')
                # I will pick a timestep assuming that I can only refine the input head data or upscale (upscale means sample every n points)
                A_tmp = dt_optimal / dt_aquitard_tmp
                if A_tmp > 2:
                    dt_tmp = 2**np.floor(np.log2(A_tmp)) * dt_aquitard_tmp # Round down to the nearest multiple of two. If not a multiple of two, you might have head time series with prime numbers i.e. cannot add to a different clay.

                elif A_tmp > 1: # due to the else, this means A is between 1 and 2
                    dt_tmp = dt_aquitard_tmp
                elif A_tmp < 1:
                    dt_tmp = 1/(np.ceil(1/A_tmp))
                print('\t\t\t\t\tSelected to use timestep of %.4f.' % dt_tmp)
            else: 
                print('ERROR: terminal. Sse greater than ssv. Code not yet set up to calculate the CTL condition in that cirumstance. EXITING!.')
                sys.exit()

            if dt_tmp > maxdt_ManualCap:
                print('WARNING: selected dt exceeds maxdt input variable. Solving at the smaller, maxdt input instead.')
                dt_tmp = maxdt_ManualCap

            if dt_tmp > max_solver_dt:
                max_solver_dt = dt_tmp
            
            dt_head_solver[layer] = dt_tmp
            dz_head_solver[layer] = dz_tmp

            if not all(t_top == t_bot):
                print('\t\t\tSolving head series error: TERMINAL. Time series in %s and %s aquifers have different dates.' % (top_boundary,bot_boundary))
                print(t_top)
                print(t_bot)
                sys.exit(1)
            else:
                print('\t\t\tTime series found with correct dt and over same timespan.')          

            t_in = np.linspace(min(t_top),min(t_top) + int((max(t_top)-min(t_top))/dt_tmp) * dt_tmp, num = int((max(t_top)-min(t_top))/dt_tmp) + 1)
            f_tmp_top = scipy.interpolate.interp1d(t_top,head_data[top_boundary][:,1])
            top_head_tmp = np.array([t_in,f_tmp_top(t_in)]).T
            f_tmp_bot = scipy.interpolate.interp1d(t_bot,head_data[bot_boundary][:,1])
            bot_head_tmp = np.array([t_in,f_tmp_bot(t_in)]).T

            
            t_gwflow[layer]=t_in
                        
            if overburden_stress_gwflow:
                if layer != unconfined_aquifer_name:
                    print('\t\tPreparing overburden stress.')
                    if len(overburden_dates) != len(t_in):
                        print('\t\t\tInterpolating overburden stress.')
                        f_tmp = scipy.interpolate.interp1d(overburden_dates,overburden_data)
                        overburden_data_tmp = f_tmp(t_in)
                        overburden_dates_tmp = t_in
                    else:
                        overburden_data_tmp = overburden_data
                else: 
                    overburden_data_tmp = np.zeros_like(t_in)
            else:
                overburden_data_tmp = [0]
            
            if initial_stress_type[layer]=='initial_constant':
                preset_initial_maxstress=False
                print('\t\Initial stress type is initial_constant, so a constant stress initial condition will be applied. Constant value is mean of surrounding aquifers.')
                initial_maxstress[layer]=np.array([])
                initial_condition_tmp = ((top_head_tmp[0,1] + bot_head_tmp[0,1]) / 2) * np.ones_like(z_tmp) 
                print('\t\t\tinitial head value is %.2f' % initial_condition_tmp[0])
            elif initial_stress_type[layer]=='initial_equilibrium':
                print("\tInitial stress type is equilibrium_initial, so a linear initial stress condition will be applied.")
                # print('top head is %.2f' % top_head_tmp[0,1] + initial_stress_offset[top_boundary])                                 
                # print('bot head is %.2f' % bot_head_tmp[0,1]          
                # print('m is %.2f' % ((bot_head_tmp[0,1] - top_head_tmp[0,1]) / np.max(z)))                       

                initial_condition_tmp = ((bot_head_tmp[0,1] + initial_stress_offset[bot_boundary] - top_head_tmp[0,1] - initial_stress_offset[top_boundary]) / np.max(z_tmp)) * z_tmp + top_head_tmp[0,1] + initial_stress_offset[top_boundary]
            if MODE=='resume':
                print('\t\tMode is resume. Looking for initial condition in directory %s/head_outputs.' % resume_directory)
                if os.path.isfile("%s/head_outputs/%s_head_data.nc" % (resume_directory,layer.replace(' ','_'))):
                    print('Head found as .nc file. Reading.')
                    Dat = Dataset("%s/head_outputs/%s_head_data.nc" % (resume_directory,layer.replace(' ','_')), "r", format="CF-1.7")
                    time_bc_tmp = Dat.variables['time'][:]
                    head_bc_tmp = Dat.variables['z'][:]
                    idx_bc_tmp = np.argmin(np.abs(date2num(resume_date) - time_bc_tmp))
                    if np.min(np.abs(date2num(resume_date) - time_bc_tmp))>=1:
                        print('\tNote that you are are resuming with the initial head condition from %s, but the specified resume date was %s.' % (time_bc_tmp[idx_bc_tmp],resume_date))
                elif os.path.isfile("%s/head_outputs/%s_head_data.csv" % (resume_directory,layer.replace(' ','_'))):
                    print('Head found as .csv file. Reading.')
                    head_bc_tmp=np.genfromtxt("%s/head_outputs/%s_head_data.csv" % (resume_directory,layer.replace(' ','_')),delimiter=',') 
                    time_bc_tmp1 = np.core.defchararray.rstrip(np.genfromtxt('%s/head_outputs/%s_groundwater_solution_dates.csv' % (resume_directory,layer.replace(' ','_')),dtype=str,delimiter=','))
                    if time_bc_tmp1[0][-1] =='M': # This means it was done on Bletchley in xterm, so %c won't work as it whacks AM or PM on the end.
                        time_bc_tmp = date2num([dt.strptime(string, '%a %d %b %Y %I:%M:%S %p') for string in time_bc_tmp1])
                    else:
                        time_bc_tmp = date2num([dt.strptime(string, '%c') for string in time_bc_tmp1])
                    idx_bc_tmp = np.argmin(np.abs(date2num(resume_date) - time_bc_tmp))

                else:
                    print('\tUnable to find head file as .nc or .csv. Something has gone wrong; aborting.')
                    sys.exit(1)    

                print('\t\tNow looking for effective stress initial condition in directory %s.' % resume_directory)
                if os.path.isfile("%s/%seffective_stress.nc" % (resume_directory,layer.replace(' ','_'))):
                    print('t_eff found as .nc file. Reading.')
                    Dat = Dataset("%s/%seffective_stress.nc" % (resume_directory,layer.replace(' ','_')), "r", format="CF-1.7")
                    time_teff_tmp = Dat.variables['time'][:]
                    teff_bc_tmp = Dat.variables['z'][:] / (rho_w*g) # Get it into units of head
                    idx_teff_bc_tmp = np.argmin(np.abs(date2num(resume_date) - time_teff_tmp))
                    if np.min(np.abs(date2num(resume_date) - time_teff_tmp))>=1:
                        print('\tNote that you are are resuming with the initial t_eff condition from %s, but the specified resume date was %s.' % (time_bc_tmp[idx_bc_tmp],resume_date))
                elif os.path.isfile("%s/%sclayeffective_stress.csv" % (resume_directory,layer.replace(' ','_'))):                   
                    print('t_eff found as .csv file. Reading.')
                    teff_bc_tmp=np.genfromtxt("%s/%seffective_stress.csv" % (resume_directory,layer.replace(' ','_')),delimiter=',') 
                    teff_bc_tmp = teff_bc_tmp/(rho_w*g)
                    time_teff_tmp1 = np.core.defchararray.rstrip(np.genfromtxt('%s/head_outputs/%s_groundwater_solution_dates.csv' % (resume_directory,layer.replace(' ','_')),dtype=str,delimiter=','))
                    if time_bc_tmp1[0][-1] =='M': # This means it was done on Bletchley in xterm, so %c won't work as it whacks AM or PM on the end.
                        time_teff_tmp = date2num([dt.strptime(string, '%a %d %b %Y %I:%M:%S %p') for string in time_teff_tmp1])
                    else:
                        time_teff_tmp = date2num([dt.strptime(string, '%c') for string in time_teff_tmp1])
                    idx_teff_bc_tmp = np.argmin(np.abs(date2num(resume_date) - teff_bc_tmp))

                    print("T_eff resume not yet coded for a csv file - abort.")
                    sys.exit(1)
                else:
                    print('\tUnable to find t_eff file as .nc or .csv. Something has gone wrong; aborting.')
                    sys.exit(1)    


                
                initial_condition_tmp = head_bc_tmp[:,idx_bc_tmp][::-1] # The [::-1] is because aquitard head is saved upside down
                initial_maxstress[layer] = np.max(teff_bc_tmp[:,:idx_teff_bc_tmp+1],axis=1)[::-1]

            
            if groundwater_flow_solver_type[layer] == 'singlevalue':
                hmat_tmp=solve_head_equation_singlevalue(dt_master[layer],t_in,dz_clays[layer],z_tmp,np.vstack((top_head_tmp[::spacing_top,1],top_head_tmp[::spacing_bot,1])),initial_condition_tmp,vertical_conductivity[layer]/(clay_Ssk[layer]+compressibility_of_water))
            elif groundwater_flow_solver_type[layer] == 'elastic-inelastic':
                t1_start = process_time() 
                # hmat_tmp,inelastic_flag_tmp=solve_head_equation_elasticinelastic(dt_master[layer],t_interp_new,dz_clays[layer],z_tmp,np.vstack((h_aquifer_tmp_interpolated[:,1],h_aquifer_tmp_interpolated[:,1])),initial_condition_tmp,vertical_conductivity[layer]/clay_Sse[layer],vertical_conductivity[layer]/clay_Ssv[layer],overburdenstress=overburden_stress_gwflow,overburden_data=1/(rho_w * g) * np.array(overburden_data_tmp))
                hmat_tmp,inelastic_flag_tmp=solve_head_equation_elasticinelastic(dt_tmp,t_in,dz_tmp,z_tmp,np.vstack((top_head_tmp[:,1],bot_head_tmp[:,1])),initial_condition_tmp,vertical_conductivity[layer]/clay_Sse[layer],vertical_conductivity[layer]/clay_Ssv[layer],overburdenstress=overburden_stress_gwflow,overburden_data=1/(rho_w * g) * np.array(overburden_data_tmp),preset_initial_maxstress=preset_initial_maxstress,initial_maxstress=-initial_maxstress[layer]) 
                t1_stop = process_time() 
                print("\t\t\tElapsed time in seconds:",  t1_stop-t1_start)  

            head_series[layer]=hmat_tmp
            inelastic_flag[layer] = inelastic_flag_tmp
            effective_stress[layer] = np.tile(overburden_data_tmp, (np.shape(hmat_tmp)[0],1)) -  rho_w*g*hmat_tmp 

            if save_output_head_timeseries:            
                if np.size(inelastic_flag_tmp) >= 3e6:
                    if gmt:
                        print('\t\t\tInelastic flag gwflow has more than 3 million entries; saving as signed char.')
                        inelastic_flag_tmp.astype(np.byte).tofile('%s/head_outputs/%sinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_')))
                        print('\t\t\t\tConverting to netCDF format. Command is:')
                        #cmd_tmp="gmt xyz2grd %s/head_outputs/%sinelastic_flag_GWFLOW -G%s/head_outputs/%sinelastic_flag_GWFLOW.nb -I%i+n/%i+n -R%.3ft/%.3ft/%.3f/%.3f -ZTLc" % (outdestination, layer.replace(' ','_'),outdestination, layer.replace(' ','_'),len(t_gwflow[layer]),n_z,np.min(t_gwflow[layer]),np.max(t_gwflow[layer]),np.min(Z[layer]),np.max(Z[layer]))
                        #print(cmd_tmp)
                        #subprocess.call(cmd_tmp,shell=True)
                        
                        pygmt.xyz2grd(data='%s/head_outputs/%sinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_')),outgrid='%s/head_outputs/%sinelastic_flag_GWFLOW.nb' % (outdestination, layer.replace(' ','_')),spacing='%i+n/%i+n' % (len(t_gwflow[layer]),n_z),region=[str(np.min(t_gwflow[layer]))+'t',str(np.max(t_gwflow[layer]))+'t',np.min(Z[layer]),np.max(Z[layer])],convention='TLc')
                        
                        os.remove('%s/head_outputs/%sinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_')))
    
                    else:
                        print('\t\t\tInelastic flag gwflow has more than 3 million entries; saving as signed char.')
                        inelastic_flag_tmp.astype(np.byte).tofile('%s/head_outputs/%sinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_')))
    
    
                else:
                    with open('%s/head_outputs/%sinelastic_flag_GWFLOW.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                        csvWriter = csv.writer(myCsv, delimiter=',')
                        csvWriter.writerows(inelastic_flag_tmp)





            #dateslist = [x.strftime('%d-%b-%Y') for x in num2date(t_in)]
            groundwater_solution_dates[layer]=t_in

            if overburden_stress_gwflow:
                effective_stress[layer] = np.tile(overburden_data_tmp, (np.shape(hmat_tmp)[0],1)) -  rho_w*g*hmat_tmp 
            else:
                effective_stress[layer] = np.zeros_like(hmat_tmp) -  rho_w*g*hmat_tmp                             


            if save_effective_stress:
                print('\t\tSaving effective stress and overburden stress outputs.')
                if overburden_stress_gwflow:
                    if len(overburden_data_tmp) * len(z_tmp) >= 1e6:
                        print('\t\t\tOverburden stress has more than 1 million entries; saving as 32 bit floats.')
                        overburden_tmp_tosave = np.tile(overburden_data_tmp, (len(z_tmp),1))
                        overburden_tmp_tosave.astype(np.single).tofile('%s/%s_overburden_stress' % (outdestination, layer.replace(' ','_'))) 
                        if gmt:
                            print('\t\t\t\tConverting to netCDF format. Command is:')
                            cmd_tmp="gmt xyz2grd %s/%s_overburden_stress -G%s/%s_overburden_stress.nc -I%i+n/%i+n -R%.3ft/%.3ft/%.3f/%.3f -ZTLf" % (outdestination, layer.replace(' ','_'),outdestination, layer.replace(' ','_'),len(t_gwflow[layer]),n_z,np.min(t_gwflow[layer]),np.max(t_gwflow[layer]),np.min(Z[layer]),np.max(Z[layer]))
                            
                            print(cmd_tmp)
                            subprocess.call(cmd_tmp,shell=True)
                            os.remove('%s/%s_overburden_stress' % (outdestination, layer.replace(' ','_')))                                                         
                    else:                                
                        with open('%s/%s_overburden_stress.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                            csvWriter = csv.writer(myCsv, delimiter=',')
                            csvWriter.writerows(np.tile(overburden_data_tmp, (len(z_tmp),1)))
                else:
                    with open('%s/%s_overburden_stress.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                        csvWriter = csv.writer(myCsv, delimiter=',')
                        csvWriter.writerows(np.zeros_like(hmat_tmp))

                    
                    
                    
                    


                if np.size(effective_stress[layer]) >= 1e6:
                        print('\t\t\tEffective stress has more than 1 million entries; saving as 32 bit floats.')
                        effective_stress[layer].astype(np.single).tofile('%s/%seffective_stress' % (outdestination, layer.replace(' ','_'))) 
                        if gmt:
                            print('\t\t\t\tConverting to netCDF format. Command is:')
                            cmd_tmp="gmt xyz2grd %s/%seffective_stress -G%s/%seffective_stress.nc -I%i+n/%i+n -R%.3ft/%.3ft/%.2f/%.2f -ZTLf" % (outdestination, layer.replace(' ','_'),outdestination, layer.replace(' ','_'),len(t_gwflow[layer]),n_z,np.min(t_gwflow[layer]),np.max(t_gwflow[layer]),np.min(Z[layer]),np.max(Z[layer]))
                            print(cmd_tmp)
                            subprocess.call(cmd_tmp,shell=True)
                            os.remove('%s/%seffective_stress' % (outdestination, layer.replace(' ','_')))                                                         
                else:
                    with open('%s/%seffective_stress.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                        csvWriter = csv.writer(myCsv, delimiter=',')
                        csvWriter.writerows(effective_stress[layer])






        elif layer_types[layer]=='Aquifer':
            head_series[layer]=dict([('Interconnected matrix',head_data[layer])])
            inelastic_flag[layer]={}
            groundwater_solution_dates[layer]={}
            Z[layer]={}
            t_gwflow[layer]={}
            effective_stress[layer]={}
            initial_maxstress[layer]={}
            dt_head_solver[layer]={}
            dz_head_solver[layer] = {}

            print('\t\t%s is an aquifer.' % layer)
            if interbeds_switch[layer]:
                interbeds_tmp=interbeds_distributions[layer]
                bed_thicknesses_tmp=list(interbeds_tmp.keys())
                print('\t\t%s is an aquifer with interbedded clays. Thicknesses of clays to be solved are %s' % (layer,bed_thicknesses_tmp))
                for thickness in bed_thicknesses_tmp:
                    print('')
                    print('\t\tSolving for thickness %.2f.' % thickness)
                    t_aquifer_tmp=head_data[layer][:,0] # This is raw input timeseries for the boundary condition.
                    dt_aquifer_tmp = dt_headseries[layer] # This is the raw input timestep for the boundary condition.
                    h_aquifer_tmp=head_data[layer][:,1]
                    z_tmp = np.linspace(0,thickness,n_z)
                    dz_tmp = np.diff(z_tmp)[0]
                    dz_head_solver[layer]['%.2f clays' % thickness] = dz_tmp

                    # At this point, calculate the timestep to give a CTL condition just greater than the CTL value.
                    # Assume it will only be the elastic that matters.
                    if clay_Sse[layer] < clay_Ssv[layer]:
                        print('\t\t\tFinding a good timestep at which to solve...')
                        k_tmp = vertical_conductivity[layer]/clay_Sse[layer] 
                        dt_optimal = CTL_value * dz_tmp**2 / k_tmp # I must pick a timestep equal to or smaller than this
                        print('\t\t\t\tOptimal (or maximum convergent) timestep is %.4f to four decimal places.' % dt_optimal)
                        print('\t\t\t\tRaw timestep given by the input head data is %.2f.' % dt_aquifer_tmp)
                        print('\t\t\t\tSelecting timestep to avoid complicated resampling of input data...')
                        # I will pick a timestep assuming that I can only refine the input head data or upscale (upscale means sample every n points)
                        A_tmp = dt_optimal / dt_aquifer_tmp
                        if A_tmp > 2:
                            dt_tmp = 2**np.floor(np.log2(A_tmp)) * dt_aquifer_tmp # Round down to the nearest multiple of two. If not a multiple of two, you might have head time series with prime numbers i.e. cannot add to a different clay.
                        elif A_tmp > 1: # due to the else, this means A is between 1 and 2
                            dt_tmp = dt_aquifer_tmp
                        elif A_tmp < 1:
                            if np.ceil(1/A_tmp) == 77: 
# Why is this here? Well, there are issues with very small A_tmp. To explore this, you can play with the following code snippet, which shows how np.linspace fails due to rounding errors for some: 
#for denom_tmp in np.arange(50,100):
#    print(denom_tmp)
#    dt_tmp = 1./denom_tmp
#    t_interp_new = np.linspace(-6501.,-6501. + int((17453.--6501.)/dt_tmp)*dt_tmp,num=int((17453.--6501.)/dt_tmp)+1)
#    print(t_interp_new[-1])
                                dt_tmp=1./78
                            elif np.ceil(1/A_tmp) >90:
                                dt_tmp =1./90
                                print("WARNING: timestep is tiny, such that rounding errors are likely. Continuing but setting dt_tmp to %.10f as the smallest value which doesn't result in rounding errors." % dt_tmp)
                            else:
                                dt_tmp = 1/(np.ceil(1/A_tmp))
                        print('\t\t\t\t\tSelected to use timestep of %.4f.' % dt_tmp)
                    else: 
                        print('ERROR: terminal. Sse greater than ssv. Code not yet set up to calculate the CTL condition in that cirumstance. EXITING!.')
                        sys.exit()
                                    
                    if dt_tmp > maxdt_ManualCap:
                        print('WARNING: selected dt exceeds maxdt input variable. Solving at the smaller, maxdt input instead.')
                        dt_tmp = maxdt_ManualCap


                    if dt_tmp > max_solver_dt:
                        max_solver_dt = dt_tmp                        
    
                    # B_tmp =  (max(t_aquifer_tmp)-min(t_aquifer_tmp))/ dt_tmp
                    # if not B_tmp.is_integer():
                    #     print("\t\tIF THIS ISN'T A WHOLE NUMBER, written as a float, YOU MAY BE IN TROUBLE! WARNING!" )
                    #     print((max(t_aquifer_tmp)-min(t_aquifer_tmp))/ dt_tmp)
                    t_interp_new = np.linspace(min(t_aquifer_tmp),min(t_aquifer_tmp) + int((max(t_aquifer_tmp)-min(t_aquifer_tmp))/ dt_tmp) * dt_tmp ,num=int((max(t_aquifer_tmp)-min(t_aquifer_tmp))/ dt_tmp)+1)
                    f_tmp = scipy.interpolate.interp1d(t_aquifer_tmp,h_aquifer_tmp) # linear interpolation of input head
                    h_aquifer_tmp_interpolated = np.array([t_interp_new,f_tmp(t_interp_new)]).T

                    if initial_stress_type[layer]=='initial_equilibrium' or initial_stress_type[layer]=='initial_constant':
                        print('\t\tHead initial condition is %s, so a constant head initial condition will be applied. Offset is %.2f.' % (initial_stress_type[layer],initial_stress_offset[layer]))
                        preset_initial_maxstress=False
                        initial_condition_tmp = h_aquifer_tmp[0] * np.ones_like(z_tmp) + initial_stress_offset[layer]
                        initial_maxstress[layer]['%.2f clays' % thickness]=np.array([])
                        print('\t\tinitial head value is %.2f' % initial_condition_tmp[0])
                    elif initial_stress_type[layer]=='past_stepdrop':
                        print('\t\tInitial stress condition is past_stepdrop; calculating initial stress accordingly.')
                        preset_initial_maxstress=False
                        initial_condition_tmp = hoffman_subsidence(initial_stress_paststepdrop_time[layer],thickness,h_aquifer_tmp[0] + initial_stress_paststepdrop_size[layer],-initial_stress_paststepdrop_size[layer],clay_Ssv[layer],vertical_conductivity[layer],n_z=n_z)
                        initial_maxstress[layer]['%.2f clays' % thickness]=np.array([])
                        print('Initial condition precalculated to be ',initial_condition_tmp)
                    if MODE=='resume':
                        preset_initial_maxstress=True
                        print('\t\tMode is resume. Looking for initial condition in directory %s/head_outputs.' % resume_directory)
                        if os.path.isfile("%s/head_outputs/%s_%sclay_head_data.nc" % (resume_directory,layer.replace(' ','_'),'%.2f' % thickness)):
                            print('Head found as .nc file. Reading.')
                            Dat = Dataset("%s/head_outputs/%s_%sclay_head_data.nc" % (resume_directory,layer.replace(' ','_'),'%.2f' % thickness), "r", format="CF-1.7")
                            time_bc_tmp = Dat.variables['time'][:]
                            head_bc_tmp = Dat.variables['z'][:]
                            idx_bc_tmp = np.argmin(np.abs(date2num(resume_date) - time_bc_tmp))
                            if np.min(np.abs(date2num(resume_date) - time_bc_tmp))>=1:
                                print('\tNote that you are are resuming with the initial head condition from %s, but the specified resume date was %s.' % (time_bc_tmp[idx_bc_tmp],resume_date))
                        elif os.path.isfile("%s/head_outputs/%s_%sclay_head_data.csv" % (resume_directory,layer.replace(' ','_'),'%.2f' % thickness)):
                            print('Head found as .csv file. Reading.')
                            head_bc_tmp=np.genfromtxt("%s/head_outputs/%s_%sclay_head_data.csv" % (resume_directory,layer.replace(' ','_'),'%.2f' % thickness),delimiter=',') 
                            time_bc_tmp1 = np.core.defchararray.rstrip(np.genfromtxt('%s/head_outputs/%s_groundwater_solution_dates.csv' % (resume_directory,layer.replace(' ','_')),dtype=str,delimiter=','))
                            if time_bc_tmp1[0][-1] =='M': # This means it was done on Bletchley in xterm, so %c won't work as it whacks AM or PM on the end.
                                time_bc_tmp = date2num([dt.strptime(string, '%a %d %b %Y %I:%M:%S %p') for string in time_bc_tmp1])
                            else:
                                time_bc_tmp = date2num([dt.strptime(string, '%c') for string in time_bc_tmp1])
                            idx_bc_tmp = np.argmin(np.abs(date2num(resume_date) - time_bc_tmp))

                        else:
                            print('\tUnable to find head file as .nc or .csv. Something has gone wrong; aborting.')
                            sys.exit(1)    

                        print('\t\tNow looking for effective stress initial condition in directory %s.' % resume_directory)
                        if os.path.isfile("%s/%s_%sclayeffective_stress.nc" % (resume_directory,layer.replace(' ','_'),'%.2f' % thickness)):
                            print('t_eff found as .nc file. Reading.')
                            Dat = Dataset("%s/%s_%sclayeffective_stress.nc" % (resume_directory,layer.replace(' ','_'),'%.2f' % thickness), "r", format="CF-1.7")
                            time_teff_tmp = Dat.variables['time'][:]
                            teff_bc_tmp = Dat.variables['z'][:] / (rho_w*g) # Get it into units of head
                            idx_teff_bc_tmp = np.argmin(np.abs(date2num(resume_date) - time_teff_tmp))
                            if np.min(np.abs(date2num(resume_date) - time_teff_tmp))>=1:
                                print('\tNote that you are are resuming with the initial t_eff condition from %s, but the specified resume date was %s.' % (time_bc_tmp[idx_bc_tmp],resume_date))
                        elif os.path.isfile("%s/%s_%sclayeffective_stress.csv" % (resume_directory,layer.replace(' ','_'),'%.1f' % thickness)):
                            print('t_eff found as .csv file. Reading.')
                            teff_bc_tmp=np.genfromtxt("%s/%s_%sclayeffective_stress.csv" % (resume_directory,layer.replace(' ','_'),'%.1f' % thickness),delimiter=',') 
                            teff_bc_tmp = teff_bc_tmp/(rho_w*g)
                            time_teff_tmp1 = np.core.defchararray.rstrip(np.genfromtxt('%s/head_outputs/%s_groundwater_solution_dates.csv' % (resume_directory,layer.replace(' ','_')),dtype=str,delimiter=','))
                            if time_teff_tmp1[0][-1] =='M': # This means it was done on Bletchley in xterm, so %c won't work as it whacks AM or PM on the end.
                                time_teff_tmp = date2num([dt.strptime(string, '%a %d %b %Y %I:%M:%S %p') for string in time_teff_tmp1])
                            else:
                                time_teff_tmp = date2num([dt.strptime(string, '%c') for string in time_teff_tmp1])
                            idx_teff_bc_tmp = np.argmin(np.abs(date2num(resume_date) - teff_bc_tmp))
                        else:
                            print('\tUnable to find t_eff file as .nc or .csv. Something has gone wrong; aborting.')
                            sys.exit(1)    


                        
                        initial_condition_tmp = head_bc_tmp[:,idx_bc_tmp]
                        initial_maxstress[layer]['%.2f clays' % thickness] = np.max(teff_bc_tmp[:,:idx_teff_bc_tmp+1],axis=1)
                        
                    if overburden_stress_gwflow:
                        if layer != unconfined_aquifer_name:
                            if len(overburden_dates) != len(t_interp_new):
                                print('\t\t\tInterpolating overburden stress.')
                                f_tmp = scipy.interpolate.interp1d(overburden_dates,overburden_data)
                                overburden_data_tmp = f_tmp(t_interp_new)
                                overburden_dates_tmp = t_interp_new
                            else:
                                overburden_data_tmp = overburden_data
                        else:
                            print('\t\t\tThis is the unconfined aquifer; overburden still being included.')
                            if len(overburden_dates) != len(t_interp_new):
                                print('\t\t\tInterpolating overburden stress.')
                                f_tmp = scipy.interpolate.interp1d(overburden_dates,overburden_data)
                                overburden_data_tmp = f_tmp(t_interp_new)
                                overburden_dates_tmp = t_interp_new
                            else:
                                overburden_data_tmp = overburden_data
   
                            #overburden_data_tmp = np.zeros_like(t_interp_new)
                    else:
                        overburden_data_tmp=[0]

                    t1_start = process_time() 
                    hmat_tmp,inelastic_flag_tmp=solve_head_equation_elasticinelastic(dt_tmp,t_interp_new,dz_tmp,z_tmp,np.vstack((h_aquifer_tmp_interpolated[:,1],h_aquifer_tmp_interpolated[:,1])),initial_condition_tmp,vertical_conductivity[layer]/clay_Sse[layer],vertical_conductivity[layer]/clay_Ssv[layer],overburdenstress=overburden_stress_gwflow,overburden_data=1/(rho_w * g) * np.array(overburden_data_tmp),preset_initial_maxstress=preset_initial_maxstress,initial_maxstress=-initial_maxstress[layer]['%.2f clays' % thickness])
                    t1_stop = process_time() 
                    print("\t\t\tElapsed time in seconds:",  t1_stop-t1_start)  
                    head_series[layer]['%.2f clays' % thickness]=hmat_tmp
                    inelastic_flag[layer]['%.2f clays' % thickness] = inelastic_flag_tmp
                    t_gwflow[layer]['%.2f clays' % thickness] = t_interp_new
                    Z[layer]['%.2f clays' % thickness]=z_tmp
                    dt_head_solver[layer]['%.2f clays' % thickness] = dt_tmp

                    if save_output_head_timeseries:
                        if np.size(inelastic_flag_tmp) >= 3e6:
                            if gmt:
                                print('\t\t\tInelastic flag gwflow has more than 3 million entries; saving as signed char.')
                                inelastic_flag_tmp.astype(np.byte).tofile('%s/head_outputs/%s_%sclayinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_'),'%.2f' % thickness))
                                #print('\t\t\t\tConverting to netCDF format. Command is:')
                                # cmd_tmp="gmt xyz2grd %s/head_outputs/%s_%sclayinelastic_flag_GWFLOW -G%s/head_outputs/%s_%sclayinelastic_flag_GWFLOW.nb -I%i+n/%i+n -R%.3ft/%.3ft/%.3f/%.3f -ZTLc" % (outdestination, layer.replace(' ','_'),'%.2f' % thickness,outdestination, layer.replace(' ','_'),'%.2f' % thickness,len(t_gwflow[layer]['%.2f clays' % thickness]),n_z,np.min(t_gwflow[layer]['%.2f clays' % thickness]),np.max(t_gwflow[layer]['%.2f clays' % thickness]),np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness]))
                                # print(cmd_tmp)
                                # subprocess.call(cmd_tmp,shell=True)
                                
                                pygmt.xyz2grd(data='%s/head_outputs/%s_%sclayinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_'),'%.2f' % thickness),outgrid='%s/head_outputs/%s_%sclayinelastic_flag_GWFLOW.nb' % (outdestination, layer.replace(' ','_'),'%.2f' % thickness),spacing='%i+n/%i+n' % (len(t_gwflow[layer]['%.2f clays' % thickness]),n_z),region=[str(np.min(t_gwflow[layer]['%.2f clays' % thickness]))+'t',str(np.max(t_gwflow[layer]['%.2f clays' % thickness]))+'t',np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness])],convention='TLc',D='+xtime+yz+vflag') # UNRESOLVED: for certain file names, this fails with a bizarre 'invalid characters' error. I investigated and couldn't figure it out at all. The D option doesn't work weirdly enough, so the variable name is 'W.nb'... peculiar.

                                
                                os.remove('%s/head_outputs/%s_%sclayinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_'),'%.2f' % thickness))

                            else:
                                print('\t\t\tInelastic flag gwflow has more than 3 million entries; saving as signed char.')
                                inelastic_flag_tmp.astype(np.byte).tofile('%s/head_outputs/%s_%sclayinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_'),'%.2f' % thickness))


                        else:
                            with open('%s/head_outputs/%s_%sclayinelastic_flag_GWFLOW.csv' % (outdestination, layer.replace(' ','_'),'%.2f' % thickness), "w+") as myCsv:
                                csvWriter = csv.writer(myCsv, delimiter=',')
                                csvWriter.writerows(inelastic_flag_tmp)

                    #dateslist = [x.strftime('%d-%b-%Y') for x in num2date(t_interp_new)]
                    groundwater_solution_dates[layer]['%.2f clays' % thickness]=t_interp_new
                    
                    if overburden_stress_gwflow:
                        effective_stress[layer]['%.2f clays' % thickness] = np.tile(overburden_data_tmp, (np.shape(hmat_tmp)[0],1)) -  rho_w*g*hmat_tmp 
                    else:
                        effective_stress[layer]['%.2f clays' % thickness] = np.zeros_like(hmat_tmp) -  rho_w*g*hmat_tmp                             

                    if save_effective_stress:
                        print('\t\tSaving effective stress and overburden stress outputs.')
                        if np.size(effective_stress[layer]['%.2f clays' % thickness]) >= 1e6:
                            print('\t\t\tEffective stress has more than 1 million entries; saving as 32 bit floats.')
                            effective_stress[layer]['%.2f clays' % thickness].astype(np.single).tofile('%s/%s_%sclayeffective_stress' % (outdestination, layer.replace(' ','_'), thickness)) 
                            if gmt:
                                print('\t\t\t\tConverting to netCDF format. Command is:')
                                cmd_tmp="gmt xyz2grd %s/%s_%sclayeffective_stress -G%s/%s_%sclayeffective_stress.nc -I%i+n/%i+n -R%.3ft/%.3ft/%.2f/%.2f -ZTLf" % (outdestination, layer.replace(' ','_'),thickness,outdestination, layer.replace(' ','_'),'%.2f' % thickness,len(t_gwflow[layer]['%.2f clays' % thickness]),n_z,np.min(t_gwflow[layer]['%.2f clays' % thickness]),np.max(t_gwflow[layer]['%.2f clays' % thickness]),np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness]))
                        
                                print(cmd_tmp)
                                subprocess.call(cmd_tmp,shell=True)
                                os.remove('%s/%s_%sclayeffective_stress' % (outdestination, layer.replace(' ','_'),thickness))                                                         
                        else:
                            with open('%s/%s_%sclayeffective_stress.csv' % (outdestination, layer.replace(' ','_'),thickness), "w+") as myCsv:
                                csvWriter = csv.writer(myCsv, delimiter=',')
                                csvWriter.writerows(effective_stress[layer]['%.2f clays' % thickness])

                        if overburden_stress_gwflow:
                            if len(overburden_data_tmp) * len(z_tmp) >= 1e6:
                                print('\t\t\tOverburden stress has more than 1 million entries; saving as 32 bit floats.')
                                overburden_tmp_tosave = np.tile(overburden_data_tmp, (len(z_tmp),1))
                                overburden_tmp_tosave.astype(np.single).tofile('%s/%s_%sclay_overburden_stress' % (outdestination, layer.replace(' ','_'), thickness)) 
                                if gmt:
                                    print('\t\t\t\tConverting to netCDF format. Command is:')
                                    cmd_tmp="gmt xyz2grd %s/%s_%sclay_overburden_stress -G%s/%s_%sclay_overburden_stress.nc -I%i+n/%i+n -R%.3ft/%.3ft/%.3f/%.3f -ZTLf" % (outdestination, layer.replace(' ','_'),thickness,outdestination, layer.replace(' ','_'),'%.2f' % thickness,len(t_gwflow[layer]['%.2f clays'] % thickness),n_z,np.min(t_gwflow[layer]['%.2f clays' % thickness]),np.max(t_gwflow[layer]['%.2f clays' % thickness]),np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness]))
                        
                                    print(cmd_tmp)
                                    subprocess.call(cmd_tmp,shell=True)
                                    os.remove('%s/%s_%sclay_overburden_stress' % (outdestination, layer.replace(' ','_'),thickness))                                                         
                            else:

                            
                                with open('%s/%s_%sclay_overburden_stress.csv' % (outdestination, layer.replace(' ','_'),thickness), "w+") as myCsv:
                                    csvWriter = csv.writer(myCsv, delimiter=',')
                                    csvWriter.writerows(np.tile(overburden_data_tmp, (len(z_tmp),1)))

                        
            else:
                print('\t%s is an aquifer with no interbedded clays. No solution required.')

else:
    print('No layers require head time series solutions; skipping solving for head time series in clay layers.')
    
solving_head_stop = process_time()
solving_head_time = solving_head_stop - solving_head_start
#%% New section, saving head outputs.
print()
print()
print(''.center(80, '*'))
print('  SAVING HEAD TIMESERIES OUTPUTS  '.center(80, '*'))
print(''.center(80, '*'))
print()
saving_head_start = process_time()
time.sleep(internal_time_delay)


if save_output_head_timeseries:
    print('save_output_head_timeseries = True. Saving head timeseries for all layers.')
    for layer in layers_requiring_solving:
        print('\tSaving head timeseries for %s.' % layer)
        if layer_types[layer]=='Aquifer':
            dates_str = [x.strftime('%d-%b-%Y') for x in num2date(head_data[layer][:,0])]
            np.savetxt('%s/head_outputs/%s_head_data.csv' % (outdestination, layer.replace(' ','_')),np.column_stack((dates_str,head_data[layer][:,1])),fmt="%s")
            if interbeds_switch[layer]:
                interbeds_tmp=interbeds_distributions[layer]
                for thickness in list(interbeds_tmp.keys()):
                    if np.size(head_series[layer]['%.2f clays' % thickness]) >= 1e6:
                        if gmt:
                            print('\t\t\tHead has more than 1 million entries; saving as 32 bit floats.')
                            head_series[layer]['%.2f clays' % thickness].astype(np.single).tofile('%s/head_outputs/%s_%sclay_head_data' % (outdestination, layer.replace(' ','_'),thickness))
                            print('\t\t\t\tConverting to netCDF format. Command is:')
                            print()
                            # cmd_tmp="gmt xyz2grd %s/head_outputs/%s_%sclay_head_data -G%s/head_outputs/%s_%sclay_head_data.nc -I%i+n/%i+n -R%.3ft/%.3ft/%.3f/%.3f -ZTLf" % (outdestination, layer.replace(' ','_'),thickness,outdestination, layer.replace(' ','_'),'%.2f' % thickness,len(t_gwflow[layer]['%.2f clays' % thickness]),n_z,np.min(t_gwflow[layer]['%.2f clays' % thickness]),np.max(t_gwflow[layer]['%.2f clays' % thickness]),np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness]))
                            # print(cmd_tmp)
                            # subprocess.call(cmd_tmp,shell=True)
                            pygmt.xyz2grd(data='%s/head_outputs/%s_%sclay_head_data' % (outdestination, layer.replace(' ','_'),thickness),outgrid='%s/head_outputs/%s_%sclay_head_data.nc' % (outdestination, layer.replace(' ','_'),'%.2f' % thickness),spacing='%i+n/%i+n' % (len(t_gwflow[layer]['%.2f clays' % thickness]),n_z),region=[str(np.min(t_gwflow[layer]['%.2f clays' % thickness]))+'t',str(np.max(t_gwflow[layer]['%.2f clays' % thickness]))+'t',np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness])],convention='TLf')

                            
                            
                            os.remove('%s/head_outputs/%s_%sclay_head_data' % (outdestination, layer.replace(' ','_'),thickness))
                        else:
                            print('\t\t\tHead has more than 1 million entries; saving as 16 bit floats.')
                            head_series[layer]['%.2f clays' % thickness].astype(np.half).tofile('%s/head_outputs/%s_%sclay_head_data' % (outdestination, layer.replace(' ','_'),thickness))

                    else:
                        with open('%s/head_outputs/%s_%sclay_head_data.csv' % (outdestination, layer.replace(' ','_'),'%.2f' % thickness), "w+") as myCsv:
                            csvWriter = csv.writer(myCsv, delimiter=',')
                            csvWriter.writerows(head_series[layer]['%.2f clays' % thickness])
            
                with open('%s/head_outputs/%s_%sclay_groundwater_solution_dates.csv' % (outdestination, layer.replace(' ','_'),'%.2f' % thickness), 'w') as myfile:
                    wr = csv.writer(myfile)
                    #res = list(groundwater_solution_dates[layer].keys())[0] 
                    wr.writerow([x.strftime('%c') for x in num2date(groundwater_solution_dates[layer]['%.2d clays' % thickness])])

            
        if layer_types[layer]=='Aquitard':
            if np.size(head_series[layer]) >= 1e6:
                if gmt:
                    print('\t\t\tHead has more than 1 million entries; saving as 32 bit floats.')
                    head_series[layer].astype(np.single).tofile('%s/head_outputs/%s_head_data' % (outdestination, layer.replace(' ','_')))
                    #print('\t\t\t\tConverting to netCDF format. Command is:')
                    # cmd_tmp="gmt xyz2grd %s/head_outputs/%s_head_data -G%s/head_outputs/%s_head_data.nc -I%i+n/%i+n -R%.3ft/%.3ft/%.3f/%.3f -ZTLf" % (outdestination, layer.replace(' ','_'),outdestination, layer.replace(' ','_'),len(t_gwflow[layer]),n_z,np.min(t_gwflow[layer]),np.max(t_gwflow[layer]),np.min(Z[layer]),np.max(Z[layer]))
                    # print(cmd_tmp)
                    # subprocess.call(cmd_tmp,shell=True)
                    pygmt.xyz2grd(data='%s/head_outputs/%s_head_data' % (outdestination, layer.replace(' ','_')),outgrid='%s/head_outputs/%s_head_data.nc' % (outdestination, layer.replace(' ','_')),spacing='%i+n/%i+n' % (len(t_gwflow[layer]),n_z),region=[str(np.min(t_gwflow[layer]))+'t',str(np.max(t_gwflow[layer]))+'t',np.min(Z[layer]),np.max(Z[layer])],convention='TLf')

                    os.remove('%s/head_outputs/%s_head_data' % (outdestination, layer.replace(' ','_')))
                else:
                    print('\t\t\tHead has more than 1 million entries; saving as 16 bit floats.')
                    head_series[layer].astype(np.half).tofile('%s/head_outputs/%s_head_data' % (outdestination, layer.replace(' ','_')))

            else:
                with open('%s/head_outputs/%s_head_data.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                    csvWriter = csv.writer(myCsv, delimiter=',')
                    csvWriter.writerows(head_series[layer])
            with open('%s/head_outputs/%s_groundwater_solution_dates.csv' % (outdestination, layer.replace(' ','_')), 'w') as myfile:
                wr = csv.writer(myfile)
                wr.writerow([x.strftime('%c') for x in num2date(groundwater_solution_dates[layer])])

for layer in layers_requiring_solving:
    if create_output_head_video[layer]:
        print('create_output_head_video = True. Creating head timeseries video for specified layers. Note: requires ffmpeg installed.')
        if layer_types[layer]=='Aquitard':
                print('\tCreating video for %s.' % layer)
                hmat_tmp = head_series[layer]
                inelastic_flag_vid = inelastic_flag[layer]
                inelastic_flag_vid = inelastic_flag_vid==1            
                dates_str = [x.strftime('%d-%b-%Y') for x in num2date(groundwater_solution_dates[layer])]
                if max_solver_dt < 100: # If the timesteps would be less than 100, make the video timestep longer to be almost 100
                    viddelt = int(100/max_solver_dt) * max_solver_dt
                else:
                    viddelt = max_solver_dt
                print('\t\tUsing dt = %i.' % viddelt)
                create_head_video_elasticinelastic(hmat_tmp,Z[layer],inelastic_flag_vid,dates_str,outdestination+'/figures',layer,delt=viddelt)

        if layer_types[layer]=='Aquifer':
            print('\tCreating video for %s.' % layer)
            if max_solver_dt < 100: # If the timesteps would be less than 100, make the video timestep longer to be almost 100
                viddelt = int(100/max_solver_dt) * max_solver_dt
            else:
                viddelt = max_solver_dt
            print('\tUsing dt = %i.' % viddelt)
            interbeds_tmp=interbeds_distributions[layer]
            bed_thicknesses_tmp=list(interbeds_tmp.keys())
            print('\t\t%s is an aquifer with interbedded clays. Thicknesses of clays to make videos are %s' % (layer,bed_thicknesses_tmp))
            for thickness in bed_thicknesses_tmp:
                print('\t\tCreating video for %s_%.2f clays' % (layer,thickness))
                hmat_tmp = head_series[layer]['%.2f clays' % thickness]
                inelastic_flag_vid = inelastic_flag[layer]['%.2f clays' % thickness]
                inelastic_flag_vid = inelastic_flag_vid==1            
                dates_str = [x.strftime('%d-%b-%Y') for x in num2date(groundwater_solution_dates[layer]['%.2f clays' % thickness])]
                create_head_video_elasticinelastic(hmat_tmp,Z[layer]['%.2f clays' % thickness],inelastic_flag_vid,dates_str,outdestination+'/figures','%s_%.2f_clays' % (layer, thickness),delt=viddelt)

saving_head_stop = process_time()
saving_head_time = saving_head_stop - saving_head_start

#%% New section, compaction solver.
print()
print()
print(''.center(80, '*'))
print('  SOLVING COMPACTION EQUATION  '.center(80, '*'))
print(''.center(80, '*'))
print()
solving_compaction_start = process_time()
time.sleep(internal_time_delay)

#deformation_series={}
#deformation_series_elastic={}
#deformation_series_inelastic={}
#deformation_series_sand={}
deformation={}
db={}
deformation_OUTPUT={}
compacting_layers = [name for name,value in layer_compaction_switch.items() if value==True]


for layer in layer_names:
    if layer_types[layer]=='Aquifer':
        if layer_compaction_switch[layer]:
            print()
            print('%s is an Aquifer. Solving for layer compaction.' % layer)
            deformation[layer]={}
            db[layer]={}
            inelastic_flag_compaction[layer]={}
            if layer_thickness_types[layer]=='constant':
                layer_sand_thickness_tmp = layer_thicknesses[layer] - np.sum([list(interbeds_distributions[layer].keys())[i] * list(interbeds_distributions[layer].values())[i] for i in range(len(interbeds_distributions[layer]))])
                print('\tTotal sand thickness in aquifer is %.2f m.' % layer_sand_thickness_tmp)
            elif layer_thickness_types[layer]=='step_changes':
                layer_sand_thickness_tmp = initial_thicknesses[layer] - np.sum([list(interbeds_distributions[layer].keys())[i] * list(interbeds_distributions[layer].values())[i] for i in range(len(interbeds_distributions[layer]))])
                print('\tInitial total sand thickness in aquifer is %.2f m.' % layer_sand_thickness_tmp)
            deformation[layer]['Interconnected matrix']=[layer_sand_thickness_tmp*(sand_Sse[layer]-compressibility_of_water)*(head_series[layer]['Interconnected matrix'][i,1] - head_series[layer]['Interconnected matrix'][0,1]) for i in range(len(head_series[layer]['Interconnected matrix'][:,1]))]
            
            # Now find the timestep at which to solve
            dt_sand_tmp=np.diff(head_data[layer][:,0])[0]
            if dt_sand_tmp > max_solver_dt:
                max_solver_dt = dt_sand_tmp
            print('\tCalculating deformation for layer %s. dt is %.2f.' % (layer, max_solver_dt))
            dt_max_tmp = max_solver_dt
            t_total_tmp = np.linspace(np.min(head_data[layer][:,0]),np.min(head_data[layer][:,0]) + int((np.max(head_data[layer][:,0])-np.min(head_data[layer][:,0]))/dt_max_tmp) * dt_max_tmp,num=int((np.max(head_data[layer][:,0])-np.min(head_data[layer][:,0]))/dt_max_tmp)+1) # this is the master t which will apply for all the solved layers. Note that, for certain dt_max_tmp which are larger than 1, the finish value may be smaller than the greatest input time.
            
            interbeds_tmp=interbeds_distributions[layer]
            bed_thicknesses_tmp=list(interbeds_tmp.keys())
            print('\t\t%s is an aquifer with interbedded clays. Thicknesses of clays to solve compaction are %s' % (layer,bed_thicknesses_tmp))


            for thickness in bed_thicknesses_tmp:
                print('\t\t\tSolving for thickness %.2f.' % thickness)
                
                if overburden_stress_compaction:
                    unconfined_tmp = unconfined_aquifer_name==layer
                    if len(overburden_dates) != len(groundwater_solution_dates[layer]['%.2f clays' % thickness]):
                        print('\t\t\tOverburden series is %i long whereas head series is %i long. Interpolating overburden stress.' % (len(overburden_dates),len(groundwater_solution_dates[layer]['%.2f clays' % thickness])))
                        f_tmp = scipy.interpolate.interp1d(overburden_dates,overburden_data)
                        overburden_data_tmp = f_tmp(groundwater_solution_dates[layer]['%.2f clays' % thickness])
                        overburden_dates_tmp = groundwater_solution_dates[layer]['%.2f clays' % thickness]
                    else:
                        overburden_data_tmp = overburden_data
                    
                    deformation[layer]['total_%.2f clays' % thickness],inelastic_flag_compaction[layer]['elastic_%.2f clays' % thickness]=subsidence_solver_aquitard_elasticinelastic(head_series[layer]['%.2f clays' % thickness][:,np.isin(t_gwflow[layer]['%.2f clays' % thickness],t_total_tmp)],(clay_Sse[layer]-compressibility_of_water),(clay_Ssv[layer]-compressibility_of_water),np.max(Z[layer]['%.2f clays' % thickness]),n_z,unconfined=unconfined_tmp,overburden=overburden_stress_compaction,overburden_data=1/(rho_w * g) * np.array(overburden_data_tmp)[np.isin(t_gwflow[layer]['%.2f clays' % thickness],t_total_tmp)],endnodes=compaction_solver_debug_include_endnodes,preset_initial_maxstress=preset_initial_maxstress,ic_maxstress=initial_maxstress[layer]['%.2f clays' % thickness])
                else:
                    deformation[layer]['total_%.2f clays' % thickness],inelastic_flag_compaction[layer]['elastic_%.2f clays' % thickness]=subsidence_solver_aquitard_elasticinelastic(head_series[layer]['%.2f clays' % thickness][np.where(np.isin(head_series[layer]['%.2f clays' % thickness],t_total_tmp))],(clay_Sse[layer]-compressibility_of_water),(clay_Ssv[layer]-compressibility_of_water),np.max(Z[layer]['%.2f clays' % thickness]),n_z,endnodes=compaction_solver_debug_include_endnodes,preset_initial_maxstress=preset_initial_maxstress,ic_maxstress=initial_maxstress[layer]['%.2f clays' % thickness])
                deformation[layer]['total_%.2f clays' % thickness] = interbeds_distributions[layer][thickness] * deformation[layer]['total_%.2f clays' % thickness] 
                # deformation[layer]['elastic_%.2f clays' % thickness] = interbeds_distributions[layer][thickness] * deformation[layer]['elastic_%.2f clays' % thickness]
                # deformation[layer]['inelastic_%.2f clays' % thickness]= interbeds_distributions[layer][thickness] * deformation[layer]['elastic_%.2f clays' % thickness]
            # Now collect the results at the max timestep


            deformation_OUTPUT_tmp={}
            deformation_OUTPUT_tmp['dates'] = [x.strftime('%d-%b-%Y') for x in num2date(t_total_tmp)]

            deformation_OUTPUT_tmp['Interconnected Matrix'] = np.array(deformation[layer]['Interconnected matrix'])[np.where(np.isin(head_data[layer][:,0],t_total_tmp))]
            def_tot_tmp = np.zeros_like(t_total_tmp,dtype='float')               
            
            def_tot_tmp += np.array(deformation[layer]['Interconnected matrix'])[np.where(np.isin(head_data[layer][:,0],t_total_tmp))]
            for thickness in bed_thicknesses_tmp:
                def_tot_tmp += np.array(deformation[layer]['total_%.2f clays' % thickness])
                deformation_OUTPUT_tmp['total_%.2f clays' % thickness]=np.array(deformation[layer]['total_%.2f clays' % thickness])
            
            
            deformation[layer]['total'] = np.array([t_total_tmp,def_tot_tmp])
            deformation_OUTPUT_tmp['total']=def_tot_tmp
            deformation_OUTPUT[layer] = pd.DataFrame(deformation_OUTPUT_tmp)
            


    if layer_types[layer]=='Aquitard':
        print()
        print('%s is an Aquitard. Solving for layer compaction.' % layer)
        if layer_compaction_switch[layer]:
            if compaction_solver_compressibility_type[layer]=='elastic-inelastic':
                deformation[layer]={}
                # Find the timestep at which to solve
                print('\tCalculating deformation for layer %s. dt is %.2f.' % (layer, max_solver_dt))
                dt_max_tmp = max_solver_dt
                t_total_tmp = np.linspace(np.min(t_gwflow[layer]),np.min(t_gwflow[layer]) + int((np.max(t_gwflow[layer])-np.min(t_gwflow[layer]))/dt_max_tmp) * dt_max_tmp,num=int((np.max(t_gwflow[layer])-np.min(t_gwflow[layer]))/dt_max_tmp)+1) # this is the master t which will apply for all the solved layers. Note that, for certain dt_max_tmp which are larger than 1, the finish value may be smaller than the greatest input time.

                if overburden_stress_compaction:
                    unconfined_tmp = unconfined_aquifer_name==layer
                    print('UNCONFINED STATUS = %s' % unconfined_tmp)
                    if len(overburden_dates) != len(groundwater_solution_dates[layer]):
                        print('\t\t\tOverburden series is %i long whereas head series is %i long. Interpolating overburden stress.' % (len(overburden_dates),len(groundwater_solution_dates[layer])))
                        f_tmp = scipy.interpolate.interp1d(overburden_dates,overburden_data)
                        overburden_data_tmp = f_tmp(groundwater_solution_dates[layer])
                        overburden_dates_tmp = groundwater_solution_dates[layer]
                    else:
                        overburden_data_tmp = overburden_data

                    totdeftmp,inelastic_flag_compaction[layer]=subsidence_solver_aquitard_elasticinelastic(head_series[layer][:,np.isin(t_gwflow[layer],t_total_tmp)],(clay_Sse[layer]-compressibility_of_water),(clay_Ssv[layer]-compressibility_of_water),np.max(Z[layer]),n_z,unconfined=unconfined_tmp,overburden=overburden_stress_compaction,overburden_data=1/(rho_w * g) * np.array(overburden_data_tmp)[np.isin(t_gwflow[layer],t_total_tmp)],preset_initial_maxstress=preset_initial_maxstress,ic_maxstress=initial_maxstress[layer])
                    deformation[layer]['total'] = np.array([t_total_tmp,totdeftmp])

                else:
                    totdeftmp,inelastic_flag_compaction[layer]=subsidence_solver_aquitard_elasticinelastic(head_series[layer][:,np.isin(t_gwflow[layer],t_total_tmp)],(clay_Sse[layer]-compressibility_of_water),(clay_Ssv[layer]-compressibility_of_water),np.max(Z[layer]),n_z,preset_initial_maxstress=preset_initial_maxstress,ic_precons=initial_maxstress[layer])
                    deformation[layer]['total'] = np.array([t_total_tmp,totdeftmp])          

if MODE=='Normal': # If we are resuming, we do not scale layer thicknesses by default.
    if len(layers_var_thickness)>=1:
        print('')
        print('Scaling layer outputs by temporally varying layer thicknesses.')
        for layer in layers_var_thickness:
            print('\tScaling TOTAL outputs for %s.' % layer)
            prekeyname = np.array(list(layer_thicknesses[layer].keys()))[np.where(['pre' in key for key in list(layer_thicknesses[layer].keys())])[0][0]]
            datetimedates = num2date(deformation[layer]['total'][0,:])
    #        datetimedates = [dt.strptime(d,'%d-%b-%Y') for d in deformation_OUTPUT[layer]['dates'].values]
            logicaltmp = [datetimedate <= dt(int('%s' % prekeyname.split('-')[1]) ,9,1,tzinfo=datetime.timezone.utc) for datetimedate in datetimedates]
            scaling_factor_tmp = layer_thicknesses[layer][prekeyname]/initial_thicknesses[layer] 
            deformation_scaled_tmp =  deformation[layer]['total'][1,:][logicaltmp] * scaling_factor_tmp
    
            nonprekeynames = np.array(list(layer_thicknesses[layer].keys()))[np.where(['pre' not in key for key in list(layer_thicknesses[layer].keys())])[0]]
            nonprekeynames.sort()
            for key in nonprekeynames:
                if not key.endswith('-'):
                    scaling_factor_tmp = layer_thicknesses[layer][key]/initial_thicknesses[layer] 
                    years_tmp=key.split('-')
                    print('\t\tScaling years', years_tmp,'by ',scaling_factor_tmp)
                    logicaltmp = [(datetimedate <= dt(int('%s' % years_tmp[1]) ,9,1,tzinfo=datetime.timezone.utc)) and (datetimedate > dt(int('%s' % years_tmp[0]) ,9,1,tzinfo=datetime.timezone.utc))  for datetimedate in datetimedates]
                    deformation_scaled_tmp=np.append(deformation_scaled_tmp,( deformation[layer]['total'][1,:][logicaltmp]- deformation[layer]['total'][1,:][np.where(logicaltmp)[0][0]-1])  * scaling_factor_tmp +  deformation_scaled_tmp[-1])
            for key in nonprekeynames:
                if key.endswith('-'):
                    scaling_factor_tmp = layer_thicknesses[layer][key]/initial_thicknesses[layer] 
                    years_tmp=key.split('-')
                    print('\t\tScaling years', years_tmp,'by ',scaling_factor_tmp)
                    logicaltmp = [datetimedate > dt(int('%s' % years_tmp[0]) ,9,1,tzinfo=datetime.timezone.utc) for datetimedate in datetimedates]
                    deformation_scaled_tmp=np.append(deformation_scaled_tmp,( deformation[layer]['total'][1,:][logicaltmp]- deformation[layer]['total'][1,:][np.where(logicaltmp)[0][0]-1])  * scaling_factor_tmp + deformation_scaled_tmp[-1])
            deformation[layer]['total'][1,:]=deformation_scaled_tmp
            print('\tScaling SUBOUTPUTS for %s.' % layer)
            
            scaling_factor_tmp = layer_thicknesses[layer][prekeyname]/initial_thicknesses[layer] 
            deformation_scaled_tmp =  deformation[layer]['total'][1,:][logicaltmp] * scaling_factor_tmp

solving_compaction_stop = process_time()
solving_compaction_time = solving_compaction_stop - solving_compaction_start

#%% New section, saving compaction outputs.
print()
print()
print(''.center(80, '*'))
print('  SAVING COMPACTION SOLVER OUTPUTS  '.center(80, '*'))
print(''.center(80, '*'))
print()
saving_compaction_start = process_time()
time.sleep(internal_time_delay)


for layer in layer_names:
    if layer_types[layer]=='Aquitard':
        if layer_compaction_switch[layer]:
            print('Saving data for aquitard layer %s.' % layer)
            print('\tMaking deformation figure')
            sns.set_style('darkgrid')
            sns.set_context('talk')            
            # plt.figure(figsize=(18,12))
            # plt.plot_date(deformation[layer]['total'][0],deformation[layer]['total'][1])
            # # plt.plot_date(t,deformation[layer]['elastic'],label='elastic')
            # # plt.plot_date(t,deformation[layer]['inelastic'],label='inelastic')
            # plt.savefig('%s/figures/compaction_%s.png' % (outdestination, layer.replace(' ','_')),bbox_inches='tight')
            # plt.legend()
            # plt.xlim(date2num([date(2015,1,1),date(2020,1,1)]))
            # plt.close()
            # plt.savefig('%s/figures/compaction_%s_20152020.png' % (outdestination, layer.replace(' ','_')),bbox_inches='tight')
            
    
            if save_s:        
                print('\tSaving s timeseries')
                np.savetxt('%s/%s_s.csv' % (outdestination, layer.replace(' ','_')),deformation[layer]['total'])
    
    if layer_types[layer]=='Aquifer':
        if layer_compaction_switch[layer]:
            if save_internal_compaction:   
                    if np.size(inelastic_flag_compaction[layer]['elastic_%.2f clays' % thickness]) >= 3e6:
                        if gmt:
                            print('\t\t\tInelastic flag has more than 3 million entries; saving as signed char.')
                            inelastic_flag_compaction[layer]['elastic_%.2f clays' % thickness].astype(np.byte).tofile('%s/s_outputs/%s_%sclayinelastic_flag_COMPACTION' % (outdestination, layer.replace(' ','_'),thickness))
                            print('\t\t\t\tConverting to netCDF format. Command is:')
                            cmd_tmp="gmt xyz2grd %s/s_outputs/%s_%sclayinelastic_flag_COMPACTION -G%s/s_outputs/%s_%sclayinelastic_flag_COMPACTION.nb -I%i+n/%i+n -R%.3ft/%.3ft/%.3f/%.3f -ZTLc" % (outdestination, layer.replace(' ','_'),thickness,outdestination, layer.replace(' ','_'),'%.2f' % thickness,n_z,len(t_gwflow[layer]['%.2f clays' % thickness]),np.min(t_gwflow[layer]['%.2f clays' % thickness]),np.max(t_gwflow[layer]['%.2f clays' % thickness]),np.min(Z[layer]['%.2f clays' % thickness])+ np.diff(Z[layer]['%.2f clays' % thickness])[0]/2,np.max(Z[layer]['%.2f clays' % thickness])-np.diff(Z[layer]['%.2f clays' % thickness])[0]/2)
                            
                            print(cmd_tmp)
                            subprocess.call(cmd_tmp,shell=True)
                            os.remove('%s/s_outputs/%s_%sclayinelastic_flag_COMPACTION' % (outdestination, layer.replace(' ','_'),thickness))
        
                        else:
                            print('\t\t\tInelastic flag gwflow has more than 3 million entries; saving as signed char.')
                            inelastic_flag_tmp.astype(np.byte).tofile('%s/s_outputs/%s_%sclayinelastic_flag_COMPACTION' % (outdestination, layer.replace(' ','_'),thickness))
        
        
                    else:
                        with open('%s/s_outputs/%s_%sclayinelastic_flag_COMPACTION.csv' % (outdestination, layer.replace(' ','_'),thickness), "w+") as myCsv:
                            csvWriter = csv.writer(myCsv, delimiter=',')
                            csvWriter.writerows(inelastic_flag_tmp)

            
            print('Saving figures and data for aquifer layer %s.' % layer)
            if not os.path.isdir('%s/figures/%s' % (outdestination,layer)):
                os.mkdir('%s/figures/%s' % (outdestination,layer))
                           
            print('\tSaving layer data.')
            print('')
            print('printing def output thing')
            print(deformation_OUTPUT[layer])
            if len(layers_var_thickness)>=1:
                if layer in layers_var_thickness:
                    print('')
                    print('\tScaling sublayer outputs by temporally varying layer thicknesses.')
                    interbeds_tmp=interbeds_distributions[layer]
                    bed_thicknesses_tmp=list(interbeds_tmp.keys())
                    print('\t\t%s is an aquifer with interbedded clays. Scaling thickness of clays: %s' % (layer,bed_thicknesses_tmp))
                    print('\t\tFirst scaling the interbeds....')
                    for thickness in bed_thicknesses_tmp:                
                        prekeyname = np.array(list(layer_thicknesses[layer].keys()))[np.where(['pre' in key for key in list(layer_thicknesses[layer].keys())])[0][0]]
                        datetimedates = [dt.strptime(d,'%d-%b-%Y') for d in deformation_OUTPUT[layer]['dates'].values]
                        logicaltmp = [datetimedate <= dt(int('%s' % prekeyname.split('-')[1]) ,9,1) for datetimedate in datetimedates]
                        
                        scaling_factor_tmp = layer_thicknesses[layer][prekeyname]/initial_thicknesses[layer] 
                        deformation_scaled_tmp =  deformation_OUTPUT[layer]['total_%.2f clays' % thickness][logicaltmp].values * scaling_factor_tmp
                        nonprekeynames = np.array(list(layer_thicknesses[layer].keys()))[np.where(['pre' not in key for key in list(layer_thicknesses[layer].keys())])[0]]
                        nonprekeynames.sort()
                        for key in nonprekeynames:
                            if not key.endswith('-'):
                                scaling_factor_tmp = layer_thicknesses[layer][key]/initial_thicknesses[layer] 
                                years_tmp=key.split('-')
                                print('\t\tScaling years', years_tmp,'by ',scaling_factor_tmp)
                                logicaltmp = [(datetimedate <= dt(int('%s' % years_tmp[1]) ,9,1,)) and (datetimedate > dt(int('%s' % years_tmp[0]) ,9,1))  for datetimedate in datetimedates]
                                deformation_scaled_tmp=np.append(deformation_scaled_tmp,( deformation_OUTPUT[layer]['total_%.2f clays' % thickness][logicaltmp]- deformation_OUTPUT[layer]['total_%.2f clays' % thickness][np.where(logicaltmp)[0][0]-1])  * scaling_factor_tmp +  deformation_scaled_tmp[-1])
                        for key in nonprekeynames:
                            if key.endswith('-'):
                                scaling_factor_tmp = layer_thicknesses[layer][key]/initial_thicknesses[layer] 
                                years_tmp=key.split('-')
                                print('\t\tScaling years', years_tmp,'by ',scaling_factor_tmp)
                                logicaltmp = [datetimedate > dt(int('%s' % years_tmp[0]) ,9,1) for datetimedate in datetimedates]
                                deformation_scaled_tmp=np.append(deformation_scaled_tmp,( deformation_OUTPUT[layer]['total_%.2f clays' % thickness][logicaltmp]- deformation_OUTPUT[layer]['total_%.2f clays' % thickness][np.where(logicaltmp)[0][0]-1])  * scaling_factor_tmp + deformation_scaled_tmp[-1])
                        deformation_OUTPUT[layer]['total_%.2f clays' % thickness]=deformation_scaled_tmp

                    print('\t\tNow scaling the elastic, and recomputing the total')
                    scaling_factor_tmp = layer_thicknesses[layer][prekeyname]/initial_thicknesses[layer] 
                    logicaltmp = [datetimedate <= dt(int('%s' % prekeyname.split('-')[1]) ,9,1) for datetimedate in datetimedates]
                    deformation_scaled_tmp =  deformation_OUTPUT[layer]['Interconnected Matrix'][logicaltmp].values * scaling_factor_tmp
                    nonprekeynames = np.array(list(layer_thicknesses[layer].keys()))[np.where(['pre' not in key for key in list(layer_thicknesses[layer].keys())])[0]]
                    nonprekeynames.sort()
                    for key in nonprekeynames:
                        if not key.endswith('-'):
                            scaling_factor_tmp = layer_thicknesses[layer][key]/initial_thicknesses[layer] 
                            years_tmp=key.split('-')
                            print('\t\tScaling years', years_tmp,'by ',scaling_factor_tmp)
                            logicaltmp = [(datetimedate <= dt(int('%s' % years_tmp[1]) ,9,1,)) and (datetimedate > dt(int('%s' % years_tmp[0]) ,9,1))  for datetimedate in datetimedates]
                            deformation_scaled_tmp=np.append(deformation_scaled_tmp,( deformation_OUTPUT[layer]['Interconnected Matrix'][logicaltmp]- deformation_OUTPUT[layer]['Interconnected Matrix'][np.where(logicaltmp)[0][0]-1])  * scaling_factor_tmp +  deformation_scaled_tmp[-1])
                    for key in nonprekeynames:
                        if key.endswith('-'):
                            scaling_factor_tmp = layer_thicknesses[layer][key]/initial_thicknesses[layer] 
                            years_tmp=key.split('-')
                            print('\t\tScaling years', years_tmp,'by ',scaling_factor_tmp)
                            logicaltmp = [datetimedate > dt(int('%s' % years_tmp[0]) ,9,1) for datetimedate in datetimedates]
                            deformation_scaled_tmp=np.append(deformation_scaled_tmp,( deformation_OUTPUT[layer]['Interconnected Matrix'][logicaltmp]- deformation_OUTPUT[layer]['Interconnected Matrix'][np.where(logicaltmp)[0][0]-1])  * scaling_factor_tmp + deformation_scaled_tmp[-1])
                    deformation_OUTPUT[layer]['Interconnected Matrix']=deformation_scaled_tmp

                    def_tot_tmp = np.zeros_like(deformation_scaled_tmp,dtype='float')               
            
                    def_tot_tmp += np.array(deformation_OUTPUT[layer]['Interconnected Matrix'])
                    for thickness in bed_thicknesses_tmp:
                        def_tot_tmp += np.array(deformation_OUTPUT[layer]['total_%.2f clays' % thickness])
                        
                    deformation_OUTPUT[layer]['total']=def_tot_tmp
            
            deformation_OUTPUT[layer].to_csv('%s/%s_Total_Deformation_Out.csv' % (outdestination,layer),index=False)
                
            
            print('\tSaving layer compaction figure')
            sns.set_style('whitegrid')
            sns.set_context('talk')
            l_aqt=[]

            plt.figure(figsize=(18,12))
            line_tmp, = plt.plot_date(deformation[layer]['total'][0,:],deformation_OUTPUT[layer]['Interconnected Matrix'],'-',label='Interconnected matrix')
            l_aqt.append(line_tmp)

            interbeds_tmp=interbeds_distributions[layer]
            bed_thicknesses_tmp=list(interbeds_tmp.keys())

            for thickness in bed_thicknesses_tmp:
                line_tmp, = plt.plot_date([dt.strptime(d,'%d-%b-%Y') for d in deformation_OUTPUT[layer]['dates'].values],deformation_OUTPUT[layer]['total_%.2f clays' % thickness],'-',label='%s_%ix%.2f clays' % (layer,interbeds_distributions[layer][thickness],thickness))
                l_aqt.append(line_tmp)

            line_tmp, = plt.plot_date(deformation[layer]['total'][0,:],deformation_OUTPUT[layer]['total'],'-',label='total')
            l_aqt.append(line_tmp)
            plt.title('%s' % layer)
            plt.ylabel('Deformation (m)')
            plt.legend()
            plt.savefig('%s/figures/%s/overall_compaction_%s.png' % (outdestination,layer,layer),bbox_inches='tight')
#            if np.min(line_tmp.get_xdata()) <= date2num(date(2015,1,1)):
#                for line in l_aqt:
#                    line.set_ydata(np.array(line.get_ydata()) - np.array(line.get_ydata())[np.array(line.get_xdata())==date2num(date(2015,1,1))])
#    
#                plt.xlim(date2num([date(2015,1,1),date(2020,1,1)]))
#                
#                plt.savefig('%s/figures/%s/overall_compaction_%s_201520.png' % (outdestination,layer,layer),bbox_inches='tight')
            plt.close() 
            
            if save_s:
                print('\tSaving s (interconnected matrix)')
                
                if np.size(deformation[layer]['Interconnected matrix']) >= 3e6:
                    print('\t\t\ts (interconnected matrix) has more than 1 million entries; saving as 32 bit floats.')
                    deformation[layer]['Interconnected matrix'].astype(np.single).tofile('%s/s_outputs/%s_s_matrix' % (outdestination, layer.replace(' ','_')))
    
    
                else:            
                    np.savetxt('%s/s_outputs/%s_s_matrix.csv' % (outdestination, layer.replace(' ','_')),deformation[layer]['Interconnected matrix'])
                
                print('\tSaving s (clay layers)')
                for thickness in bed_thicknesses_tmp:
                    print('\t\t%.2f' % thickness)
                    if np.size(deformation[layer]['total_%.2f clays' % thickness]) >= 1e6:
                        print('\t\t\ts (clay) has more than 1 million entries; saving as 32 bit floats.')
                        deformation[layer]['total_%.2f clays' % thickness].astype(np.single).tofile('%s/s_outputs/%s_s_%.2fclays' % (outdestination, layer.replace(' ','_'),thickness))
                    else:                   
                        np.savetxt('%s/s_outputs/%s_s_%.2fclays.csv' % (outdestination, layer.replace(' ','_'),thickness),deformation[layer]['total_%.2f clays' % thickness])

          
print("Creating overall compaction plot and saving deformation series")
sns.set_style('whitegrid')
plt.figure(figsize=(18,12))
#dt_master_compacting_layers = {key:dt_master[key] for key in compacting_layers} # OLD kept for "just in case"
#maxdt = max(dt_master_compacting_layers.values())
maxdt = max_solver_dt
print('\tmax dt (output dt) = %.2f' % maxdt)
deformation_OUTPUT2 = {}
t_total_tmp = np.arange(starttime,endtime,maxdt) # This might lead to problems if maxdt is small, i.e. much less than 1, due to rounding errors which occur in np.arange.
deformation_OUTPUT2['dates']=[x.strftime('%d-%b-%Y') for x in num2date(t_total_tmp)]
t_overall = np.zeros_like(t_total_tmp,dtype=float)
l_aqt=[]
l_aqf=[]
for layer in layer_names:
    if layer_compaction_switch[layer]:
        if layer_types[layer]=='Aquitard':
            #dt_tmp = int(maxdt/dt_master[layer])
#            t_overall = t_overall + deformation[layer]['total'][1,::dt_tmp]
            l_tmp, = plt.plot_date(deformation[layer]['total'][0,:],deformation[layer]['total'][1,:],'-',label='%s' % layer)
            l_aqt.append(l_tmp)
            deformation_OUTPUT2[layer]=deformation[layer]['total'][1,:][np.isin(deformation[layer]['total'][0,:],t_total_tmp)]
             
        if layer_types[layer]=='Aquifer':
            l_tmp, = plt.plot_date(deformation[layer]['total'][0,:],deformation[layer]['total'][1,:],'-',label='%s' % layer)
            l_aqf.append(l_tmp)
            #dt_tmp = int(maxdt/dt_master[layer])
            deformation_OUTPUT2[layer]=deformation[layer]['total'][1,:][np.isin(deformation[layer]['total'][0,:],t_total_tmp)]

# Add up all the deformations from each layer
            
def_tot_tmp = np.zeros_like(t_total_tmp,dtype='float')
for layer in layer_names:
    if layer_compaction_switch[layer]:
        newtot = np.array(deformation[layer]['total'][1,:])[np.isin(deformation[layer]['total'][0,:],t_total_tmp)]
        t_overall = t_overall + newtot

deformation_OUTPUT2['Total']=t_overall
def_out = pd.DataFrame(deformation_OUTPUT2)
def_out.to_csv('%s/Total_Deformation_Out.csv' % outdestination,index=False)

l3, = plt.plot_date(t_total_tmp,t_overall,label='TOTAL def')

plt.ylabel('Z (m)')
plt.legend()
plt.savefig('%s/figures/total_deformation_figure.png' % outdestination,bbox_inches='tight')

saving_compaction_stop = process_time()
saving_compaction_time = saving_compaction_stop - saving_compaction_start

t_total_stop = process_time()
t_total = t_total_stop - t_total_start

plt.figure(figsize=(18,18))
plt.pie(np.abs([param_read_time, solving_compaction_time, saving_head_time, reading_head_time, solving_head_time, saving_compaction_time,t_total - np.sum([param_read_time, solving_compaction_time, saving_head_time, reading_head_time, solving_head_time, saving_compaction_time])]),labels=['param_read_time', 'solving_compaction_time', 'saving_head_time', 'reading_inputs_time', 'solving_head_time', 'saving_compaction_time','misc'],autopct=lambda p : '{:.2f}%  ({:,.0f})'.format(p,p * t_total/100))
plt.title('Total runtime = %i seconds' % t_total)
plt.savefig('%s/figures/runtime_breakdown.png' % outdestination,bbox_inches='tight')
plt.close()


print('Model Run Complete')
