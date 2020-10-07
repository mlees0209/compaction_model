#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### This is the main model script. This is where all the "under the hood" operations occur. Do not edit unless you know what you are doing. The correct call should be: python execute_model.py parameter_file.par

#from __future__ import print_function
from model_functions import *
from matplotlib.dates import num2date
import datetime
import csv
from datetime import datetime as dt
import seaborn as sns
import subprocess
import distutils.util as ut
from time import process_time 
import matplotlib.dates as mdates
import scipy 
from datetime import date
import matplotlib.colors as colors
import operator
from astropy.convolution import convolve

t_total_start = process_time()

gmt=True # If this is true, head outputs can be saved as .netCDF grid files which are tiny. It uses gmt xyz2grd command, hence requires gmt to be installed.

output_folder='.'
run_name='.'
sys.stdout = Logger(output_folder,run_name)

if len(sys.argv) != 2:
    print('Execute model error; terminal. Incorrect number of input arguments. Correct usage: python execute_model.py parameter_file.par')
    sys.exit(1)

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

internal_time_delay = read_parameter('internal_time_delay',float,1,paramfilelines)
overwrite=read_parameter('overwrite',bool,1,paramfilelines)
run_name = read_parameter('run_name',str,1,paramfilelines)
output_folder = read_parameter('output_folder',str,1,paramfilelines)
outdestination="%s/%s" % (output_folder,run_name)

if not os.path.isdir(outdestination):
    print("\t\tMaking output directory at %s. This stdout will now be directed to a log in that folder as well as displayed here." % (outdestination))
    os.mkdir(outdestination)
    os.mkdir('%s/figures' % outdestination)
    os.mkdir('%s/s_outputs' % outdestination)
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
        else:
            print('\t\tNot overwriting. Aborting.' % check)
            sys.exit(1)
    #    OVERWRITE = input("\t\tOutput directory %s already exists. Do you want to overwrite this directory? WARNING: may delete existing data." % (outdestination))


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
layer_thicknesses=read_parameter('layer_thicknesses',float,no_layers,paramfilelines)
layer_compaction_switch=read_parameter('layer_compaction_switch',bool,no_layers,paramfilelines)
interbeds_switch=read_parameter('interbeds_switch',bool,list(layer_types.values()).count('Aquifer'),paramfilelines)
#interbeds_type=read_parameter('interbeds_type',str,list(layer_types.values()).count('Aquifer'),paramfilelines)
# Import interbeds_distributions -- an awkward parameter as its a dictionary of dictionaries!
interbeds_distributions1=read_parameter('interbeds_distributions',dict,sum(interbeds_switch.values()),paramfilelines)
interbeds_distributions1=np.array(interbeds_distributions1)
if np.shape(interbeds_distributions1)[0]==1:
    interbeds_distributions1=interbeds_distributions1[0]
    minidics = [dict([(float(re.split(',|:',interbeds_distributions1[2*i + 1])[2*j]),int( re.split(',|:',interbeds_distributions1[2*i + 1])[2*j+1])) for j in range(int(len( re.split(',|:',interbeds_distributions1[2*i + 1]))/2))]) for i in range(sum(interbeds_switch.values()))]
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
dt_master = read_parameter('dt_master',float,no_layers_containing_clay,paramfilelines)
dz_clays = read_parameter('dz_clays',float,no_layers_containing_clay,paramfilelines)
vertical_conductivity = read_parameter('vertical_conductivity',float,len(layers_requiring_solving),paramfilelines)
overburden_stress_compaction = read_parameter('overburden_stress_compaction',bool,1,paramfilelines)
#overburden_compaction = read_parameter('overburden_compaction',bool,1,paramfilelines)
if overburden_stress_gwflow or overburden_stress_compaction: # Only used if we're doing overburden anywhere
    specific_yield = read_parameter('specific_yield',float,1,paramfilelines)


param_read_stop = process_time()
param_read_time = param_read_start - param_read_stop
### Next section is "READING IN HEAD MODULE"
print()
print()
print(''.center(80, '*'))
print('  READING HEAD DATA  '.center(80, '*'))
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
        os.mkdir('%s/input_head_data' % outdestination)
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
            copy2(fileloc,'%s/input_head_data/%s' % (outdestination,fileloc.split('/')[-1]))
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
            if overburden_stress_gwflow or overburden_stress_compaction:
                if aquifer == layer_names[0]:
                    print('\t%s is the uppermost aquifer; therefore it is unconfined. Calculating overburden stress from changing water levels in this aquifer.' % aquifer)
                    unconfined_aquifer_name = aquifer
                    overburden_dates = dates
                    overburden_data = [specific_yield * rho_w * g * (data[i] - data[0]) for i in range(len(data))]
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
for aquifer in all_aquifers_needing_head_data:
    idx_to_keep = ((head_data[aquifer][:,0]>=starttime) & (head_data[aquifer][:,0]<=endtime)) 
    datesnew = head_data[aquifer][:,0][idx_to_keep]
    datanew = head_data[aquifer][:,1][idx_to_keep]
    head_data[aquifer]=np.array([datesnew,datanew]).T
    with open('%s/input_head_data/input_time_series_%s.csv' % (outdestination, aquifer.replace(' ','_')), "w+") as myCsv:
        csvWriter = csv.writer(myCsv, delimiter=',')
        csvWriter.writerows(head_data[aquifer])
print('Clipping done.')

if overburden_stress_gwflow or overburden_stress_compaction:
    print('Clipping overburden stress.')
    idx_to_keep = ((overburden_dates>=starttime) & (overburden_dates<=endtime)) 
    overburden_dates = overburden_dates[idx_to_keep]
    overburden_data = np.array(overburden_data)[idx_to_keep]
    print('Clipping done.')
    plt.plot_date(overburden_dates,overburden_data)
    plt.savefig('%s/input_head_data/overburden_stress_series.png' % outdestination)
    print('\tOverburden stress calculated and saved in input_head_data.')                    
effective_stress={}


plt.figure(figsize=(18,12))
sns.set_style('darkgrid')
for aquifer in all_aquifers_needing_head_data:
    plt.plot_date(head_data[aquifer][:,0],head_data[aquifer][:,1],label='%s' % aquifer)
plt.ylabel('Head (masl)')
plt.legend()
plt.savefig('%s/input_head_data/inputtimeseries.png' % outdestination)
plt.close()
sns.set_style('white')

reading_head_stop = process_time()
reading_head_time = reading_head_stop - reading_head_start

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
Z={}
t_gwflow={}

head_series=copy.deepcopy(head_data)


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
            print('\t\tHead time series required for overlying layer %s and lower layer %s.' % (top_boundary,bot_boundary))
            if top_boundary in head_data.keys() and bot_boundary in head_data.keys():
                print('\t\t\tHead time series found.')

            # check if the in-series for these layers are equally spaced with dt = dt_master
            if layer not in dt_master.keys():
                print('\t\t\tSolving head series error: TERMINAL. dt_master not specified for layer %s. EXITING.' % layer)
                sys.exit(1)
                
            t_top=head_data[top_boundary][:,0]
            dt_tmp = np.diff(t_top)
            test1 = [n.is_integer() for n in dt_master[layer]/dt_tmp]
            test2 = [n.is_integer() for n in dt_tmp/dt_master[layer]]
            if ((not np.all(test1) and not np.all(test2)) or (not list(dt_tmp).count(dt_tmp[0])) == len(dt_tmp)):
                print('\t\t\tSolving head series error: TERMINAL. dt_master not compatible with dt in the %s input series. dt_master must be an integer multiple of dt in the time series, and dt in the time series must be constant. EXITING.' % top_boundary)
                sys.exit(1)
            if dt_master[layer]>=dt_tmp[0]:
                spacing_top = int(dt_master[layer]/dt_tmp[0])
            else:
                print('\t\t\tNOTE: dt_master < dt_data. Linear resampling of input head series occuring.')
                t_interp_new = np.arange(min(t_top),max(t_top)+0.00001,dt_master[layer])
                f_tmp = scipy.interpolate.interp1d(t_top,head_data[top_boundary][:,1])
                head_data[top_boundary] = np.array([t_interp_new,f_tmp(t_interp_new)]).T
                spacing_top=1

            t_bot=head_data[bot_boundary][:,0]
            dt_tmp = np.diff(t_bot)
            test1 = [n.is_integer() for n in dt_master[layer]/dt_tmp]
            test2 = [n.is_integer() for n in dt_tmp/dt_master[layer]]
            if ((not np.all(test1) and not np.all(test2)) or (not list(dt_tmp).count(dt_tmp[0])) == len(dt_tmp)):
                print('\t\t\tSolving head series error: TERMINAL. dt_master not compatible with dt in the %s input series. dt_master must be an integer multiple of dt in the time series. EXITING.' % top_boundary)
                sys.exit(1)
            if dt_master[layer]>=dt_tmp[0]:
                spacing_bot = int(dt_master[layer]/dt_tmp[0])
            else:
                print('\t\t\tNOTE: dt_master < dt_data. Linear resampling of input head series occuring.')
                t_interp_new = np.arange(min(t_bot),max(t_bot)+0.000001,dt_master[layer])
                f_tmp = scipy.interpolate.interp1d(t_top,head_data[bot_boundary][:,1])
                head_data[bot_boundary] = np.array([t_interp_new,f_tmp(t_interp_new)]).T
                spacing_bot=1

            if not all(t_top == t_bot):
                print('\t\t\tSolving head series error: TERMINAL. Time series in %s and %s aquifers have different dates.' % (top_boundary,bot_boundary))
                sys.exit(1)
            else:
                print('\t\t\tTime series found with correct dt and over same timespan.')          

            z=np.arange(0,layer_thicknesses[layer]+0.00001,dz_clays[layer]) # 0.000001 to include the stop value.
            Z[layer]=z
            
            t_in = np.arange(np.min(t_top),np.max(t_top)+0.0001,dt_master[layer]) # 0.000001 to include the stop value.
            
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
            
            initial_condition_tmp = (head_data[top_boundary][0,1] + head_data[bot_boundary][0,1]) / 2
            
            if groundwater_flow_solver_type[layer] == 'singlevalue':
                hmat=solve_head_equation_singlevalue(dt_master[layer],t_in,dz_clays[layer],z,np.vstack((head_data[top_boundary][::spacing_top,1],head_data[bot_boundary][::spacing_bot,1])),initial_condition_tmp,vertical_conductivity[layer]/(clay_Ssk[layer]+compressibility_of_water))
            elif groundwater_flow_solver_type[layer] == 'elastic-inelastic':
                t1_start = process_time() 
                hmat,inelastic_flag_tmp=solve_head_equation_elasticinelastic(dt_master[layer],t_in,dz_clays[layer],z,np.vstack((head_data[top_boundary][::spacing_top,1],head_data[bot_boundary][::spacing_bot,1])),initial_condition_tmp,vertical_conductivity[layer]/clay_Sse[layer],vertical_conductivity[layer]/clay_Ssv[layer],overburdenstress=overburden_stress_gwflow,overburden_data=1/(rho_w * g) * np.array(overburden_data_tmp))
                t1_stop = process_time() 
                print("\t\t\tElapsed time in seconds:",  t1_stop-t1_start)  

            head_series[layer]=hmat
            inelastic_flag[layer] = inelastic_flag_tmp
            effective_stress[layer] = np.tile(overburden_data_tmp, (np.shape(hmat)[0],1)) -  rho_w*g*hmat 

            with open('%s/%s_inelastic_flag_DEBUG.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                csvWriter = csv.writer(myCsv, delimiter=',')
                csvWriter.writerows(inelastic_flag_tmp)
            #dateslist = [x.strftime('%d-%b-%Y') for x in num2date(t_in)]
            groundwater_solution_dates[layer]=t_in

            if overburden_stress_gwflow:
                effective_stress[layer] = np.tile(overburden_data_tmp, (np.shape(hmat)[0],1)) -  rho_w*g*hmat 
            else:
                effective_stress[layer] = np.zeros_like(hmat) -  rho_w*g*hmat                             


            if save_effective_stress:
                print('\t\tSaving effective stress and overburden stress outputs.')
                if overburden_stress_gwflow:
                    with open('%s/%s_overburden_stress.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                        csvWriter = csv.writer(myCsv, delimiter=',')
                        csvWriter.writerows(np.tile(overburden_data_tmp, (len(z_tmp),1)))
                else:
                    with open('%s/%s_overburden_stress.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                        csvWriter = csv.writer(myCsv, delimiter=',')
                        csvWriter.writerows(np.zeros_like(hmat))

                if np.size(effective_stress[layer]) >= 5e6:
                    print('\t\t\tEffective stress has more than 5 million entries; saving as 16 bit floats.')
                    effective_stress[layer].astype(half).tofile('%s/%s_effective_stress' % (outdestination, layer.replace(' ','_')))

                else:
                    
                    with open('%s/%s_effective_stress.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                        csvWriter = csv.writer(myCsv, delimiter=',')
                        csvWriter.writerows(effective_stress[layer])



        else:
            if layer_types[layer]=='Aquifer':
                head_series[layer]=dict([('Interconnected matrix',head_data[layer])])
                inelastic_flag[layer]={}
                groundwater_solution_dates[layer]={}
                Z[layer]={}
                t_gwflow[layer]={}
                effective_stress[layer]={}
                print('\t\t%s is an aquifer.' % layer)
                if interbeds_switch[layer]:
                    interbeds_tmp=interbeds_distributions[layer]
                    bed_thicknesses_tmp=list(interbeds_tmp.keys())
                    print('\t\t%s is an aquifer with interbedded clays. Thicknesses of clays to be solved are %s' % (layer,bed_thicknesses_tmp))
                    for thickness in bed_thicknesses_tmp:
                        print('')
                        print('\t\tSolving for thickness %.2f.' % thickness)
                        # This bit interpolated if dt_master < dt_boundarycondition 
                        t_aquifer_tmp=head_data[layer][:,0]
                        h_aquifer_tmp=head_data[layer][:,1]
                        z_tmp = np.arange(0,thickness+0.00001,dz_clays[layer]) # 0.000001 to include the stop value.
                        t_interp_new = 0.0001*np.arange(10000*min(t_aquifer_tmp),10000*max(t_aquifer_tmp)+1,10000*dt_master[layer]) # The ridiculous factor of 1/10000 is to ensure the np.arange function takes integer steps. Else, floating point precision issues mean that things go wrong for timesteps less than <0.1 seconds.
                        f_tmp = scipy.interpolate.interp1d(t_aquifer_tmp,h_aquifer_tmp)
                        h_aquifer_tmp_interpolated = np.array([t_interp_new,f_tmp(t_interp_new)]).T

                        initial_condition_tmp=h_aquifer_tmp[0]
                        
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
                                print('\t\t\tThis is the unconfined aquifer; overburden=0.')
                                overburden_data_tmp = np.zeros_like(t_interp_new)
                        else:
                            overburden_data_tmp=[0]

                        t1_start = process_time() 
                        hmat_tmp,inelastic_flag_tmp=solve_head_equation_elasticinelastic(dt_master[layer],t_interp_new,dz_clays[layer],z_tmp,np.vstack((h_aquifer_tmp_interpolated[:,1],h_aquifer_tmp_interpolated[:,1])),initial_condition_tmp,vertical_conductivity[layer]/clay_Sse[layer],vertical_conductivity[layer]/clay_Ssv[layer],overburdenstress=overburden_stress_gwflow,                        overburden_data=1/(rho_w * g) * np.array(overburden_data_tmp))
                        t1_stop = process_time() 
                        print("\t\t\tElapsed time in seconds:",  t1_stop-t1_start)  
                        head_series[layer]['%.2f clays' % thickness]=hmat_tmp
                        inelastic_flag[layer]['%.2f clays' % thickness] = inelastic_flag_tmp
                        t_gwflow[layer]['%.2f clays' % thickness] = t_interp_new
                        Z[layer]['%.2f clays' % thickness]=z_tmp

                        if np.size(inelastic_flag_tmp) >= 3e6:
                            if gmt:
                                print('\t\t\tInelastic flag gwflow has more than 3 million entries; saving as signed char.')
                                inelastic_flag_tmp.astype(np.byte).tofile('%s/%s_%sclayinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_'),thickness))
                                print('\t\t\t\tConverting to netCDF format. Command is:')
                                cmd_tmp="gmt xyz2grd %s/%s_%sclayinelastic_flag_GWFLOW -G%s/%s_%sclayinelastic_flag_GWFLOW.nb -I%.2f/%.5f -R%.2ft/%.2ft/%.2f/%.2f -ZTLc" % (outdestination, layer.replace(' ','_'),thickness,outdestination, layer.replace(' ','_'),thickness,dt_master[layer],np.diff(Z[layer]['%.2f clays' % thickness])[0],np.min(t_gwflow[layer]['%.2f clays' % thickness]),np.max(t_gwflow[layer]['%.2f clays' % thickness]),np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness]))
                                
                                print(cmd_tmp)
                                subprocess.call(cmd_tmp,shell=True)
                                os.remove('%s/%s_%sclayinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_'),thickness))

                            else:
                                print('\t\t\tInelastic flag gwflow has more than 3 million entries; saving as signed char.')
                                inelastic_flag_tmp.astype(np.byte).tofile('%s/%s_%sclayinelastic_flag_GWFLOW' % (outdestination, layer.replace(' ','_'),thickness))


                        else:
                            with open('%s/%s_%sclayinelastic_flag_GWFLOW.csv' % (outdestination, layer.replace(' ','_'),thickness), "w+") as myCsv:
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
                                    cmd_tmp="gmt xyz2grd %s/%s_%sclayeffective_stress -G%s/%s_%sclayeffective_stress.nc -I%.2f/%.5f -R%.2ft/%.2ft/%.2f/%.2f -ZTLf" % (outdestination, layer.replace(' ','_'),thickness,outdestination, layer.replace(' ','_'),thickness,dt_master[layer],np.diff(Z[layer]['%.2f clays' % thickness])[0],np.min(t_gwflow[layer]['%.2f clays' % thickness]),np.max(t_gwflow[layer]['%.2f clays' % thickness]),np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness]))
                            
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
                                        cmd_tmp="gmt xyz2grd %s/%s_%sclay_overburden_stress -G%s/%s_%sclay_overburden_stress.nc -I%.2f/%.5f -R%.2ft/%.2ft/%.2f/%.2f -ZTLf" % (outdestination, layer.replace(' ','_'),thickness,outdestination, layer.replace(' ','_'),thickness,dt_master[layer],np.diff(Z[layer]['%.2f clays' % thickness])[0],np.min(t_gwflow[layer]['%.2f clays' % thickness]),np.max(t_gwflow[layer]['%.2f clays' % thickness]),np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness]))
                            
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
    os.mkdir('%s/head_outputs' % outdestination)
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
                            cmd_tmp="gmt xyz2grd %s/head_outputs/%s_%sclay_head_data -G%s/head_outputs/%s_%sclay_head_data.nc -I%.2f/%.5f -R%.2ft/%.2ft/%.2f/%.2f -ZTLf" % (outdestination, layer.replace(' ','_'),thickness,outdestination, layer.replace(' ','_'),thickness,dt_master[layer],np.diff(Z[layer]['%.2f clays' % thickness])[0],np.min(t_gwflow[layer]['%.2f clays' % thickness]),np.max(t_gwflow[layer]['%.2f clays' % thickness]),np.min(Z[layer]['%.2f clays' % thickness]),np.max(Z[layer]['%.2f clays' % thickness]))
                            
                            print(cmd_tmp)
                            subprocess.call(cmd_tmp,shell=True)
                            os.remove('%s/head_outputs/%s_%sclay_head_data' % (outdestination, layer.replace(' ','_'),thickness))
                        else:
                            print('\t\t\tHead has more than 1 million entries; saving as 16 bit floats.')
                            head_series[layer]['%.2f clays' % thickness].astype(np.half).tofile('%s/head_outputs/%s_%sclay_head_data' % (outdestination, layer.replace(' ','_'),thickness))

                    else:
                        with open('%s/head_outputs/%s_%sclay_head_data.csv' % (outdestination, layer.replace(' ','_'),thickness), "w+") as myCsv:
                            csvWriter = csv.writer(myCsv, delimiter=',')
                            csvWriter.writerows(head_series[layer]['%.2f clays' % thickness])
            with open('%s/head_outputs/%s_groundwater_solution_dates.csv' % (outdestination, layer.replace(' ','_')), 'w') as myfile:
                wr = csv.writer(myfile)
                res = list(groundwater_solution_dates[layer].keys())[0] 
                wr.writerow([x.strftime('%c') for x in num2date(groundwater_solution_dates[layer][res])])

            
        if layer_types[layer]=='Aquitard':
            with open('%s/head_outputs/%s_head_data.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
                csvWriter = csv.writer(myCsv, delimiter=',')
                csvWriter.writerows(head_series[layer])
            with open('%s/head_outputs/%s_groundwater_solution_dates.csv' % (outdestination, layer.replace(' ','_')), 'w') as myfile:
                wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
                wr.writerow([x.strftime('%d-%b-%Y') for x in num2date(groundwater_solution_dates[layer])])

for layer in layers_requiring_solving:
    if create_output_head_video[layer]:
        print('create_output_head_video = True. Creating head timeseries video for specified layers. Note: requires ffmpeg installed.')
        if layer_types[layer]=='Aquitard':
                print('\tCreating video for %s.' % layer)
                hmat_tmp = head_series[layer]
                inelastic_flag_vid = inelastic_flag[layer]
                inelastic_flag_vid = inelastic_flag_vid==1            
                dates_str = [x.strftime('%d-%b-%Y') for x in num2date(groundwater_solution_dates[layer])]
                create_head_video_elasticinelastic(hmat_tmp,Z[layer],inelastic_flag_vid,dates_str,outdestination+'/figures',layer)

        if layer_types[layer]=='Aquifer':
            print('\tCreating video for %s.' % layer)
            interbeds_tmp=interbeds_distributions[layer]
            bed_thicknesses_tmp=list(interbeds_tmp.keys())
            print('\t\t%s is an aquifer with interbedded clays. Thicknesses of clays to make videos are %s' % (layer,bed_thicknesses_tmp))
            for thickness in bed_thicknesses_tmp:
                print('\t\tCreating video for %s_%.2f clays' % (layer,thickness))
                hmat_tmp = head_series[layer]['%.2f clays' % thickness]
                inelastic_flag_vid = inelastic_flag[layer]['%.2f clays' % thickness]
                inelastic_flag_vid = inelastic_flag_vid==1            
                dates_str = [x.strftime('%d-%b-%Y') for x in num2date(groundwater_solution_dates[layer]['%.2f clays' % thickness])]
                create_head_video_elasticinelastic(hmat_tmp,Z[layer]['%.2f clays' % thickness],inelastic_flag_vid,dates_str,outdestination+'/figures','%s_%.2f_clays' % (layer, thickness),delt=30)

saving_head_stop = process_time()
saving_head_time = saving_head_stop - saving_head_start

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
            layer_sand_thickness_tmp = layer_thicknesses[layer] - np.sum([list(interbeds_distributions[layer].keys())[i] * list(interbeds_distributions[layer].values())[i] for i in range(len(interbeds_distributions[layer]))])
            print('\tTotal sand thickness in aquifer is %.2f m.' % layer_sand_thickness_tmp)
            deformation[layer]['Interconnected matrix']=[layer_sand_thickness_tmp*(sand_Sse[layer]-compressibility_of_water)*(head_series[layer]['Interconnected matrix'][i,1] - head_series[layer]['Interconnected matrix'][0,1]) for i in range(len(head_series[layer]['Interconnected matrix'][:,1]))]
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

                    db[layer]['total_%.2f clays' % thickness],deformation[layer]['total_%.2f clays' % thickness],deformation[layer]['elastic_%.2f clays' % thickness],deformation[layer]['inelastic_%.2f clays' % thickness]=subsidence_solver_aquitard_elasticinelastic(head_series[layer]['%.2f clays' % thickness],inelastic_flag[layer]['%.2f clays' % thickness],(clay_Sse[layer]-compressibility_of_water),(clay_Ssv[layer]-compressibility_of_water),dz_clays[layer],unconfined=unconfined_tmp,overburden=overburden_stress_compaction,overburden_data=1/(rho_w * g) * np.array(overburden_data_tmp),endnodes=compaction_solver_debug_include_endnodes)
                else:
                    db[layer]['total_%.2f clays' % thickness],deformation[layer]['total_%.2f clays' % thickness],deformation[layer]['elastic_%.2f clays' % thickness],deformation[layer]['inelastic_%.2f clays' % thickness]=subsidence_solver_aquitard_elasticinelastic(head_series[layer]['%.2f clays' % thickness],inelastic_flag[layer]['%.2f clays' % thickness],(clay_Sse[layer]-compressibility_of_water),(clay_Ssv[layer]-compressibility_of_water),dz_clays[layer],endnodes=compaction_solver_debug_include_endnodes)
                deformation[layer]['total_%.2f clays' % thickness] = interbeds_distributions[layer][thickness] * deformation[layer]['total_%.2f clays' % thickness] 
                deformation[layer]['elastic_%.2f clays' % thickness] = interbeds_distributions[layer][thickness] * deformation[layer]['elastic_%.2f clays' % thickness]
                deformation[layer]['inelastic_%.2f clays' % thickness]= interbeds_distributions[layer][thickness] * deformation[layer]['elastic_%.2f clays' % thickness]
            # Now collect the results at the max timestep
            dt_sand_tmp=np.diff(head_data[layer][:,0])[0]
            print('\tSumming deformation for layer %s. dts are %.2f and %.2f.' % (layer, dt_sand_tmp,dt_master[layer]))
            dt_max_tmp = np.max([dt_sand_tmp,dt_master[layer]])

            t_total_tmp = 0.001 * np.arange(1000*np.min(head_data[layer][:,0]),1000*np.max(head_data[layer][:,0]+0.000001),1000*dt_max_tmp) # this is the master dt which will apply for all the sublayers within layer


            deformation_OUTPUT_tmp={}
            deformation_OUTPUT_tmp['dates'] = [x.strftime('%d-%b-%Y') for x in num2date(t_total_tmp)]

            deformation_OUTPUT_tmp['Interconnected Matrix'] = np.array(deformation[layer]['Interconnected matrix'])[np.where(np.isin(head_data[layer][:,0],t_total_tmp))]
            def_tot_tmp = np.zeros_like(t_total_tmp,dtype='float')               
            
            def_tot_tmp += np.array(deformation[layer]['Interconnected matrix'])[np.where(np.isin(head_data[layer][:,0],t_total_tmp))]
            for thickness in bed_thicknesses_tmp:
                def_tot_tmp += np.array(deformation[layer]['total_%.2f clays' % thickness])[np.isin(0.0001*np.arange(10000*np.min(head_data[layer][:,0]),10000*(np.max(head_data[layer][:,0])+0.0001),10000*dt_master[layer]),t_total_tmp)]
                ting=np.array(deformation[layer]['total_%.2f clays' % thickness])[np.isin(0.0001*np.arange(10000*np.min(head_data[layer][:,0]),10000*(np.max(head_data[layer][:,0])+0.0001),10000*dt_master[layer]),t_total_tmp)]
                deformation_OUTPUT_tmp['total_%.2f clays' % thickness]=ting

            deformation[layer]['total'] = np.array([t_total_tmp,def_tot_tmp])
            deformation_OUTPUT_tmp['total']=def_tot_tmp
            deformation_OUTPUT[layer] = pd.DataFrame(deformation_OUTPUT_tmp)

    if layer_types[layer]=='Aquitard':
        print()
        print('%s is an Aquitard. Solving for layer compaction.' % layer)
        if layer_compaction_switch[layer]:
            if compaction_solver_compressibility_type[layer]=='elastic-inelastic':
                deformation[layer]={}
                
                if overburden_stress_compaction:
                    unconfined_tmp = unconfined_aquifer_name==layer
                    if len(overburden_dates) != len(groundwater_solution_dates[layer]):
                        print('\t\t\tOverburden series is %i long whereas head series is %i long. Interpolating overburden stress.' % (len(overburden_dates),len(groundwater_solution_dates[layer])))
                        f_tmp = scipy.interpolate.interp1d(overburden_dates,overburden_data)
                        overburden_data_tmp = f_tmp(groundwater_solution_dates[layer])
                        overburden_dates_tmp = groundwater_solution_dates[layer]
                    else:
                        overburden_data_tmp = overburden_data

                    db[layer],totdeftmp,deformation[layer]['elastic'],deformation[layer]['inelastic']=subsidence_solver_aquitard_elasticinelastic(head_series[layer],inelastic_flag[layer],(clay_Sse[layer]-compressibility_of_water),(clay_Ssv[layer]-compressibility_of_water),dz_clays[layer],unconfined=unconfined_tmp,overburden=overburden_stress_compaction,overburden_data=1/(rho_w * g) * np.array(overburden_data_tmp))
                    deformation[layer]['total'] = np.array([groundwater_solution_dates[layer],totdeftmp])

                else:
                    db[layer],totdeftmp,deformation[layer]['elastic'],deformation[layer]['inelastic']=subsidence_solver_aquitard_elasticinelastic(head_series[layer],inelastic_flag[layer],(clay_Sse[layer]-compressibility_of_water),(clay_Ssv[layer]-compressibility_of_water),dz_clays[layer])
                    deformation[layer]['total'] = np.array([groundwater_solution_dates[layer],totdeftmp])          
  
solving_compaction_stop = process_time()
solving_compaction_time = solving_compaction_stop - solving_compaction_start

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
            print('Saving figures and data for aquitard layer %s.' % layer)
            print('\tMaking deformation figure')
            sns.set_style('darkgrid')
            sns.set_context('talk')
            plt.figure(figsize=(18,12))
            t = groundwater_solution_dates[layer]
            
            plt.plot_date(t,deformation[layer]['total'][1,:])
            plt.plot_date(t,deformation[layer]['elastic'],label='elastic')
            plt.plot_date(t,deformation[layer]['inelastic'],label='inelastic')
            plt.legend()
            plt.savefig('%s/figures/compaction_%s.png' % (outdestination, layer.replace(' ','_')),bbox_inches='tight')
            plt.xlim([date.toordinal(date(2015,1,1)),date.toordinal(date(2020,1,1))])
            plt.savefig('%s/figures/compaction_%s_20152020.png' % (outdestination, layer.replace(' ','_')),bbox_inches='tight')
            plt.close()
            
            # print('\tSaving db')
            # with open('%s/%s_db.csv' % (outdestination, layer.replace(' ','_')), "w+") as myCsv:
            #     csvWriter = csv.writer(myCsv, delimiter=',')
            #     csvWriter.writerows(np.array(db[layer]).T)
            
            print('\tSaving s_elastic timeseries')
            np.savetxt('%s/%s_s_elastic.csv' % (outdestination, layer.replace(' ','_')),deformation[layer]['elastic'])
    
            print('\tSaving s timeseries')
            np.savetxt('%s/%s_s.csv' % (outdestination, layer.replace(' ','_')),deformation[layer]['total'])
    
            if save_internal_compaction:

                print('\tSaving plots of compaction within layer %s.' % layer)
                t = groundwater_solution_dates[layer]
                x_lims = list(map(dt.fromordinal,[int(min(t)),int(max(t))]))
                x_lims = date2num(x_lims)
                y_lims=[min(Z[layer]),max(Z[layer])]
                
                sns.set_style('white')
                sns.set_context('talk')
                plt.figure(figsize=(18,12))
                plt.imshow(np.array(db[layer]).T,aspect='auto',cmap='RdBu',vmin=-np.max(np.abs(np.array(db[layer])[5:,:])),vmax=np.max(np.abs(np.array(db[layer])[5:,:])),extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]]) # note the min/max are set starting at the 5th timestep because the early timesteps can have large changes due to the initial condition and the boundary condition being discontinuous at these times
                plt.gca().xaxis_date()
                date_format = mdates.DateFormatter('%Y')
                plt.gca().xaxis.set_major_formatter(date_format)
                plt.gcf().autofmt_xdate()
                plt.colorbar(label='db (m)')
                plt.ylabel('Z (m)')
                plt.savefig('%s/figures/%s_compaction_internal.png' % (outdestination, layer.replace(' ','_')),bbox_inches='tight')
                plt.close()
                
                plt.figure(figsize=(18,12))
                plt.imshow(np.array(db[layer]).T,aspect='auto',cmap='RdBu',norm=colors.TwoSlopeNorm(vmin=-np.max(np.abs(np.array(db[layer])[5:,:])), vcenter=-0.1*np.max(np.abs(np.array(db[layer])[5:,:])), vmax=0),extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]]) # note the min/max are set starting at the 5th timestep because the early timesteps can have large changes due to the initial condition and the boundary condition being discontinuous at these times 
                plt.gca().xaxis_date()
                date_format = mdates.DateFormatter('%Y')
                plt.gca().xaxis.set_major_formatter(date_format)
                plt.gcf().autofmt_xdate()
                plt.colorbar(label='db (m)')
                plt.ylabel('Z (m)')
                plt.savefig('%s/figures/%s_compaction_internal_highconstast.png' % (outdestination, layer.replace(' ','_')),bbox_inches='tight')
                plt.close()
#            
    if layer_types[layer]=='Aquifer':
        if layer_compaction_switch[layer]:
            print('Saving figures and data for aquifer layer %s.' % layer)
            if not os.path.isdir('%s/figures/%s' % (outdestination,layer)):
                os.mkdir('%s/figures/%s' % (outdestination,layer))
                           
            print('\tSaving layer data.')
            deformation_OUTPUT[layer].to_csv('%s/%s_Total_Deformation_Out.csv' % (outdestination,layer),index=False)
                
            print('\tSaving layer compaction figure')
            sns.set_style('whitegrid')
            sns.set_context('talk')
            plt.figure(figsize=(18,12))
            plt.plot_date(head_series[layer]['Interconnected matrix'][:,0],deformation[layer]['Interconnected matrix'],'-',label='Interconnected matrix')

            interbeds_tmp=interbeds_distributions[layer]
            bed_thicknesses_tmp=list(interbeds_tmp.keys())

            for thickness in bed_thicknesses_tmp:
                plt.plot_date(groundwater_solution_dates[layer]['%.2f clays' % thickness],deformation[layer]['total_%.2f clays' % thickness],'-',label='%s_%ix%.2f clays' % (layer,interbeds_distributions[layer][thickness],thickness))

            plt.plot_date(deformation[layer]['total'][0,:],deformation[layer]['total'][1,:],'-',label='total')
            plt.title('%s' % layer)
            plt.ylabel('Deformation (m)')
            plt.legend()
            plt.savefig('%s/figures/%s/overall_compaction_%s.png' % (outdestination,layer,layer),bbox='tight')
            plt.xlim([date.toordinal(date(2015,1,1)),date.toordinal(date(2020,1,1))])
            plt.savefig('%s/figures/%s/overall_compaction_%s_201520.png' % (outdestination,layer,layer),bbox='tight')
            plt.close() 
            
            
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

            
            if save_internal_compaction:
                print('\tSaving internal compaction plots for clays of thickness ',bed_thicknesses_tmp)
                for thickness in bed_thicknesses_tmp:
                    print('\t\t',thickness)
                    plt.figure(figsize=(18,12))
                    t = groundwater_solution_dates[layer]['%.2f clays' % thickness]
                    x_lims = list(map(dt.fromordinal,[int(min(t)),int(max(t))]))
                    x_lims = date2num(x_lims)
                    y_lims=[min(Z[layer]['%.2f clays' % thickness]),max(Z[layer]['%.2f clays' % thickness])]
                
                    sns.set_style('white')
                    sns.set_context('talk')
                    plt.figure(figsize=(18,12))
                    plt.imshow(np.array(db[layer]['total_%.2f clays' % thickness]).T,aspect='auto',cmap='RdBu',vmin=-np.max(np.abs(np.array(db[layer]['total_%.2f clays' % thickness])[5:,:])),vmax=np.max(np.abs(np.array(db[layer]['total_%.2f clays' % thickness])[5:,:])),extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]]) # note the min/max are set starting at the 5th timestep because the early timesteps can have large changes due to the initial condition and the boundary condition being discontinuous at these times
                    plt.gca().xaxis_date()
                    date_format = mdates.DateFormatter('%Y')
                    plt.gca().xaxis.set_major_formatter(date_format)
                    plt.gcf().autofmt_xdate()
                    plt.colorbar(label='db (m)')
                    plt.ylabel('Z (m)')
                    plt.savefig('%s/figures/%s/%sclay_compaction_internal.png' % (outdestination, layer,thickness),bbox_inches='tight')
                    plt.close()
                    plt.figure(figsize=(18,12))
                    plt.imshow(np.array(db[layer]['total_%.2f clays' % thickness]).T,aspect='auto',cmap='RdBu',norm=colors.TwoSlopeNorm(vmin=-np.max(np.abs(np.array(db[layer]['total_%.2f clays' % thickness])[5:,:])), vcenter=-0.1*np.max(np.abs(np.array(db[layer]['total_%.2f clays' % thickness])[5:,:])), vmax=0)
    ,extent = [x_lims[0], x_lims[1],  y_lims[0], y_lims[1]]) # note the min/max are set starting at the 5th timestep because the early timesteps can have large changes due to the initial condition and the boundary condition being discontinuous at these times
                    plt.gca().xaxis_date()
                    date_format = mdates.DateFormatter('%Y')
                    plt.gca().xaxis.set_major_formatter(date_format)
                    plt.gcf().autofmt_xdate()
                    plt.colorbar(label='db (m)')
                    plt.ylabel('Z (m)')
                    plt.savefig('%s/figures/%s/%sclay_compaction_internalHIGCONTRAST.png' % (outdestination, layer,thickness),bbox_inches='tight')
                    plt.close()

print("Creating overall compaction plot and saving deformation series")
sns.set_style('whitegrid')
plt.figure(figsize=(18,12))
dt_master_compacting_layers = {key:dt_master[key] for key in compacting_layers}
maxdt = max(dt_master_compacting_layers.values())
maxdtlayer = max(dt_master_compacting_layers.items(), key=operator.itemgetter(1))[0]
dt_interconnecteds = [dt_headseries[layer] for layer in layer_names if (layer in compacting_layers) and (layer_types[layer]=='Aquifer')]
if np.max(dt_interconnecteds)>maxdt:
    maxdt=np.max(dt_interconnecteds)
deformation_OUTPUT = {}
t_total_tmp = 0.0001*np.arange(10000*np.min(deformation[maxdtlayer]['total'][0,:]),10000*(np.max(deformation[maxdtlayer]['total'][0,:]+0.001)),10000*maxdt)
deformation_OUTPUT['dates']=[x.strftime('%d-%b-%Y') for x in num2date(t_total_tmp)]
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
            deformation_OUTPUT[layer]=deformation[layer]['total'][1,:][np.where(np.isin(deformation[layer]['total'][0,:],t_total_tmp))]


        if layer_types[layer]=='Aquifer':
            l_tmp, = plt.plot_date(deformation[layer]['total'][0,:],deformation[layer]['total'][1,:],'-',label='%s' % layer)
            l_aqf.append(l_tmp)
            #dt_tmp = int(maxdt/dt_master[layer])
            deformation_OUTPUT[layer]=deformation[layer]['total'][1,:][np.where(np.isin(deformation[layer]['total'][0,:],t_total_tmp))]

# Add up all the deformations from each layer
            
def_tot_tmp = np.zeros_like(t_total_tmp,dtype='float')
for layer in layer_names:
    if layer_compaction_switch[layer]:
        newtot = np.array(deformation[layer]['total'][1,:])[np.isin(deformation[layer]['total'][0,:],t_total_tmp)]
        t_overall = t_overall + newtot

deformation_OUTPUT['Total']=t_overall
#print(deformation_OUTPUT)
def_out = pd.DataFrame(deformation_OUTPUT)
def_out.to_csv('%s/Total_Deformation_Out.csv' % outdestination,index=False)

l3, = plt.plot_date(t_total_tmp,t_overall,label='TOTAL def')

plt.ylabel('Z (m)')
plt.legend()
plt.savefig('%s/figures/total_deformation_figure.png' % outdestination,bbox_inches='tight')
# Rezero on jan 2015
plt.xlim([date.toordinal(date(2015,1,1)),date.toordinal(date(2020,1,1))])
for line in l_aqt:
    line.set_ydata(np.array(line.get_ydata()) - np.array(line.get_ydata())[np.array(line.get_xdata())==date2num(date(2015,1,1))])
plt.savefig('%s/figures/total_deformation_figure_20152020.png' % outdestination,bbox_inches='tight')
plt.close()

saving_compaction_stop = process_time()
saving_compaction_time = saving_compaction_stop - saving_compaction_start

t_total_stop = process_time()
t_total = t_total_stop - t_total_start

plt.figure(figsize=(18,18))
plt.pie(np.abs([param_read_time, solving_compaction_time, saving_head_time, reading_head_time, solving_head_time, saving_compaction_time,t_total - np.sum([param_read_time, solving_compaction_time, saving_head_time, reading_head_time, solving_head_time, saving_compaction_time])]),labels=['param_read_time', 'solving_compaction_time', 'saving_head_time', 'reading_head_time', 'solving_head_time', 'saving_compaction_time','misc'],autopct=lambda p : '{:.2f}%  ({:,.0f})'.format(p,p * t_total/100))
plt.title('Total runtime = %i seconds' % t_total)
plt.savefig('%s/figures/runtime_breakdown.png' % outdestination,bbox_inches='tight')
plt.close()
#for line in l_aqf:
#    line.set_ydata(np.array(line.get_ydata()) - np.array(line.get_ydata())[np.array(line.get_xdata())==date.toordinal(date(2015,1,1))])
#
#l3.set_ydata(np.array(l3.get_ydata()) - np.array(l3.get_ydata())[np.array(l3.get_xdata())==date.toordinal(date(2015,1,1))])
#
#plt.ylim([2*np.array(l3.get_ydata())[np.array(l3.get_xdata())==date.toordinal(date(2020,1,1))],-0.2*np.array(l3.get_ydata())[np.array(l3.get_xdata())==date.toordinal(date(2020,1,1))]])
#plt.savefig('%s/figures/total_deformation_figure_20152020zoom.png' % outdestination,bbox_inches='tight')


#print("Creating overall compaction % plot")
#overall_pcs={}
#
#
#for layer in layer_names:
#    if layer_compaction_switch[layer]:
#        dt_tmp = int(maxdt/dt_master[layer])
#        overall_pcs[layer] = 100*np.array([(deformation[layer]['total'][1,::dt_tmp][i+1] - deformation[layer]['total'][1,::dt_tmp][i])/(t_overall[i+1] - t_overall[i]) for i in range(len(t_overall)-1)])
#sns.set_style('white')
#
#plt.figure(figsize=(18,12))
#for layer in layer_names:
#    if layer_compaction_switch[layer]:
#        plt.plot_date(date2num(np.array([dt.strptime(date, '%d-%b-%Y').date() for date in groundwater_solution_dates[maxdtlayer]]))[:-1][np.abs(overall_pcs[layer])<=200],np.array(overall_pcs[layer])[np.abs(overall_pcs[layer])<=200],'-',label='%s' % layer)
#
#ax1 = plt.gca()
#plt.ylabel('instantaneous % deformation')
#ax2 = ax1.twinx()
#line, = ax2.plot_date(date2num([dt.strptime(date, '%d-%b-%Y').date() for date in groundwater_solution_dates[maxdtlayer]]),t_overall,'k--',label='TOTAL def')
#
#lines, labels = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax2.legend(lines + lines2, labels + labels2,fancybox=True)
#plt.savefig('%s/figures/total_deformation_percentage_figure.png' % outdestination,bbox_inches='tight')
#
#plt.xlim([date.toordinal(date(2015,1,1)),date.toordinal(date(2020,1,1))])
## Rezero on jan 2015
#line.set_ydata(np.array(line.get_ydata()) - np.array(line.get_ydata())[np.array(line.get_xdata())==date.toordinal(date(2015,1,1))])
#plt.ylim((1.2*np.min(np.array(line.get_ydata())[np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2015,1,1)))[0][0]:np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2020,1,1)))[0][0]]),1.2*np.max(np.array(line.get_ydata())[np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2015,1,1)))[0][0]:np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2020,1,1)))[0][0]])))
#plt.ylabel('Total Deformation (m)')
#plt.savefig('%s/figures/total_deformation_percentage_figure_20152020zoom.png' % outdestination,bbox_inches='tight')
#
## Now repeat with 31 day smoothing
#plt.figure(figsize=(18,12))
#for layer in layer_names:
#    if layer_compaction_switch[layer]:
#        plt.plot_date(date2num(np.array([dt.strptime(date, '%d-%b-%Y').date() for date in groundwater_solution_dates[maxdtlayer]]))[:-1][np.abs(overall_pcs[layer])<=200],convolve(np.array(overall_pcs[layer])[np.abs(overall_pcs[layer])<=200],np.ones((31,))/31),'-',label='%s' % layer)
#
#ax1 = plt.gca()
#plt.ylabel('instantaneous % deformation')
#ax2 = ax1.twinx()
#line, = ax2.plot_date(date2num([dt.strptime(date, '%d-%b-%Y').date() for date in groundwater_solution_dates[maxdtlayer]]),t_overall,'k--',label='TOTAL def')
#
#lines, labels = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax2.legend(lines + lines2, labels + labels2,fancybox=True)
#
#plt.xlim([date.toordinal(date(2015,1,1)),date.toordinal(date(2020,1,1))])
## Rezero on jan 2015
#line.set_ydata(np.array(line.get_ydata()) - np.array(line.get_ydata())[np.array(line.get_xdata())==date.toordinal(date(2015,1,1))])
#plt.ylim((1.2*np.min(np.array(line.get_ydata())[np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2015,1,1)))[0][0]:np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2020,1,1)))[0][0]]),1.2*np.max(np.array(line.get_ydata())[np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2015,1,1)))[0][0]:np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2020,1,1)))[0][0]])))
#plt.ylabel('Total Deformation (m)')
#plt.savefig('%s/figures/total_deformation_percentage_figure_20152020zoom_31daysmooth.png' % outdestination,bbox_inches='tight')
#
#
## Finally do it cumulative
#
#for layer in layer_names:
#    if layer_compaction_switch[layer]:
#        dt_tmp = int(maxdt/dt_master[layer])
#        overall_pcs[layer] = 100*np.array([(deformation[layer]['total'][1,::dt_tmp][i] - deformation[layer]['total'][1,::dt_tmp][0])/(t_overall[i] - t_overall[0]) for i in range(len(t_overall))])
#sns.set_style('white')
#
#plt.figure(figsize=(18,12))
#for layer in layer_names:
#    if layer_compaction_switch[layer]:
#        plt.plot_date(date2num(np.array([dt.strptime(date, '%d-%b-%Y').date() for date in groundwater_solution_dates[maxdtlayer]]))[:],np.array(overall_pcs[layer]),'-',label='%s' % layer)
#
#plt.ylim([-10,110])
#ax1 = plt.gca()
#plt.ylabel('cumulative % deformation')
#ax2 = ax1.twinx()
#line, = ax2.plot_date(date2num([dt.strptime(date, '%d-%b-%Y').date() for date in groundwater_solution_dates[maxdtlayer]]),t_overall,'k--',label='TOTAL def')
#
#lines, labels = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax2.legend(lines + lines2, labels + labels2,fancybox=True)
#plt.savefig('%s/figures/total_deformation_cum_percentage_figure.png' % outdestination,bbox_inches='tight')
#
## Rezero on jan 2015
#for layer in layer_names:
#    if layer_compaction_switch[layer]:
#        dt_tmp = int(maxdt/dt_master[layer])
#        arg2015 = np.argwhere(date2num(np.array([dt.strptime(date, '%d-%b-%Y').date() for date in groundwater_solution_dates[maxdtlayer]]))==date.toordinal(date(2015,1,1)))[0][0]
#        overall_pcs[layer] = 100*np.array([(deformation[layer]['total'][1,::dt_tmp][i] - deformation[layer]['total'][1,::dt_tmp][arg2015])/(t_overall[i] - t_overall[arg2015]) for i in range(len(t_overall))])
#
#
#plt.figure(figsize=(18,12))
#for layer in layer_names:
#    if layer_compaction_switch[layer]:
#        plt.plot_date(date2num(np.array([dt.strptime(date, '%d-%b-%Y').date() for date in groundwater_solution_dates[maxdtlayer]]))[:],np.array(overall_pcs[layer]),'-',label='%s' % layer)
#
#plt.ylim([-10,110])
#ax1 = plt.gca()
#plt.ylabel('cumulative % deformation')
#ax2 = ax1.twinx()
#line, = ax2.plot_date(date2num([dt.strptime(date, '%d-%b-%Y').date() for date in groundwater_solution_dates[maxdtlayer]]),t_overall,'k--',label='TOTAL def')
#
#lines, labels = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax2.legend(lines + lines2, labels + labels2,fancybox=True)
#plt.savefig('%s/figures/total_deformation_cum_percentage_figure.png' % outdestination,bbox_inches='tight')
#
#
#line.set_ydata(np.array(line.get_ydata()) - np.array(line.get_ydata())[np.array(line.get_xdata())==date.toordinal(date(2015,1,1))])
#plt.xlim([date.toordinal(date(2015,1,1)),date.toordinal(date(2020,1,1))])
#plt.ylim((1.2*np.min(np.array(line.get_ydata())[np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2015,1,1)))[0][0]:np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2020,1,1)))[0][0]]),1.2*np.max(np.array(line.get_ydata())[np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2015,1,1)))[0][0]:np.argwhere(np.array(line.get_xdata())==date.toordinal(date(2020,1,1)))[0][0]])))
#plt.ylabel('Total Deformation (m)')
#plt.savefig('%s/figures/total_deformation_cum_percentage_figure_20152020zoom.png' % outdestination,bbox_inches='tight')



