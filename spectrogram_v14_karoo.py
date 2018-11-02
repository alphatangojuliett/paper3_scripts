#!/usr/bin/python
from fft_accum_avg_p1_rts_2048 import fft_accum_avg_p1_rts_2048
import matplotlib.pyplot as plt
import numpy as np
import scipy
import json 
import io
import datetime 
import argparse

parser = argparse.ArgumentParser(description='Adjust spectrogram sweep settings and outputs.')

parser.add_argument('--out_data', help='Boolean: Do you want to output 2D array of data for specpeak.py?', default = 1)
parser.add_argument('--out_ridm', help = 'Boolean: Do you want to output a .ridm file of instrument/measurement metadata?', default = 0)
parser.add_argument('--cal_internal', help = 'Boolean: Do you want your calibration to avoid the effects of FFT power when nothing is connected to your SDR?', default=1)
parser.add_argument('--samp_rate', help='Sample Rate (MHz)', default=10.0)
parser.add_argument('--fmin', help = 'The minimum value (MHz) of the SDR center frequency (f_c). Frequencies to be measured and recorded per spectra are are in range f_c +/ (SAMP_RATE - BW_OVERLAP/2) ', default =29. )
parser.add_argument('--fmax', help = 'The maximum value (MHz) of the SDR center frequency (f_c). Frequencies to be measured and recorded per spectra are are in range f_c +/ (SAMP_RATE - BW_OVERLAP/2) ', default =499. )
parser.add_argument('--overlap',help = 'The amount of bandwidth (MHz) to be ignored during each FFT, to avoid the effects of analog filtes and DSP at band edges. BW_OVERLAP/2 is ignored from each edge of FFT', default =0)
parser.add_argument('--dt',help = 'Exectute a new sweep once every (this many) seconds. Should not be < 17s if integration time is 0.2s', default = 30)
parser.add_argument('--t_total',help = 'End this script (stop sweeping) after this many seconds.',default = 1 )
parser.add_argument('--ant', help = 'This is the name of the antenna connected to the SDR',default=str('Bicolog 30100'))

args = parser.parse_args()

#GLOBAL PARAMS
REMOVE_CLOCK_AND_SIDE_CHANNELS = False #A very crude way to eliminate certain channels from contaminating data. Only to be used for troubleshooting. Do not use when recording data.
CALIBRATE_OUT_SDR_NO_CONNECTION = bool(int(args.cal_internal)) # A measurement calibration, accounting for fact that the SDR has nonzero FFT power even when nothing is plugged in. 
REMOVE_CLOCK_SPURIOUS = False #Deletes data at various channels known to cause spurs. Should only be used while troubleshooting, data should not saved while this is True.
OUTPUT_DATA = bool(int(args.out_data))
OUTPUT_RIDM = bool(int(args.out_ridm))
CAL_EACH_SPECTRA = True # Every SAMP_RATE-wide spectra of the sweep must have its own absolutely-calibrated point. Must be on to make accurate masurement in that spectrum. Turn off only if troubleshooting with fewer cal points.
PLOT_EACH_SWEEP = False #Should only be used as a diagnostic. 

SAMP_RATE = float(float(args.samp_rate)*1.e6) 
SIZE_FFT = 2048 # Number of channels in each FFT. If this is changed, then 'vec_FFT_len' must be manually (couldn't make it work automatically) changed in the gnuradio module.
FMIN = float(float(args.fmin)*1.e6)
FMAX = float(float(args.fmax)*1.e6)
BW_OVERLAP = float(float(args.overlap)*1.e6)
SPURIOUS_CHAN_TO_DEL = 11 # Used only when troubleshooting via REMOVE_CLOCK_AND_SIDE_CHANNELS. The number channels, center on spur, that should be ignored. Must be odd number. 
MOD_SEC_SCHEDULE = float(args.dt) # A new sweep spectrum is triggered every (start_time = datetime.datetime.now) % MOD_SEC_SCHEDULE 
DT_EXIT = float(args.t_total) # End this script after this many seconds.
T_INT= 0.2 #If this changes, you have to recalibrate the SDR. Keep at 0.2 for now. 

n_to_avg = int(np.round((SAMP_RATE/SIZE_FFT)*T_INT))# *** should this have int() around it?

STR_ANT_NAME = str(args.ant)

STR_FILE_TEMP = '/home/hera/SDR_RFI_monitoring/fft_accum_temp/fft_accum_temp.txt'
STR_FILE_MIN = '/home/hera/SDR_RFI_monitoring/fft_accum_temp/fft_accum_min.txt'
STR_FILE_MAX = '/home/hera/SDR_RFI_monitoring/fft_accum_temp/fft_accum_max.txt'
STR_FILE_MEAN = '/home/hera/SDR_RFI_monitoring/fft_accum_temp/fft_accum_mean.txt'
STR_FILE_CAL_NOISE = '/home/hera/SDR_RFI_monitoring/25-500MHz_2kHzRBW_500points_noise_source.csv'#noise_source_cal_patch_2_files_and_running_mean_1st.csv'#noise_source_cal_20kHz_RBW_250_points_prep.csv'#spectrum analysis of noise source, measured by keysight
STR_FILE_CAL_SDR = '/home/hera/SDR_RFI_monitoring/airspyR2_noise_source_dot2s_10MSPS_29-499MHz_0MHZ_ovlap.txt'#HackRF_noise_source_gains_0_20_20_BW_20_sizefft_2048.txt'#These are the SDR amplitude (=sqrt power) values when the noise source is connected, measured via spectrogram_v4_working
STR_FILE_SDR_NO_CONNECTION = '/home/hera/SDR_RFI_monitoring/airspy_nothing_connected.txt'#SDR_nothing_connected_1s.txt'
STR_CALIBRATION = '/home/hera/SDR_RFI_monitoring/airspy_calibration_v1.csv'#latest_calibration_patched_mean.csv'



arr_latest_cal = np.loadtxt(STR_CALIBRATION,delimiter=',')


if (np.array(arr_latest_cal).ndim > 1):
	ARR_F_CAL = arr_latest_cal[:,1]*1.e6
	ARR_P_CAL_SDR = arr_latest_cal[:,2]
	ARR_P_CAL_VNA = arr_latest_cal[:,0]
else:
	ARR_F_CAL = arr_latest_cal[1]
	ARR_P_CAL_SDR = arr_latest_cal[2]
	ARR_P_CAL_VNA = arr_latest_cal[0]

n_FFTs = int(((FMAX-FMIN)/(SAMP_RATE- BW_OVERLAP) +1))
n_spec_cal = np.floor( (ARR_F_CAL - (FMIN - (SAMP_RATE- BW_OVERLAP)/2)) / (SAMP_RATE - BW_OVERLAP) ).astype(int)
#print('FMIN = '+str(FMIN))
#print('SAMP_RATE = '+str(SAMP_RATE))
#print('BW_OVERLAP = '+str(BW_OVERLAP))
#print('ARR_F_CAL: '+str(ARR_F_CAL))
#print('n_spec_cal: '+str(n_spec_cal))
SIDE_CHAN_TO_DEL=0
CENTER_CHAN_TO_DEL=2

ii_first_keep = int( ((BW_OVERLAP/2)/SAMP_RATE)*SIZE_FFT )
ii_last_keep = int( (1-(BW_OVERLAP/2)/SAMP_RATE)*SIZE_FFT )


###### Functions
def remove_clock_effects(arr, freqs):
	if (np.size(arr) > 0):
		spurious_mask = np.arange(10.0*1.e6, (FMAX+(SAMP_RATE - BW_OVERLAP)/2), 5.0*1.e6)
		spurious_indices = np.searchsorted(freqs,spurious_mask)
		for ii in spurious_indices: #index of the integer multiples of 5MHz
			if (ii == 0): spurious_indices = np.append(spurious_indices,np.arange(ii, ii + SPURIOUS_CHAN_TO_DEL/2 +1,1))
			spurious_indices = np.append(spurious_indices, np.arange(ii - SPURIOUS_CHAN_TO_DEL/2, ii + SPURIOUS_CHAN_TO_DEL + 1,1))
		spurious_indices=np.unique(spurious_indices)
		arr[spurious_indices] = np.nan
	return arr

def find_index_of_closest(arr, val):
    arr = np.asarray(arr)
    ii = (np.abs(arr - val)).argmin()
    return ii
    #return arr[ii]
def find_index_of_closest_given_array(arr, values):
    arr = np.asarray(arr)
    ii=0
    arr_indices = np.zeros(np.size(values))
    for val in values:
    	arr_indices[int(ii)] = (np.abs(arr - val)).argmin()
    	ii=ii+1
    return arr_indices

def find_nearest(arr, val):
    arr = np.asarray(arr)
    ii = (np.abs(arr - val)).argmin()
    return arr[ii]


def analysis_of_sample():
	file = scipy.fromfile(open(STR_FILE_TEMP), dtype=scipy.float32)
	dat1D = np.array(file)
	#print('size of 1D is: '+str(np.size(dat1D)))
	dat2D = np.reshape( file,(np.size(file)/SIZE_FFT,SIZE_FFT) ) #order them by every size_fft cols to get rid of edge and clock peaks
	if REMOVE_CLOCK_AND_SIDE_CHANNELS:
		elims = np.arange(int(SIZE_FFT/2 - CENTER_CHAN_TO_DEL), int(SIZE_FFT/2 +CENTER_CHAN_TO_DEL +1),1) #elimiate the clock peaks and the N channels on each side
		elims = np.concatenate((elims, np.arange(0,int(SIDE_CHAN_TO_DEL)),np.arange(int(SIZE_FFT - SIDE_CHAN_TO_DEL),int(SIZE_FFT))))
		elims = np.unique(elims)

		for ii in elims: #set the power at each clock peak and edge (+/- side_chan_to_del channels) to zero
			dat2D[:,ii]=0

	n_to_avg_effective = int(np.size(dat1D)/2048) # this should be exactly the same as n_to_avg if the SDR works, otherwise the SDR may not be recording FFTs despite recording a head output of length vec_fft_len
	#print('n_to_avg_effective = '+str(n_to_avg_effective))
	#link for following line: https://stackoverflow.com/questions/30379311/fast-way-to-take-average-of-every-n-rows-in-a-npy-array
	#note, the power(,2) below is to convert from measured amplitude to measured power. Confirmed in by doing input power tests.
	dat2D_mean = dat2D.transpose().reshape(-1,n_to_avg_effective).mean(1).reshape(SIZE_FFT,-1).transpose()# This line avgs every R rows together in a 2D array with C cols
	dat2D_min = dat2D.transpose().reshape(-1,n_to_avg_effective).min(1).reshape(SIZE_FFT,-1).transpose()
	dat2D_max = dat2D.transpose().reshape(-1,n_to_avg_effective).max(1).reshape(SIZE_FFT,-1).transpose()
	#print('Shape of data after transposes: '+str(np.shape(dat2D_mean)))
	f_handle_mean = open(STR_FILE_MEAN,'ab')
	f_handle_min = open(STR_FILE_MIN,'ab')
	f_handle_max = open(STR_FILE_MAX,'ab')
	np.savetxt(f_handle_mean,np.array(np.ravel(dat2D_mean))[ii_first_keep:ii_last_keep])#,dtype=scipy.float32)
	np.savetxt(f_handle_min,np.array(np.ravel(dat2D_min))[ii_first_keep:ii_last_keep])#,dtype=scipy.float32)
	np.savetxt(f_handle_max,np.array(np.ravel(dat2D_max))[ii_first_keep:ii_last_keep])#,dtype=scipy.float32)
	f_handle_mean.close()
	f_handle_min.close()
	f_handle_max.close()
	open(STR_FILE_TEMP, 'w').close() #delete contents of temporary file

def file_flush():
	print('flushing temporary and output files...')
	open(STR_FILE_TEMP, 'w').close()
	open(STR_FILE_MEAN, 'w').close()
	open(STR_FILE_MIN, 'w').close()
	open(STR_FILE_MAX, 'w').close()
	print('done.')

def calibrate_noise_source(freqs):
	print('Calibrating noise source (one-time)...')
	SDR_no_connection = np.power(np.loadtxt(STR_FILE_SDR_NO_CONNECTION),2)
	SDR_noise_source = np.power(np.loadtxt(STR_FILE_CAL_SDR),2)
	
	noise_source_dbm_ref=np.loadtxt(STR_FILE_CAL_NOISE, delimiter=',')#MAY HAVE DIFFERENT NUMBER OF ELEMENTS THAN sdr_noise_source


	#add calibration arrays here
	arr_cal = np.zeros((6,np.size(freqs)))# rows, 0: freqs, 1: f_cal, 2: i_cal, 3: P_CAL_VNA, 4: P_CAL_SDR, 5: SDR_noise_source 
	arr_cal[0,:] = freqs
	arr_ii_of_length_SDR = find_index_of_closest_given_array(noise_source_dbm_ref[:,0],freqs)
	for ii in np.arange(np.size(freqs)):
		n_spec_freq = np.floor( (freqs[ii] - (FMIN - (SAMP_RATE- BW_OVERLAP)/2)) / (SAMP_RATE - BW_OVERLAP) ).astype(int)#determine which of the n independent spectra this value belongs to
		if (n_spec_freq == int(n_FFTs)): n_spec_freq = int(n_spec_freq-1) # This stops last index in freqs from throwing error
		if CAL_EACH_SPECTRA: arr_cal[1,ii] = find_nearest(ARR_F_CAL[np.asarray(n_spec_cal == n_spec_freq)],freqs[ii])
		else: arr_cal[1,ii] = find_nearest(ARR_F_CAL,freqs[ii])
		arr_cal[2,ii] = int( find_index_of_closest(freqs,arr_cal[1,ii]) )
		i_P_cal =  int(np.where(ARR_F_CAL==arr_cal[1,ii])[0])
		arr_cal[3,ii] = ARR_P_CAL_VNA[i_P_cal]
		arr_cal[4,ii] = ARR_P_CAL_SDR[i_P_cal]
		arr_cal[5,ii] = SDR_noise_source[int(arr_cal[2,ii])]

	#calibrate-out the SDR behaviour when nothing is connected, but only measurement power > no_connection power (to avoid subtracting noise floor):
	if CALIBRATE_OUT_SDR_NO_CONNECTION:
		SDR_noise_source[SDR_noise_source>SDR_no_connection]=np.array(SDR_noise_source - SDR_no_connection)[SDR_noise_source>SDR_no_connection]
		np.array(arr_cal[4,:])[np.array(arr_cal[4,:]) > SDR_no_connection] = np.array( arr_cal[4,:] - np.take(SDR_no_connection,np.array(arr_cal[2,:]).astype(int)) )[np.array(arr_cal[4,:]) > SDR_no_connection]
	


	cal_noise_source_linear = SDR_noise_source- (arr_cal[5,:]-arr_cal[4,:])
	print('Done.')
	np.save('arr_ii_of_length_SDR_airspy.npy',arr_ii_of_length_SDR)
	np.save('noise_source_dbm_ref_airspy.npy',noise_source_dbm_ref)
	np.save('cal_noise_source_linear_airspy.npy',cal_noise_source_linear)
	np.save('SDR_no_connection_airspy.npy',SDR_no_connection)
	return arr_ii_of_length_SDR, noise_source_dbm_ref, cal_noise_source_linear, SDR_no_connection

def calibrate_spectra(freqs,SDR_meas_uncal, arr_ii_of_length_SDR,noise_source_dbm_ref,cal_noise_source_linear,SDR_no_connection ):

	if CALIBRATE_OUT_SDR_NO_CONNECTION: #Only calibrate on channels where measured power > disconnected power. 
		SDR_meas_uncal[SDR_meas_uncal>SDR_no_connection] = np.array(SDR_meas_uncal - SDR_no_connection)[SDR_meas_uncal>SDR_no_connection]
	
	abs_diff = np.abs(SDR_meas_uncal - cal_noise_source_linear)
	power_ratio_for_dbm_cal = np.divide(SDR_meas_uncal, cal_noise_source_linear) # usual case: big peak for measurement, which is above measured noise. 
	mask_negative_val = power_ratio_for_dbm_cal <0
	power_ratio_for_dbm_cal[mask_negative_val] =  np.array(np.divide((SDR_meas_uncal+abs_diff),SDR_meas_uncal))[mask_negative_val]
	num_mask_negative_val = -1*mask_negative_val # You either want to add or subtract in dBm, this is how I'm decided which operation.
	num_mask_negative_val[num_mask_negative_val == 0] = 1


	dbm_adjustment_to_meas = num_mask_negative_val*10.*np.log10(power_ratio_for_dbm_cal) 

	calibrated_SDR_meas_dbm = np.take(noise_source_dbm_ref[:,1], arr_ii_of_length_SDR.astype(int) ) + dbm_adjustment_to_meas

	return calibrated_SDR_meas_dbm



def load_spectra_and_calibrate(cal_dict,freqs,arr_ii_of_length_SDR,noise_source_dbm_ref,cal_noise_source_linear,SDR_no_connection): # 
	print('Calibrating data...')
	arr_min = np.power(np.loadtxt(STR_FILE_MIN),2)
	arr_max = np.power(np.loadtxt(STR_FILE_MAX),2)
	arr_mean = np.power(np.loadtxt(STR_FILE_MEAN),2)

	noise_source_SDR = np.power(np.loadtxt(STR_FILE_CAL_SDR),2)

	noise_source_dbm_ref=np.loadtxt(STR_FILE_CAL_NOISE, delimiter=',')

	
	cal_dict['val'] = calibrate_spectra(freqs, arr_mean,arr_ii_of_length_SDR,noise_source_dbm_ref,cal_noise_source_linear,SDR_no_connection)
	cal_dict['minhold'] = calibrate_spectra(freqs, arr_min,arr_ii_of_length_SDR,noise_source_dbm_ref,cal_noise_source_linear,SDR_no_connection)
	cal_dict['maxhold']= calibrate_spectra(freqs, arr_max,arr_ii_of_length_SDR,noise_source_dbm_ref,cal_noise_source_linear,SDR_no_connection)


	if REMOVE_CLOCK_SPURIOUS: 
		print(' Removing clock effects from all arrays in cal_dict...')
		for key in cal_dict:
			cal_dict[key]= remove_clock_effects(cal_dict[key],freqs) # you can pass array of size zero, if you do, nothing is done to the array.
		print( ' ...done')
	print('Data is calibrated.')
	return cal_dict

def plot_spectra(freqs, arr):

	plt.plot(freqs/1.e6, arr,'k', label = 'Average Calibrated Measurement')

	plt.legend()
	plt.xlabel('F (MHz)')
	plt.ylabel('Power (dBm)')
	plt.show()

def time_in_rids_fmt(datetime_time): # convert datetime.datetime.now() time into the RIDZ-requested format
	str_iso = datetime_time.isoformat(' ')
	str_time_rids = str( str_iso[0:4] + str_iso[5:7] + str_iso[8:10] + '-' + str_iso[11:13]+str_iso[14:16]+str_iso[17:19])
	return str_time_rids

def output_ridm_metadata(initial_time):
	print('Output ridm metadata...')
	metadata_dict = {'instrument':'AirspyR2 SDR','antenna':str(STR_ANT_NAME),'comment':'AirspyR2 is in container, running via paper3','time_format':'%Y%m%d-%H%M%S','channel_width':str(str((SAMP_RATE/SIZE_FFT)/1.e6)+' MHz'),'freq_unit':'MHz','time_constant':
	str(str(T_INT)+'s per spectra, '+str(n_FFTs)+'s total.') } #add all other keys you want, EXCEPT feature_sets. see github rids readme for accepted strings. 
	str_init_time= time_in_rids_fmt(initial_time)
	str_ridm_out = str('rids_output/sweep_spectrum.'+str(str_init_time)+'.ridm')
	#For Python 2/3 compatability.
	try:
		to_unicode = unicode
	except NameError:
		to_unicode = str

	#write the rids file (JSON format):
	with io.open(str_ridm_out, 'w', encoding='utf8') as outfile:
		str_ = json.dumps(metadata_dict,indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)
		outfile.write(to_unicode(str_))

	print('Done.')


def output_2D_data(start_time,freq,cal_dict): #Output a 2D data array of frequency and FFT power data (corresponding to the cal_dict keys of 'minhold', 'maxhold', and 'val')
	str_start_time = time_in_rids_fmt(start_time)
	for key in cal_dict:
		str_file_out = str('rids_output/SDR.'+str(str_start_time)+'.'+str(key)+'.unk')
		arr_2D_out = np.column_stack((freq,cal_dict[key]))
		np.savetxt(str_file_out,arr_2D_out)

# Main
freqs_fc = np.arange(FMIN,FMAX + (SAMP_RATE- BW_OVERLAP), SAMP_RATE- BW_OVERLAP) # The frequencies to which I program the SDR
freqs=np.linspace((FMIN-(SAMP_RATE- BW_OVERLAP)/2),(FMAX+(SAMP_RATE - BW_OVERLAP)/2),int(ii_last_keep-ii_first_keep)*n_FFTs) # the frequencies of every channel in my FFT
#arr_ii_of_length_SDR, noise_source_dbm_ref, cal_noise_source_linear, SDR_no_connection = calibrate_noise_source(freqs) # All calibration that isn't measurement specific is done here.
arr_ii_of_length_SDR=np.load('arr_ii_of_length_SDR_airspy.npy')
SDR_no_connection=np.load('SDR_no_connection_airspy.npy')
noise_source_dbm_ref = np.load('noise_source_dbm_ref_airspy.npy')
cal_noise_source_linear=np.load('cal_noise_source_linear_airspy.npy')
#print('Number of integrated spectra to be created per sweep: '+str(n_FFTs)) #SCH
#print('Per integrated spectra, integration time: '+str( np.round(float(n_to_avg)/ (SAMP_RATE/SIZE_FFT),decimals=1) )+'s, corresponding to '+str(n_to_avg)+' accumulated FFTs') #SCH

# since the program is not being executed by a chron job, and runs for a given amount of time, these bools guide the sweep schedling.
bool_continue_schedule = True #sch
bool_wait_min_time = False #sch
bool_sweep_has_happened = False
initial_time  = datetime.datetime.now()
cal_dict = {}
if OUTPUT_RIDM: output_ridm_metadata(initial_time)
print('Waiting for first scheduled sweep...') #SCH
while bool_continue_schedule:  # SCH
	start_time = datetime.datetime.now() #SCH
	if (( (start_time- initial_time).seconds >= DT_EXIT ) and (bool_sweep_has_happened)):
			#print('... Max scheduled time has passed. Execute no more sweeps.')
			break
	if ((start_time.second % MOD_SEC_SCHEDULE == 0) and (bool_wait_min_time == False)) : 
		#print('Start time: '+str(start_time)) # SCH,
		cal_dict.clear()
		bool_sweep_has_happened = True

		file_flush()
		#print('Sweeping...')
		FFT_to_file = fft_accum_avg_p1_rts_2048()
		FFT_to_file.set_str_file_out(STR_FILE_TEMP)
		FFT_to_file.set_samp_rate(SAMP_RATE)
		FFT_to_file.set_vec_fft_len(SIZE_FFT)

		start_time = datetime.datetime.now() # This re-definition might break something, but goal is to know time over which the sweep occured.
		for ff in freqs_fc:
			FFT_to_file.set_n_head(n_to_avg)
			#print('setting f_c to '+str(ff/1.e6)+' MHz...')
			FFT_to_file.set_f_c(ff)
			#print('f_c set to: '+str( (FFT_to_file.get_f_c()/1.e6) ) +' MHz')
			FFT_to_file.start()
			FFT_to_file.wait()
			FFT_to_file.stop()
			FFT_to_file.blocks_head_0.reset() #This resets the 'head' block counter, so that I can accumulate n_to_avg spectra at *every* f_c, not just one.
			analysis_of_sample()
	
			
		FFT_to_file.stop()
		del(FFT_to_file)


		end_time = datetime.datetime.now() #sch	
		#print('Finished sweep at '+str(end_time)) #SCH
		cal_dict =  load_spectra_and_calibrate(cal_dict,freqs,arr_ii_of_length_SDR, noise_source_dbm_ref, cal_noise_source_linear, SDR_no_connection)
		dt_init = end_time - initial_time #sch
		dt_start = end_time - start_time

		if ( ((initial_time.minute % 2)==0) and (initial_time.second <=30)): #Print out peak RFI once every two minutes. 
			ind = np.argpartition(cal_dict['val'], -15)[-15:]
			print('15 largest RFI peaks (in dBm) and their respective frequencies (MHz): (' + str(cal_dict['val'][ind]) +'); ('+str(freqs[ind]/1.e6)+')')
			ind_FM_orbcomm= np.where(((freqs/1.e6) >= 88.) & ((freqs/1.e6) <=138.))
			ind_RFI_FM_orbcomm = np.argpartition(cal_dict['val'][ind_FM_orbcomm], -10)[-10:]
			print('10 largest peaks in FM through orbcomm (88 - 138): '+ str(np.array(cal_dict['val'][ind_FM_orbcomm])[ind_RFI_FM_orbcomm]) +'); ('+str(np.array(freqs[ind_FM_orbcomm])[ind_RFI_FM_orbcomm]/1.e6)+')')
		if OUTPUT_DATA:
			#print('Outputting data to file...')
			output_2D_data(start_time,(freqs/1.e6),cal_dict)
			print('Done.')
		if PLOT_EACH_SWEEP:
			plot_spectra(freqs,cal_dict['val'])
			
		if ( dt_init.seconds >= DT_EXIT): 
			#print('DELTA_MINS_UNTIL_END have passed. break.')
			bool_continue_schedule = False #sch
			break
		#print('Waiting until next sweep.')
		if (dt_start.seconds < 1):  #sch
			bool_wait_min_time = True #sch
			#print('Two sweeps were going to execute within 1s of eachother, wait at least 1s for next unique iteration of datetime.now() % (MOD_SEC_SCHEDULE = '+str(MOD_SEC_SCHEDULE)+') = 0.') #sch
		
		#open_and_plot() #SCH
	elif ((start_time.second % MOD_SEC_SCHEDULE != 0) and (bool_wait_min_time == True)): #SCH
		bool_wait_min_time = False
		#print('... waited until next second. Boolean "bool_wait_min_time" is now reset to default value of "False".')
		if ( (datetime.datetime.now() - initial_time).seconds >= DT_EXIT ):
			#print('Max scheduled time has also passed. End sweep scheduling.')
			break

#print('End of scheduled sweep.') #SCH
print('Initial time was '+str(initial_time)+'. It is now '+str(datetime.datetime.now())+'. Finis.')






