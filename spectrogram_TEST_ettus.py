from fft_accum_avg_p1_rts_2048_ettus import fft_accum_avg_p1_rts_2048_ettus
import matplotlib.pyplot as plt
import numpy as np
import scipy
import time
import peakutils

#GLOBAL PARAMS # FOR LORDS BRIDGE TEST: SAMP,START,STOP = 25,30,490 WIH BW OVERLAP OF 5.
REMOVE_CLOCK_AND_SIDE_CHANNELS = False
SAMP_RATE = 25.0e6 
SIZE_FFT = 2048
N_TO_AVG = int((SAMP_RATE/SIZE_FFT)*0.25) #0.25 is canonical.
F_START = 30.0e6#35.0e6 
F_STOP = 490.0e6
BW_OVERLAP = 5.0e6

STR_FILE_TEMP = '/Users/josaitis/gnuradio_tutorials/fft_accum_temp/fft_accum_temp_TEST.txt'
STR_FILE_MIN = '/Users/josaitis/gnuradio_tutorials/fft_accum_temp/fft_accum_min.txt'
STR_FILE_MAX = '/Users/josaitis/gnuradio_tutorials/fft_accum_temp/fft_accum_max.txt'
STR_FILE_MEAN = '/Users/josaitis/gnuradio_tutorials/fft_accum_temp/fft_accum_mean.txt'

SIDE_CHAN_TO_DEL=0
CENTER_CHAN_TO_DEL=5
SPURIOUS_CHAN_TO_DEL = 5 # must be odd number

ii_first_keep = int( ((BW_OVERLAP/2)/SAMP_RATE)*SIZE_FFT )
ii_last_keep = int( (1-(BW_OVERLAP/2)/SAMP_RATE)*SIZE_FFT )
print('ii first and last to keep (respectively): '+str(ii_first_keep)+', '+str(ii_last_keep))


###### Functions
def remove_clock_effects(arr, freqs):
	spurious_mask = np.arange(10.0*1.e6, (F_STOP+(SAMP_RATE - BW_OVERLAP)/2), 5.0*1.e6)
	spurious_indices = np.searchsorted(freqs,spurious_mask)
	for ii in spurious_indices: #index of the integer multiples of 5MHz
		if (ii == 0): spurious_indices = np.append(spurious_indices,np.arange(ii, ii + SPURIOUS_CHAN_TO_DEL/2 +1,1))
		spurious_indices = np.append(spurious_indices, np.arange(ii - SPURIOUS_CHAN_TO_DEL/2, ii + SPURIOUS_CHAN_TO_DEL + 1,1))
	spurious_indices=np.unique(spurious_indices)
	arr[spurious_indices] = np.nan
	return arr


def analysis_of_sample():
	file = scipy.fromfile(open(STR_FILE_TEMP), dtype=scipy.float32)
	#print( 'size of temp file: ' + str(np.size(file)))
	dat1D = np.array(file)
	#print('  In analysis_of_sample: size of file: '+str(np.size(file))+', SIZE_FFT: '+str(SIZE_FFT)+', N_TO_AVG (round) is : '+str(N_TO_AVG))
	#print('size(file)/SIZE_FFT: '+str(np.size(file)/SIZE_FFT))
	dat2D = np.reshape( file,(np.size(file)/SIZE_FFT,SIZE_FFT) ) #order them by every size_fft cols to get rid of edge and clock peaks
	if REMOVE_CLOCK_AND_SIDE_CHANNELS:
		elims = np.arange(int(SIZE_FFT/2 - CENTER_CHAN_TO_DEL), int(SIZE_FFT/2 +CENTER_CHAN_TO_DEL +1),1) #elimiate the clock peaks and the N channels on each side
		elims = np.concatenate((elims, np.arange(0,int(SIDE_CHAN_TO_DEL)),np.arange(int(SIZE_FFT - SIDE_CHAN_TO_DEL),int(SIZE_FFT))))
		#elims = np.concatenate((  elims, np.arange( (((F_START/1e6)%5)/(SAMP_RATE/1.e6))*SIZE_FFT,SIZE_FFT,(5/(SAMP_RATE/1.e6))*SIZE_FFT ).astype(int)  )) # spurious signals can occur in integer multiples of 5 due to HackRF clock @ 10 MHz and 25MHz crystal on the SIxxxxxx chip.
		elims = np.unique(elims)

		for ii in elims: #set the power at each clock peak and edge (+/- side_chan_to_del channels) to zero
			dat2D[:,ii]=0


	#plt.semilogy(dat2D[100,:],'g',label = 'spectra # 100')
	#plt.semilogy(dat2D[700,:],'m',label = 'spectra # 700')
	#plt.semilogy(dat2D[1100,:],'r',label = 'spectra # 1100')
	#plt.semilogy(dat2D[5000,:],'k',label = 'spectra # 5000')
	##plt.legend()
	#plt.show()
	#link for following line: https://stackoverflow.com/questions/30379311/fast-way-to-take-average-of-every-n-rows-in-a-npy-array
	#note, the power(,2) below is to convert from measured amplitude to measured power. Confirmed in by doing input power tests.
	n_to_avg_effective = int(np.size(dat1D)/2048)
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
	open(STR_FILE_TEMP, 'w').close()
	open(STR_FILE_MEAN, 'w').close()
	open(STR_FILE_MIN, 'w').close()
	open(STR_FILE_MAX, 'w').close()

def open_and_plot():
	arr_min = np.power(np.loadtxt(STR_FILE_MIN),2)
	arr_max = np.power(np.loadtxt(STR_FILE_MAX),2)
	arr_mean = np.power(np.loadtxt(STR_FILE_MEAN),2)

	x=np.linspace((F_START-(SAMP_RATE- BW_OVERLAP)/2),(F_STOP+(SAMP_RATE - BW_OVERLAP)/2),int(ii_last_keep-ii_first_keep)*n_FFTs)


	#arr_mean = remove_clock_effects(arr_mean,x)


	plt.semilogy(x/1.e6,arr_mean,'k', label = 'Mean')
	#plt.ylim((1.e-2,1.e4))
	#plt.semilogy(x,np.ravel(arr_min),'b',label='Min')
	#plt.semilogy(x,np.ravel(arr_max),'r',label='Max')
	i_peak = peakutils.indexes(arr_mean, thres=0.3)
	print('Peak: ('+str(x[i_peak]/1.e6)+' ,'+ str(arr_mean[i_peak])+')')
	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Power (Arbitrary HackRF Units)')
	plt.title('Online FFT, Mean Power')
	plt.legend()
	#plt.savefig(str('/Users/josaitis/gnuradio_tutorials/fft_accum_temp/'+str(F_START/1.e6)+'-'+str(F_STOP/1.e6)+'MHz_samp_rate_'+str(SAMP_RATE/1.e6)+'MHz.pdf'))
	plt.show()
# Main

n_FFTs = int(((F_STOP-F_START)/(SAMP_RATE- BW_OVERLAP) +1))
freqs = np.arange(F_START,F_STOP + (SAMP_RATE- BW_OVERLAP), SAMP_RATE- BW_OVERLAP)

print('flushing temporary and output files...')
file_flush()
print('done.')

FFT_to_file = fft_accum_avg_p1_rts_2048_ettus()
FFT_to_file.set_str_file_out(STR_FILE_TEMP)
FFT_to_file.set_samp_rate(SAMP_RATE)
print('Number of final, averaged, spectra that will be made: '+str(n_FFTs))
print('Number of FFTs to average in order to make one of the above spectra '+str(N_TO_AVG))

for ff in freqs:
	#FFT_to_file = fft_accum_avg_p1_rts()
	#FFT_to_file.set_str_file_out(STR_FILE_TEMP)
	#FFT_to_file.set_samp_rate(SAMP_RATE)
	FFT_to_file.set_n_head(N_TO_AVG)
	print('setting f_c to '+str(ff/1.e6)+' MHz. N_TO_AVG is '+str(N_TO_AVG))
	FFT_to_file.set_f_c(ff)
	print('f_c set to: '+str( (FFT_to_file.get_f_c()/1.e6) ) +' MHz')
	FFT_to_file.start()
	FFT_to_file.wait()
	FFT_to_file.stop()
	print('Just before reset, head is: '+str(FFT_to_file.get_n_head()))
	FFT_to_file.blocks_head_0.reset() #Hopefully this resets the head block counter.
	#del(FFT_to_file)
	analysis_of_sample()
	

FFT_to_file.stop()
del(FFT_to_file)

open_and_plot()




