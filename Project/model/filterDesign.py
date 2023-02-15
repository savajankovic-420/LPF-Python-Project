#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy import signal
import scipy
import numpy as np
import sys
import math

# use generateSin/plotTime from the fourierTransform module
from fourierTransform import generateSin, plotTime


def lowpassfilter(fc,fs,ntaps):
	
	Normcutoff=(Fc/(Fs*2))
	h = np.zeros(ntaps)

	for i in range(ntaps):
		if i == ((ntaps-1)/2):
			h[i] = Normcutoff
		else:
			innerbracket = math.pi*Normcutoff*(i-((ntaps-1)/2))
			h[i] = Normcutoff*(np.sin(innerbracket)/innerbracket)

		if (True):
			h[i] = h[i] * (np.sin((i*math.pi)/ntaps))* (np.sin((i*math.pi)/ntaps))
			
	return h

def freqzPlot(coeff, Fs, msg):

	# find the frequency response using freqz from SciPy:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqz.html
	w, h = signal.freqz(coeff)

	# Reminder: np.pi rad/sample is actually the Nyquist frequency
	w = w * Fs/(2*np.pi) # needed to draw the frequency on the X axis

	# plots the magnitude response where the x axis is normalized in rad/sample
	fig, ax1 = plt.subplots()
	ax1.set_title('Digital filter frequency response (' + msg + ')')
	ax1.plot(w, 20 * np.log10(abs(h)), 'b')
	ax1.set_ylabel('Amplitude [dB]', color='b')
	ax1.set_xlabel('Frequency [Hz]')

	# uncomment the lines below if you wish to inspect the phase response
	# Note: as important as the phase response is for some applications,
	# it is not critical at this stage because we expect a linear phase in the passband

	# ax2 = ax1.twinx()
	# angles = np.unwrap(np.angle(h))
	# ax2.plot(w, angles, 'g')
	# ax2.set_ylabel('Angle (radians)', color='g')



def filterSin(Fs, Fc, coeff):

	# we can control the frequency relative to the filter cutoff
	time, x = generateSin(Fs, interval = 1.0, frequency = Fc* 15)
	plotTime(x, time)

	# use lfilter from SciPy for FIR filtering:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html
	fx = signal.lfilter(coeff, 1.0, x)

	# you should clearly observe the effects (attenuation, delay) introduced by the filter
	plotTime(fx, time)

def cli_error_msg():

	# error message to provide the correct command line interface (CLI) arguments
	print('Valid arguments:')
	print('\trc:  reference code')
	print('\til1: in-lab 1')
	print('\til2: in-lab 2')
	print('\tth:  take-home')
	sys.exit()

if __name__ == "__main__":

	if len(sys.argv[0:]) != 2:
		cli_error_msg()

	Fs = 100.0           # sampling rate
	Fc = 5.0           # cutoff frequency
	N_taps = 51       # number of taps for the FIR


	if (sys.argv[1] == 'rc'): # runs the reference code (rc)

		print('Reference code for the digital filter design')

		# derive filter coefficients using firwin from Scipy:
		# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin.html
		# second argument is the normalized cutoff frequency, i.e., the
		# cutoff frequency divided by Nyquist frequency (half of sampling rate)
		firwin_coeff = signal.firwin(N_taps, Fc/(Fs/2), window=('hann'))

		# plot the frequency response obtained through freqz
		freqzPlot(firwin_coeff, Fs, 'firwin for ' + str(int(Fc)) + ' Hz cutoff with ' + str(N_taps) + ' taps')
		filterSin(Fs, Fc, firwin_coeff)
	elif (sys.argv[1] == 'il1'):

		print('In-lab experiment 1 for the digital filter design')

		# implement your own method for finding the coefficients for a low pass filter
		# my_own_coeff = ... provide the following arguments: Fc, Fs and N_taps
		# compare through visual inspection the frequency response against firwin
		# freqzPlot(my_own_coeff, Fs, 'my own FIR design with ' + str(N_taps) + ' taps')
		my_own_coeff = lowpassfilter(Fc,Fs,N_taps)
		freqzPlot(my_own_coeff,Fs, 'my own FIR design with' + str(N_taps) + 'taps' )

		# comparison between the two.
		firwin_coeff = scipy.signal.firwin(N_taps, Fc/(Fs/2), window=('hann'))
		freqzPlot(firwin_coeff, Fs, 'firwin for ' + str(int(Fc)) + ' Hz cutoff with ' + str(N_taps) + ' taps')

	elif (sys.argv[1] == 'il2'):

		print('In-lab experiment 2 for the digital filter design')
		interval = 1
		#inserting the three signals from fourierTransform.py.
		time, tone1 = generateSin(Fs, interval, 20,3,0.0)      ## amplitude, frequency, step.
		time, tone2 = generateSin(Fs, interval, 30,4,0.5)
		time, tone3 = generateSin(Fs, interval, 40,4,1.0)

		sigcomb = tone1 + tone2 + tone3

		plotTime(tone1,time)  # plotting the time domain of the three signals,
		plotTime(tone2,time)
		plotTime(tone3,time)
		plotTime(sigcomb,time)        # plotting the combination of the three signals from the fouriertransform.

		scipy.signal.lfilter(lowpassfilter(Fc,Fs,N_taps),interval, sigcomb)

		filtertest = scipy.signal.lfilter(lowpassfilter(Fc,Fs,N_taps),interval, sigcomb)

		plotTime(filtertest, time)


	elif (sys.argv[1] == 'th'):
		interval = 1
		print('Take-home exercise for the digital filter design')
		lowercut = 25
		highercut = 35
		Fc = [lowercut, highercut]

		nyquisteq = (Fs/2) # definition of a nyquist rate

		cutoffeq = [Fc[0]/nyquisteq,Fc[1]/nyquisteq]

		time, tone1 = generateSin(Fs, interval, 20,3,0.0)      #frequencies being 20,30,40. The middle one is bound to stay due to the cutoff equations that filter below 25 and above 35.
		time, tone2 = generateSin(Fs, interval, 30,4,0.5)
		time, tone3 = generateSin(Fs, interval, 40,4,1.0)

		sigcomb = tone1 + tone2 + tone3
		plotTime(sigcomb,time)
		scipy.signal.firwin(N_taps,cutoffeq, window='hann', pass_zero='bandpass')
		firwin_coeff = scipy.signal.firwin(N_taps,cutoffeq, window='hann', pass_zero='bandpass')
		filtering = scipy.signal.lfilter(firwin_coeff, 1, sigcomb)



		freqzPlot(firwin_coeff, Fs,'firwin filtering with lower frequency cut of '+ str(int(lowercut)) +" and a higher frequency cut of "+ str(int(highercut)))
		plotTime(filtering,time)

	else:

		cli_error_msg()

	plt.show()
