#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
import numpy as np
import cmath, math
import sys
from scipy import signal

def plotSpectrum(x, Fs, type = 'FFT'):

	n = len(x)             # length of the signal
	df = Fs/n              # frequency increment (width of freq bin)

	# compute Fourier transform, its magnitude and normalize it before plotting
	if type == 'FFT':
		Xfreq = np.fft.fft(x)
	elif type == 'DFT':
		Xfreq = my_DFT(x)
	elif type == 'IDFT':
		Xfreq = my_IDFT(x)
	XMag = abs(Xfreq)/n

	# Note: because x is real, we keep only the positive half of the spectrum
	# Note also: half of the energy is in the negative half (not plotted)
	XMag = XMag[0:int(n/2)]

	# freq vector up to Nyquist freq (half of the sample rate)
	freq = np.arange(0, Fs/2, df)

	fig, ax = plt.subplots()
	ax.plot(freq, XMag)
	ax.set(xlabel='Frequency (Hz)', ylabel='Magnitude',
		title='Frequency domain plot')
	# fig.savefig("freq.png")
	plt.show()

def plotTime(x, time):

	fig, ax = plt.subplots()
	ax.plot(time, x)
	ax.set(xlabel='Time (sec)', ylabel='Amplitude',
			title='Time domain plot')
	# fig.savefig("time.png")
	plt.show()

def plotTime(x, time, type='IDFT'):
    x = my_IDFT(my_DFT(x))
    fig, ax = plt.subplots()
    ax.plot(time, x)
    ax.set(xlabel='Time (sec)', ylabel='Amplitude',
            title='Time domain plot')
    # fig.savefig("time.png")
    plt.show()

def generateSin(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0):

	dt = 1.0/Fs                          # sampling period (increment in time)
	time = np.arange(0, interval, dt)    # time vector over interval

	# generate the sin signal
	x = amplitude*np.sin(2*math.pi*frequency*time+phase)

	return time, x

def cli_error_msg():

	# error message to provide the correct command line interface (CLI) arguments
	print('Valid arguments:')
	print('\trc:  reference code')
	print('\til1: in-lab 1')
	print('\til2: in-lab 2')
	print('\til3: in-lab 3')
	print('\tth:  take-home')
	sys.exit()

def my_DFT(x):

	N=len(x)
	m = np.arange(N)      # setting the range of the array from 0 to 99 (100 points)
	k = m.reshape((N,1))  # changing the shape of the array.
	exp = np.exp(-2j*cmath.pi*((k*m)/N))

	dftx = np.dot(x, exp)

	return dftx


def my_IDFT(x):

	invN= len(x)
	invm = np.arange(invN)
	invk = invm.reshape((invN,1))
	invexp = np.exp(2j*cmath.pi*((invk*invm)/invN))

	idtfx = np.dot(x, invexp)

	return idtfx

if __name__ == "__main__":
	if len(sys.argv[0:]) != 2:
		cli_error_msg()

	Fs = 100.0          # sampling rate
	interval = 1.0      # set up to one full second

	if (sys.argv[1] == 'rc'): # runs the reference code (rc)

		print('Reference code for the Fourier transform')

		# generate the user-defined sin function
		time, x = generateSin(Fs, interval) ## can be changed within the function.
		# plot the signal in time domain
		plotTime(x, time)
		# plot the signal in frequency domain
		plotSpectrum(x, Fs, type = 'FFT')

	elif (sys.argv[1] == 'il1'):

		print('In-lab experiment 1 for the Fourier transform')

		# compute the spectrum with your own DFT
		# you can use cmath.exp() for complex exponentials
		# plotSpectrum(x, Fs, type = 'your DFT name')
		time, x = generateSin(Fs, interval)
		# confirm DFT/IDFT correctness by checking if x == IDFT(DFT(x))
		plotSpectrum(x,Fs,type = 'DFT')
		plotTime(x,time,type = 'IDFT')




		dftverification = np.allclose(x,my_IDFT(my_DFT(x)))
		idftverification = np.allclose(np.fft.fft(x), my_DFT(x))

		print(dftverification)
		print(idftverification)

	elif (sys.argv[1] == 'il2'):

		print('In-lab experiment 2 for the Fourier transform')

		# use np.random.randn() for randomization
		# we can overwrite the default values
		# frequency =  8.0                     # frequency of the signal
		# amplitude =  3.0                     # amplitude of the signal
		# phase = 1.0                          # phase of the signal
		# time, x = generateSin(Fs, interval, frequency, amplitude, phase)

		frequency = np.random.randn()
		amplitude = np.random.randn()
		phase = np.random.randn()
		time, x = generateSin(1000,1,frequency,amplitude,phase)
		plotTime(x,time, type = 'DFT')
		plotSpectrum(x,time, type = 'IDFT')
		plotTime(x,time, type = 'IDFT')
		plotSpectrum(x, time, type = 'DFT')

		# You should also numerically check if the signal energy
		# in time and frequency domains is identical

		# for further details, if any, check the lab document

	elif (sys.argv[1] == 'il3'):

		print('In-lab experiment 3 for the Fourier transform')

		# generate randomized multi-tone signals

		time, tone1 = generateSin(Fs, interval, 2,3,0.0)      ## amplitude, frequency, step.
		time, tone2 = generateSin(Fs, interval, 3,4,0.5)
		time, tone3 = generateSin(Fs, interval, 3,4,1.0)

		x = tone1 + tone2 + tone3

		#plot the signal in time domains
		plotTime(tone1, time)
		plotTime(tone2, time)
		plotTime(tone3, time)

		plotSpectrum(x,Fs,type="FFT")
		# for further details, if any, check the lab document

	elif (sys.argv[1] == 'th'):
				import random
				print('Take-home exercise for the Fourier transform')
				dutycycle = random.uniform(0,1)  # initializing the random value of the duty cycle
				dutycycletxt = round(dutycycle,3)  # printing the value for the user.
				print(dutycycletxt)
				#would pick a random number between 0 and 1 (that is a float).
				# This would provide a solid duty cycle.
				interval= 1
				timeint = 25*interval # repetition of the interval so the range is extended to 25.

				dt = 1/Fs     # variable is divided by the sampling rate.
				time = np.arange(0, timeint, dt) # the axes arranged for the duty cycle.

				dutycyclesquare = signal.square(time, dutycycle) 
				plotTime(dutycyclesquare, time) # generating a time plot for the square duty cycle
				plotSpectrum(dutycyclesquare, Fs, type = 'DFT') #the frequency spectrum is plotted alongside the time plot.

	else:

		cli_error_msg()

	plt.show()
