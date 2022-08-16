import inspect
import numpy as np
from scipy import interpolate

import spectral_distortions as sd
import foregrounds as fg
ndp = np.float64


class FisherEstimation:
    def __init__(self, fmin=7.5e9, fmax=3.e12, fstep=15.e9, \
                 duration=86.4, bandpass=True, fsky=0.7, mult=1., \
                 priors={'alps':0.1, 'As':0.1}, drop=0, doCO=False, doSOFTS=False, SOFTS_sensitivity_file=''):
        self.fmin = fmin
        self.fmax = fmax
        self.bandpass_step = 1.e8
        self.fstep = fstep
        self.duration = duration
        self.bandpass = bandpass
        self.fsky = fsky
        self.mult = mult
        self.priors = priors
        self.drop = drop
        self.doSOFTS = doSOFTS
        self.SOFTS_sensitivity_file = SOFTS_sensitivity_file

        self.setup()
        self.set_signals()

        if doCO:
            self.mask = ~np.isclose(115.27e9, self.center_frequencies, atol=self.fstep/2.)
        else:
            self.mask = np.ones(len(self.center_frequencies), bool)
        return

    def setup(self):
        if self.doSOFTS:
                self.set_SOFTS_frequencies()
                self.noise = self.SOFTS_sensitivity()
        else:
            self.set_frequencies()
            self.noise = self.pixie_sensitivity()
        return

    def run_fisher_calculation(self):
        N = len(self.args)
        F = self.calculate_fisher_matrix()
        for k in self.priors.keys():
            if k in self.args and self.priors[k] > 0:
                kindex = np.where(self.args == k)[0][0]
                F[kindex, kindex] += 1. / (self.priors[k] * self.argvals[k])**2
        normF = np.zeros([N, N], dtype=ndp)
        for k in range(N):
            normF[k, k] = 1. / F[k, k]
        self.cov = ((np.mat(normF, dtype=ndp) * np.mat(F, dtype=ndp)).I * np.mat(normF, dtype=ndp)).astype(ndp)
        #self.cov = np.mat(F, dtype=ndp).I
        self.F = F
        self.get_errors()
        return

    def get_errors(self):
        self.errors = {}
        for k, arg in enumerate(self.args):
            self.errors[arg] = np.sqrt(self.cov[k,k])
        return

    def print_errors(self, args=None):
        if not args:
            args = self.args
        for arg in args:
            #print arg, self.errors[arg], self.argvals[arg]/self.errors[arg]
            print(arg, self.argvals[arg]/self.errors[arg])

    def set_signals(self, fncs=None):
        if fncs is None:
            fncs = [sd.DeltaI_mu, sd.DeltaI_reltSZ_2param_yweight, sd.DeltaI_DeltaT,
                    fg.thermal_dust_rad, fg.cib_rad, fg.jens_freefree_rad,
                    fg.jens_synch_rad, fg.spinning_dust, fg.co_rad]
        self.signals = fncs
        self.args, self.p0, self.argvals = self.get_function_args()
        return

    def set_frequencies(self):
        if self.bandpass:
            self.band_frequencies, self.center_frequencies, self.binstep = self.band_averaging_frequencies()
        else:
            self.center_frequencies = np.arange(self.fmin + self.fstep/2., \
                                                self.fmax + self.fstep, self.fstep, dtype=ndp)[self.drop:]
        return

    def band_averaging_frequencies(self):
        #freqs = np.arange(self.fmin + self.bandpass_step/2., self.fmax + self.fstep, self.bandpass_step, dtype=ndp)
        freqs = np.arange(self.fmin + self.bandpass_step/2., self.fmax + self.bandpass_step + self.fmin, self.bandpass_step, dtype=ndp)
        binstep = int(self.fstep / self.bandpass_step)
        freqs = freqs[self.drop * binstep : (len(freqs) / binstep) * binstep]
        centerfreqs = freqs.reshape((len(freqs) / binstep, binstep)).mean(axis=1)
        #self.windowfnc = np.sinc((np.arange(binstep)-(binstep/2-1))/float(binstep))
        return freqs, centerfreqs, binstep

    def set_SOFTS_frequencies(self):
        if self.bandpass:
            print('no bandpass action in SOFTS yet. Exiting program...bye!')
            sys.exit()
        else:
            sdata = np.loadtxt('templates/'+self.SOFTS_sensitivity_file, dtype=ndp)
            fs = sdata[:, 0] * 1e9
            self.center_frequencies = fs
        return



    def pixie_sensitivity(self):
        sdata = np.loadtxt('templates/Sensitivities.dat', dtype=ndp)
        fs = sdata[:, 0] * 1e9
        sens = sdata[:, 1]
        template = interpolate.interp1d(np.log10(fs), np.log10(sens), bounds_error=False, fill_value="extrapolate")
        skysr = 4. * np.pi * (180. / np.pi) ** 2 * self.fsky
        if self.bandpass:
            N = len(self.band_frequencies)
            noise = 10. ** template(np.log10(self.band_frequencies)) / np.sqrt(skysr) * np.sqrt(15. / self.duration) * self.mult * 1.e26
            return (noise.reshape(( N / self.binstep, self.binstep)).mean(axis=1)).astype(ndp)
        else:
            return (10. ** template(np.log10(self.center_frequencies)) / np.sqrt(skysr) * np.sqrt(15. / self.duration) * self.mult * 1.e26).astype(ndp)

    def SOFTS_sensitivity(self):
        sdata = np.loadtxt('templates/'+self.SOFTS_sensitivity_file, dtype=ndp)
        fs = sdata[:, 0] * 1e9
        sens = sdata[:, 1]
        template = interpolate.interp1d(np.log10(fs), np.log10(sens), bounds_error=False, fill_value="extrapolate")
        skysr = 4. * np.pi * (180. / np.pi) ** 2 * self.fsky
        if self.bandpass:
            N = len(self.band_frequencies)
            noise = 10. ** template(np.log10(self.band_frequencies)) / np.sqrt(skysr) * np.sqrt(15. / self.duration) * self.mult * 1.e26
            return (noise.reshape(( N / self.binstep, self.binstep)).mean(axis=1)).astype(ndp)
        else:
            return (10. ** template(np.log10(self.center_frequencies)) / np.sqrt(skysr) * np.sqrt(15. / self.duration) * self.mult * 1.e26).astype(ndp)

    def get_function_args(self):
        targs = []
        tp0 = []
        for fnc in self.signals:
            argsp = inspect.getargspec(fnc)
            args = argsp[0][1:]
            p0 = argsp[-1]
            targs = np.concatenate([targs, args])
            tp0 = np.concatenate([tp0, p0])
        return targs, tp0, dict(zip(targs, tp0))

    def calculate_fisher_matrix(self):
        N = len(self.p0)
        F = np.zeros([N, N], dtype=ndp)
        for i in range(N):
            dfdpi = self.signal_derivative(self.args[i], self.p0[i])
            dfdpi /= self.noise
            for j in range(N):
                dfdpj = self.signal_derivative(self.args[j], self.p0[j])
                dfdpj /= self.noise
                #F[i, j] = np.dot(dfdpi, dfdpj)
                F[i, j] = np.dot(dfdpi[self.mask], dfdpj[self.mask])
        return F

    def signal_derivative(self, x, x0):
        h = 1.e-5 #default is 1E-4
        zp = 1. + h
        deriv = (self.measure_signal(**{x: x0 * zp}) - self.measure_signal(**{x: x0})) / (h * x0)
        return deriv

    def measure_signal(self, **kwarg):
        if self.bandpass:
            frequencies = self.band_frequencies
        else:
            frequencies = self.center_frequencies
        N = len(frequencies)
        model = np.zeros(N, dtype=ndp)
        for fnc in self.signals:
            argsp = inspect.getargspec(fnc)
            args = argsp[0][1:]
            if len(kwarg) and kwarg.keys()[0] in args:
                model += fnc(frequencies, **kwarg)
        if self.bandpass:
            #rmodel = model.reshape((N / self.binstep, self.binstep))
            #total = rmodel * self.windowfnc
            return model.reshape((N / self.binstep, self.binstep)).mean(axis=1)
            #return total.mean(axis=1)
        else:
            return model
