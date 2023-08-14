import numpy as np
import pandas as pd
import uncalled as unc
# import seaborn as sns
# from matplotlib import pyplot as plt
import pickle

from .utils import *

class Bin:
    def __init__(self, nbins=64, poremodel=model_6mer, 
                fixed=True, bounds=None, clip=True) -> None:
        if type(poremodel) == unc.pore_model.PoreModel:
            self.poremodel = poremodel
        else:
            self.poremodel = unc.PoreModel(df = poremodel)
        self.bounds = bounds
        self.nbins = nbins
        self.clip = clip
        if clip:
            self.clipbounds = (self.poremodel['current.mean'].min(), self.poremodel['current.mean'].max())
        if fixed:
            self._get_fixed_bins()
        else:
            self._get_density_bins()
        self.nbins = len(self.starts) + 1
        self.kmer_to_bin = np.searchsorted(self.starts, self.poremodel.means)
        self.binmodel = [(np.mean(self.poremodel.means[np.where(self.kmer_to_bin == idx)[0]]), 
                                    np.mean(self.poremodel.stdvs[np.where(self.kmer_to_bin == idx)[0]]))
                                    for idx in range(self.nbins)]
    def _get_fixed_bins(self):
        model = self.poremodel.to_df()
        # df has columns kmer, mean, stdv
        if self.bounds:
            minc, maxc = self.bounds
            curr_range = (model['current.mean'].max() - maxc) + (minc - model['current.mean'].min())
            step = curr_range / self.nbins
            lower = np.arange(model['current.mean'].min(), minc, step)
            upper = np.arange(maxc, model['current.mean'].max(), step)
            starts = np.concatenate([lower, upper])
        else:
            minc, maxc = model['current.mean'].min(), model['current.mean'].max()  
            starts = np.linspace(minc, maxc, self.nbins, endpoint=False)
        self.starts = starts
    
    def _get_density_bins(self):
        k = self.poremodel.K
        model = self.poremodel.to_df().sort_values(by='mean')
        if self.bounds:
            model = model[(model['current.mean'] <= self.bounds[0]) | (model['current.mean'] >= self.bounds[1])]
        kmers_sorted = list(model['kmer'])
        num_kmers = len(kmers_sorted)
        binsize = num_kmers // self.nbins
        self.starts = np.array([float(self.poremodel[[kmers_sorted[idx]]]) for idx in range(0, num_kmers, binsize)])

    # def viz_bins(self):
    #     fig, ax = plt.subplots()
    #     sns.histplot(self.poremodel.to_df()['mean'], ax=ax, bins=50)
    #     for s in self.starts:
    #         ax.axvline(s)
    #     return ax

    def _preprocess(self, signal, evdt=None, normalize=True):
        if evdt:
            signal = np.array(evdt.get_means(signal))
        if normalize:   
            signal, _, _ = normalize_signal(signal, self.poremodel)
        return signal
    def bin_signal(self, signal, evdt=None, normalize=True):
        signal = np.array(signal)
        signal = self._preprocess(signal, evdt=evdt, normalize=normalize)
        if self.clip:
            signal = signal[(signal >= self.clipbounds[0]) & (signal <= self.clipbounds[1])]
        if self.bounds:
            signal = signal[(signal <= self.bounds[0]) | (signal >= self.bounds[1])]
        
        return np.searchsorted(self.starts, signal)

    def bin_sequence(self, seq):
        sig = seq_to_sig(self.poremodel, seq)
        if self.bounds:
            sig = sig[(sig <= self.bounds[0]) | (sig >= self.bounds[1])]
        
        return np.searchsorted(self.starts, sig)

class HPCBin(Bin):
    def __init__(self, nbins=64, poremodel=model_6mer, bounds=None, clip=False) -> None:
        if type(poremodel) == unc.pore_model.PoreModel:
            self.poremodel = poremodel
        else:
            self.poremodel = unc.PoreModel(df = poremodel)
        self.bounds = bounds
        self.nbins = nbins
        self.clip = clip
        if clip:
            self.clipbounds = (self.poremodel['current.mean'].min(), self.poremodel['current.mean'].max())
        self.minc, self.maxc = self.poremodel['current.mean'].min(), self.poremodel['current.mean'].max()  
        self.space = (self.maxc - self.minc) / self.nbins
        self.kmer_to_bin = self.signal_to_binseq(self.poremodel['current.mean'])
        self.binmodel = [(np.mean(self.poremodel['current.mean'][np.where(self.kmer_to_bin == idx)[0]]), 
                                    np.mean(self.poremodel['current.stdv'][np.where(self.kmer_to_bin == idx)[0]]))
                                    for idx in range(self.nbins)]
    def _hpc(self, binseq):
        return binseq[np.insert(np.diff(binseq) != 0, 0, True)]
    def signal_to_binseq(self, sig):
        return np.clip(np.floor((sig - self.minc - 1e-10) / self.space), 0, self.nbins - 1).astype(np.uint8)
    def bin_signal(self, signal, evdt=None, normalize=True):
        if not isinstance(signal, np.ndarray):
            signal = np.array(signal)
        signal = self._preprocess(signal, evdt=evdt, normalize=normalize)
        if self.clip:
            signal = signal[(signal >= self.clipbounds[0]) & (signal <= self.clipbounds[1])]
        if self.bounds:
            signal = signal[(signal <= self.bounds[0]) | (signal >= self.bounds[1])]
        if len(signal) == 0:
            return []
        return self._hpc(self.signal_to_binseq(signal))
    def bin_sequence(self, seq, revcomp=False, shred_size=0):
        kmers = seq_to_kmer(self.poremodel, seq, revcomp=revcomp)
        if shred_size > 0:
            kmers = np.concatenate([s for _, s in kmers])
            for idx in range(0, len(kmers), shred_size):
                shred = kmers[idx:idx+shred_size]
                if idx + shred_size >= len(kmers):
                    yield np.append(self._hpc(self.kmer_to_bin[shred]), -1)
                else:
                    yield self._hpc(self.kmer_to_bin[shred])
        else:
            print('binning %s'%seq)
            for sid, seq in kmers:
                yield (sid, self._hpc(self.kmer_to_bin[seq]))
    # def viz_bins(self):
    #     fig, ax = plt.subplots()
    #     sns.histplot(self.poremodel.to_df()['mean'], ax=ax, bins=50)
    #     for idx in range(1, self.nbins):
    #         ax.axvline(self.minc + (idx * self.space))
    #     return ax
    def save_bins(self, fname):
        # save the defining variables
        pickle.dump((self.nbins, self.poremodel.to_df()), open(fname, 'wb'))
    @classmethod
    def from_pickle(fname):
        nbins, poremodel  = pickle.load(open(fname, 'rb'))
        return HPCBin(nbins=nbins, poremodel=unc.PoreModel(df = poremodel))
    
class SigProcHPCBin(HPCBin):
    def __init__(self, nbins=64, poremodel=model_6mer, bounds=None, clip=False) -> None:
        super(SigProcHPCBin, self).__init__(nbins=nbins, poremodel=poremodel, bounds=bounds, clip=clip)
    def _preprocess(self, signal, evdt=None, normalize=True):
        if evdt:
            signal = np.array(evdt.get_means(signal))
        mask = pd.Series(signal).rolling(25, center=True).std().rolling(25, center=True, min_periods=1).max() > 5
        signal = signal[mask.to_numpy()]
        if normalize:   
            signal, _, _ = normalize_signal(signal, self.poremodel)
        return signal