import numpy as np
import pandas as pd
import uncalled as unc
import seaborn as sns
from matplotlib import pyplot as plt
from pomegranate import NormalDistribution, HiddenMarkovModel
import subprocess as proc
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

    def viz_bins(self):
        fig, ax = plt.subplots()
        sns.histplot(self.poremodel.to_df()['mean'], ax=ax, bins=50)
        for s in self.starts:
            ax.axvline(s)
        return ax

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

    
# class HPCBin_slow(Bin):
#     def __init__(self, nbins=64, poremodel=model_6mer, fixed=True, bounds=None, clip=False) -> None:
#         super(HPCBin, self).__init__(nbins=nbins, poremodel=poremodel, fixed=fixed, bounds=bounds, clip=clip)
#     def bin_signal(self, signal, evdt=None, normalize=True):
#         return self._hpc(super(HPCBin, self).bin_signal(signal, evdt=evdt, normalize=normalize))
#     def _hpc(self, binseq):
#         return binseq[np.insert((np.where(np.diff(binseq) != 0)[0] + 1), 0, 0)]
#     def bin_sequence(self, seq):
#         return self._hpc(super(HPCBin, self).bin_sequence(seq))
#     def write_bin_sequence(self, seq, fname = 'ref.bins', revcomp = True):
#         k = self.poremodel.K
#         f = open(fname, 'wb')
#         if not os.path.exists(seq + '.pac'):
#             proc.call(['uncalled', 'index', seq, '--no-bwt'])
#         idx = unc.load_index(self.poremodel.K, seq, load_bwt=False)
#         print('loaded seq index')
#         for sid, length in idx.get_seqs():
#             kmers = np.array(idx.get_kmers(sid, 0, length))
#             kmers = kmers[kmers < 4**k]
#             bins = self.kmer_to_bin[kmers]
#             bins = self._hpc(bins)
#             # pack bins and write to file
#             charseq = bins + 3 # offset so min bin is 3 (0,1,2 are protected chars)
#             assert charseq.max() < 256
#             byteseq = charseq.astype('=B').tobytes()
#             f.write(byteseq)
#             if revcomp:
#                 revcomp_kmers = np.apply_along_axis(self.poremodel.kmer_revcomp, 0, kmers[::-1])
#                 kmers = revcomp_kmers
#                 bins = self.kmer_to_bin[kmers]
#                 bins = self._hpc(bins)
#                 # pack bins and write to file
#                 charseq = bins + 3 # offset so min bin is 3 (0,1,2 are protected chars)
#                 assert charseq.max() < 256
#                 byteseq = charseq.astype('=B').tobytes()
#                 f.write(byteseq)
#         f.close()

class HPCBin(Bin):
    def __init__(self, nbins=64, poremodel=model_6mer, fixed=True, bounds=None, clip=False) -> None:
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
        return binseq[np.insert((np.where(np.diff(binseq) != 0)[0] + 1), 0, 0)]
    def signal_to_binseq(self, sig):
        return np.clip(np.floor((sig - self.minc - 1e-10) / self.space), 0, self.nbins - 1).astype(int)
    def bin_signal(self, signal, evdt=None, normalize=True):
        signal = np.array(signal)
        signal = self._preprocess(signal, evdt=evdt, normalize=normalize)
        if self.clip:
            signal = signal[(signal >= self.clipbounds[0]) & (signal <= self.clipbounds[1])]
        if self.bounds:
            signal = signal[(signal <= self.bounds[0]) | (signal >= self.bounds[1])]
        return self._hpc(self.signal_to_binseq(signal))
    def bin_sequence(self, seq, revcomp=False, shred_size=0):
        kmers = seq_to_kmer(self.poremodel, seq, revcomp=revcomp)
        if shred_size > 0:
            shreds = []
            for idx in range(0, len(kmers), shred_size):
                shred = kmers[idx:idx+shred_size]
                shreds.append(self._hpc(self.kmer_to_bin[shred]))
            return shreds
        return self._hpc(self.kmer_to_bin[kmers])
    def viz_bins(self):
        fig, ax = plt.subplots()
        sns.histplot(self.poremodel.to_df()['mean'], ax=ax, bins=50)
        for idx in range(1, self.nbins):
            ax.axvline(self.minc + (idx * self.space))
        return ax
    def save_bins(self, fname):
        # save the defining variables
        pickle.dump((self.nbins, self.poremodel.to_df()), open(fname, 'wb'))
    @classmethod
    def from_pickle(fname):
        nbins, poremodel  = pickle.load(open(fname, 'rb'))
        return HPCBin(nbins=nbins, poremodel=unc.PoreModel(df = poremodel))
    


class DeltaBin(Bin):
    def __init__(self, nbins=64, poremodel=model_6mer, fixed=True, bounds=None, clip=True) -> None:
        if type(poremodel) == unc.pore_model.PoreModel:
            self.mainmodel = poremodel
        else:
            self.mainmodel = unc.PoreModel(df = poremodel)
        self.poremodel = model_deltas(poremodel)
        self.nbins = nbins
        self.bounds = bounds
        self.clip = clip
        if clip:
            self.clipbounds = (self.poremodel['current.mean'].min(), self.poremodel['current.mean'].max())
        if fixed:
            self._get_fixed_bins()
        else:
            self._get_density_bins()
    def _preprocess(self, signal, evdt=None, normalize=True):
        if evdt:
            events = evdt.get_events(signal)
            signal = np.array([curr.mean for curr in events])
        if normalize:
            signal, _, _ = normalize_signal(signal, self.mainmodel)
        signal = sig_to_delta(signal)
        return signal


class StayHMMBin(Bin):
    def __init__(self, nbins=64, poremodel=model_6mer, fixed=True, bounds=None, clip=False, p_stay=0.5) -> None:
        super(StayHMMBin, self).__init__(nbins=nbins, poremodel=poremodel, fixed=fixed, bounds=bounds, clip=clip)
        self.p_stay = p_stay
        self.hmm = self._hmm_stay_model()
    def bin_signal(self, signal, evdt=None, normalize=True):
        signal = np.array(signal)
        signal = self._preprocess(signal, evdt=evdt, normalize=normalize)
        if self.clip:
            signal = signal[(signal >= self.clipbounds[0]) & (signal <= self.clipbounds[1])]
        if self.bounds:
            signal = signal[(signal <= self.bounds[0]) | (signal >= self.bounds[1])]
        # bins = np.searchsorted(self.starts, signal)
        bins = self._viterbi_signal(signal)
        return self._hpc(bins)
    def _hpc(self, binseq):
        return binseq[np.insert((np.where(np.diff(binseq) != 0)[0] + 1), 0, 0)]
    def _viterbi_signal(self, signal):
        bins = self.hmm.viterbi(signal)[1]
        bins = np.array([int(state.name) for _, state in bins[1:]])
        return bins
    def write_bin_sequence(self, seq, fname = 'ref.bins', revcomp = True):
        k = self.poremodel.K
        f = open(fname, 'wb')
        if not os.path.exists(seq + '.pac'):
            proc.call(['uncalled', 'index', seq, '--no-bwt'])
        idx = unc.load_index(self.poremodel.K, seq, load_bwt=False)
        print('loaded seq index')
        for sid, length in idx.get_seqs():
            kmers = np.array(idx.get_kmers(sid, 0, length))
            kmers = kmers[kmers < 4**k]
            sig = self.poremodel[kmers]
            if self.bounds:
                sig = sig[(sig <= self.bounds[0]) | (sig >= self.bounds[1])]
            bins = self.hmm.viterbi(sig)[1]
            bins = np.array([int(state.name) for _, state in bins[1:]])
            bins = self._hpc(bins)
            # pack bins and write to file
            charseq = bins + 3 # offset so min bin is 3 (0,1,2 are protected chars)
            assert charseq.max() < 256
            byteseq = charseq.astype('=B').tobytes()
            f.write(byteseq)
            if revcomp:
                revcomp_kmers = np.apply_along_axis(self.poremodel.kmer_revcomp, 0, kmers[::-1])
                kmers = revcomp_kmers
                sig = self.poremodel[kmers]
                if self.bounds:
                    sig = sig[(sig <= self.bounds[0]) | (sig >= self.bounds[1])]
                bins = self.hmm.viterbi(sig)[1]
                bins = np.array([int(state.name) for _, state in bins[1:]])
                bins = self._hpc(bins)
                # pack bins and write to file
                charseq = bins + 3 # offset so min bin is 3 (0,1,2 are protected chars)
                assert charseq.max() < 256
                byteseq = charseq.astype('=B').tobytes()
                f.write(byteseq)
        f.close()
    def bin_sequence(self, seq):
        sig = seq_to_sig(self.poremodel, seq)
        if self.bounds:
            sig = sig[(sig <= self.bounds[0]) | (sig >= self.bounds[1])]
        # bins = np.searchsorted(self.starts, sig)
        bins = self.hmm.viterbi(sig)[1]
        bins = np.array([int(state.name) for _, state in bins[1:]])
        return self._hpc(bins)
    def _hmm_stay_model(self):
        kmer_to_bin = np.searchsorted(self.starts, self.poremodel.means)
        bin_overlap = np.zeros((self.nbins, self.nbins))
        for mer in self.poremodel.KMER_STRS:
            oldbin = kmer_to_bin[self.poremodel.str_to_kmer(mer)]
            for new in 'ACGT':
                newmer = mer[1:] + new
                newbin = kmer_to_bin[self.poremodel.str_to_kmer(newmer)]
                bin_overlap[oldbin][newbin] += 1    
        # model params
        emissions = [NormalDistribution(m, st) for m, st in self.binmodel]
        transmat = bin_overlap / bin_overlap.sum(axis=1)[:, np.newaxis]
        transmat *= (1 - self.p_stay)
        lowerright = np.identity(self.nbins) * self.p_stay
        transmat = transmat + lowerright
        # sns.heatmap(transmat)
        # plt.gcf().set_dpi(200)
        binnums = np.zeros(self.nbins)
        for k in kmer_to_bin:
            binnums[k] += 1
        starts = binnums / binnums.sum()
        model = HiddenMarkovModel.from_matrix(transmat, emissions, starts, state_names=[str(idx) for idx in range(self.nbins)], verbose=True)
        return model
    def viz_hmm(self):
        statename_to_order = {'None-start' : -1, 'None-end' : len(self.hmm.states) - 1}
        order = np.argsort(list(map(lambda x: int(x) if x not in statename_to_order else statename_to_order[x], [m.name for m in self.hmm.states])))
        _, ax = plt.subplots()
        sns.heatmap(self.hmm.dense_transition_matrix()[order, :][:, order], ax=ax)
        return ax



class DeltaStayHMMBin(DeltaBin):
    def __init__(self, nbins=64, poremodel=model_6mer, fixed=True, bounds=None, clip=False, p_stay=0.5) -> None:
        super(DeltaStayHMMBin, self).__init__(nbins=nbins, poremodel=poremodel, fixed=fixed, bounds=bounds, clip=clip)
        self.p_stay = p_stay
        self.hmm = self._hmm_stay_model()
    def _hmm_stay_model(self):
        nbins = len(self.starts) + 1
        kmer_to_bin = np.searchsorted(self.starts, self.poremodel.means)
        bin_overlap = np.zeros((nbins, nbins))
        for mer in self.poremodel.KMER_STRS:
            oldbin = kmer_to_bin[self.poremodel.str_to_kmer(mer)]
            for new in 'ACGT':
                newmer = mer[1:] + new
                newbin = kmer_to_bin[self.poremodel.str_to_kmer(newmer)]
                bin_overlap[oldbin][newbin] += 1    
        # model params
        # model params
        emissions = [NormalDistribution(np.mean(self.poremodel.means[np.where(kmer_to_bin == idx)[0]]), 
                                    np.mean(self.poremodel.stdvs[np.where(kmer_to_bin == idx)[0]]))
                                    for idx in range(nbins)]
        emissions += [NormalDistribution(0, np.mean(self.poremodel.stdvs[np.where(kmer_to_bin == idx)[0]]))
                                    for idx in range(nbins)]

        transmat = bin_overlap / bin_overlap.sum(axis=1)[:, np.newaxis]
        transmat = transmat * (1-self.p_stay)
        lowerleft = transmat
        lowerright = np.identity(nbins) * self.p_stay
        topright = np.identity(nbins) * self.p_stay
        fullmat = np.concatenate([np.concatenate([transmat, topright], axis=1),np.concatenate([lowerleft, lowerright], axis=1)], axis=0)
        transmat = fullmat
        starts = np.concatenate([np.ones(nbins)/nbins, np.zeros(nbins)])
        model = HiddenMarkovModel.from_matrix(transmat, emissions, starts, 
                    state_names=[str(idx) for idx in range(nbins)] + ['stay_' + str(idx) for idx in range(nbins)], verbose=True)
        return model
    def bin_signal(self, signal, evdt=None, normalize=True):
        signal = np.array(signal)
        signal = self._preprocess(signal, evdt=evdt, normalize=normalize)
        if self.clip:
            signal = signal[(signal >= self.clipbounds[0]) & (signal <= self.clipbounds[1])]
        if self.bounds:
            signal = signal[(signal <= self.bounds[0]) | (signal >= self.bounds[1])]
        # bins = np.searchsorted(self.starts, signal)
        bins = self.hmm.viterbi(signal)[1]
        bins = np.array([int(state.name) for _, state in bins[1:] if not state.name.startswith('stay')])
        return bins[np.insert((np.where(np.diff(bins) != 0)[0] + 1), 0, 0)]
    def bin_sequence(self, seq):
        sig = seq_to_sig(self.poremodel, seq)
        if self.bounds:
            sig = sig[(sig <= self.bounds[0]) | (sig >= self.bounds[1])]
        # bins = np.searchsorted(self.starts, sig)
        lprob, bins = self.hmm.viterbi(sig)
        print(lprob)
        bins = np.array([int(state.name) for _, state in bins[1:] if not state.name.startswith('stay')])
        return bins[np.insert((np.where(np.diff(bins) != 0)[0] + 1), 0, 0)]
