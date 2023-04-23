import collections
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import cm
import pylab
import matplotlib
import math
import os
import sys
import numpy as np
# from Bio import SeqIO
import re


class FrequencyCGR:

    def __init__(self, sequence):
        self.sequence = sequence

    def count_kmers(self, sequence, k):
        d = collections.defaultdict(int)
        for i in range(len(sequence) - (k - 1)):
            d[sequence[i:i + k]] += 1
        for key in d.keys():
            if "N" in key:
                del d[key]
        return d

    def probabilities(self, sequence, kmer_count, k):
        probabilities = collections.defaultdict(float)
        N = len(sequence)
        for key, value in kmer_count.items():
            probabilities[key] = float(value) / (N - k + 1)
        return probabilities

    def chaos_game_representation(self, probabilities, k):
        array_size = int(math.sqrt(4 ** k))
        chaos = []
        for i in range(array_size):
            chaos.append([0] * array_size)

        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1
        for key, value in probabilities.items():
            for char in key:
                if char == "U":
                    posx += maxx / 2
                elif char == "C":
                    posy += maxy / 2
                elif char == "G":
                    posx += maxx / 2
                    posy += maxy / 2
                maxx = maxx / 2
                maxy /= 2

            chaos[int(posy - 1)][int(posx - 1)] = value
            maxx = array_size
            maxy = array_size
            posx = 1
            posy = 1

        return chaos

    def save_fcgr(self, k, path):
        chaos = self.chaos_game_representation(self.probabilities(str(self.sequence), self.count_kmers(str(self.sequence), k), k), k)

        # show with
        pylab.figure(figsize=(5, 5))
        pylab.imshow(chaos, cmap=matplotlib.colormaps.get_cmap('gray'))
        pylab.plot()
        pylab.axis('off')
        pylab.xticks([])
        pylab.yticks([])
        #pylab.show()
        pylab.savefig(path, dpi=100, bbox_inches='tight', pad_inches=0)
        pylab.close()

    def list_fcgr(self, sequence, k):
        chaos = self.chaos_game_representation(self.probabilities(str(sequence), self.count_kmers(str(sequence), k), k), k)

        chaos = np.array(chaos)
        fcgr = chaos.flatten()
        np.savetxt('work/fcgrlist/%s.txt' % id, fcgr, fmt='%f', delimiter=',')

    def matrix_fcgr(self, sequence, k):
        chaos = self.chaos_game_representation(self.probabilities(str(sequence), self.count_kmers(str(sequence), k), k), k)

        chaos = np.array(chaos)
        np.savetxt('work/matrix_6/%s.txt' % id, chaos, fmt='%f', delimiter=' ')
