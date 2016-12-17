#!/usr/bin/env python3

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import multiprocessing as mulproc

def run_vmc(length):
    seed = np.random.randint(2**32 - 1)
    out = subprocess.check_output(['./vmc', '100000000', str(seed), '2', '2', '1', str(length)], universal_newlines=True)

    out.rstrip('\n')
    lines = out.split('\n')
    kinetic = lines[0].split(" ")[-2]
    potential = lines[1].split(" ")[-2]
    total = lines[2].split(" ")[-2]

    return [kinetic, potential, total]

lengths = np.linspace(0.40, 1.80, 29)
# kinetic = np.zeros_like(lengths)
# potential = np.zeros_like(lengths)
# total = np.zeros_like(lengths)

with mulproc.Pool(8) as p:
    outputs = np.array(p.map(run_vmc, lengths))

    plt.plot(lengths, outputs[:, 0], marker='.', label='kinetic')
    plt.plot(lengths, outputs[:, 1], marker='.', label='potential')
    plt.plot(lengths, outputs[:, 2], marker='.', label='total')
    plt.xlabel('b [fm]')
    plt.ylabel('$E [MeV]$')
    plt.legend()
    plt.show()
