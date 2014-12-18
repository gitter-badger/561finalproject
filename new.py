from subprocess import Popen, PIPE, STDOUT
from collections import Counter


def cleaned(seq):
    return seq.replace("-", "")


def mutated(seq):
    for i in xrange(len(seq)):

        original = seq[i].lower()
        if original == "t":
            original = "u"

        for nuc in set("acug") - {original}:
            yield seq[:i] + nuc + seq[i+1:]


def get_suboptimals(rna_string, energy_range=None, stochBT=None):
    command = ["RNAsubopt"]
    if energy_range is not None:
        command.extend(["-e", str(energy_range)])
    elif stochBT is not None:
        command.extend(["-p", str(stochBT)])

    p = Popen(command, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    stdout, _ = p.communicate(input=rna_string)
    suboptimal_lines = stdout.splitlines()[1:]
    suboptimal_structures = [l.split()[0] for l in suboptimal_lines]

    return suboptimal_structures

def get_mfe(rna_string):
    command = ["RNAfold"]
    p = Popen(command, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    stdout, _ = p.communicate(input=rna_string)
    ssline = stdout.splitlines()[1]
    ss = ssline.split()[0]
    return ss

def get_ground_truth(seq_With_gaps, cons_struct):
    str_list = []
    review = []
    # Do the following:
    # If both ends are in cons_struct, keep them
    #
    # for nuc, struc in zip(seq_With_gaps, cons_struct):

from main import convert_from_wuss, bp_distance

def most_probable_from_mutants(rna_string, samples_per_mutant=100):
    counter = Counter()
    for m in mutated(rna_string):
        counter.update(get_suboptimals(rna_string, stochBT=samples_per_mutant))
    return counter

def most_probable_from_real(rna_string, stochBT=100):
    counter = Counter()
    counter.update(get_suboptimals(rna_string, stochBT=stochBT))
    return counter

cs = "ACCUUCGGCCACACCACUUUGAAUAAGCCUGAUCUCGCCUGAUCUCGUCUGAUCUCAGAAGUGAAACAAUGUUGGGCUUGGUUAGUACUUGGAUGGGAGACCGCCUGGGAAUAUCAAGUGCUGUAGGCA"

