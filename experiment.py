import pprint
from projection import project
from new import cleaned, most_probable_from_mutants, most_probable_from_real, get_mfe
from main import bp_distance
from random import randint

f = open("5s_seed_stockholm.txt")
true_5s = "..((((((.((.(.,.,,.,<<....-<.<.<.<.<...--..........-.<<.-..........-..<<<....<.<<....______.......>>.-->>>.>-.->>---.-.................>>>>>--.>.>..<...<<..-<<---..-<-<<----..-<<____>.>.-----.>>-...>..-->>-....>>.>..).))..).))))):"

seqs = {}
for line in f:
    if line[0] == "#" or len(line) < 20:
        continue
    else:
        name, sequence = line.split()
        seqs[name] = sequence

f.close()

fname = "result_" + str(randint(1, 10000000))
print "Storing results in ", fname
with open(fname, 'w') as f:

    for one_item in seqs.items():

        seq_name= one_item[0]
        print "Running experiment on ", seq_name
        gapped_sequence = one_item[1]
        ungapped_sequnce = gapped_sequence.replace("-", "")
        projected_real_structure = project(gapped_sequence, true_5s)
        ungapped_real_structure = projected_real_structure.replace("-", "")
        assert len(ungapped_sequnce) == len(ungapped_real_structure)

        mfe_prediction = get_mfe(cleaned(gapped_sequence))
        assert len(mfe_prediction) == len(ungapped_real_structure)
        mfe_score = bp_distance(mfe_prediction, ungapped_real_structure)
        print "MFE distance: ", mfe_score

        classical_prediction = most_probable_from_real(cleaned(gapped_sequence)).most_common()[0][0]
        assert len(classical_prediction) == len(ungapped_real_structure)
        classical_score = bp_distance(classical_prediction, ungapped_real_structure)
        print "Classical BP distance: ", classical_score

        mutant_prediction = most_probable_from_mutants(cleaned(gapped_sequence)).most_common()[0][0]
        assert len(mutant_prediction) == len(ungapped_real_structure)
        mutant_score = bp_distance(mutant_prediction, ungapped_real_structure)
        print "Mutant BP distance: ", mutant_score

        print "Projected consensus structure/ mfe structure / Classical secondary struct/ Mutant secondary struct:"
        print ungapped_real_structure
        print mfe_prediction
        print classical_prediction
        print mutant_prediction

        fstring = seq_name + " " + str(mfe_score) + " " + str(classical_score) + " " + str(mutant_score) + "\n"
        f.write(fstring)
        f.flush()



