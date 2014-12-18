from subprocess import Popen, PIPE, STDOUT
import pprint
import re


def convert_from_wuss(rna_string, preserve_gaps=False):
    replacements = {
                        ',': '.',
                        '_': '.',
                        '-': '.',
                        ':': '.',
                        '.': ''
                    }

    if preserve_gaps:
        replacements["."] = "-"

    nest_map = {"(": ")", "{": "}", "<": ">", "[": "]"}

    stack = []

    list_str = list(rna_string)
    for i, c in enumerate(list_str):
        if c in nest_map.keys():
            stack.append((i, c))
        elif c in nest_map.values():
            top = stack.pop()
            if nest_map[top[1]] == c:
                list_str[top[0]] = "("
                list_str[i] = ")"
            else:
                raise Exception("Invalid structure")
        elif c in replacements.keys():
            list_str[i] = replacements[c]

    return "".join(list_str)


def get_suboptimals(rna_string, energy_range=None):
    command = ["RNAsubopt"]
    if energy_range is not None:
        command.extend(["-e", str(energy_range)])
    p = Popen(command, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    stdout, _ = p.communicate(input=rna_string)
    suboptimal_lines = stdout.splitlines()[1:]
    suboptimal_structures = [l.split()[0] for l in suboptimal_lines]

    return suboptimal_structures


def bp_distance(s1, s2):
    input_str = s1 + "\n" + s2
    p = Popen(["RNAdistance", "-DP"], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    stdout, _ = p.communicate(input=input_str)
    numbers_in_stdout = re.findall("\d+", stdout)
    if len(numbers_in_stdout) == 0:
        return "Error"

    return int(numbers_in_stdout[0])




def main():
    input_rna = raw_input("Enter RNA sequence:")
    rfam_ss = raw_input("Enter rfam secondary structure:")

    # Remove gaps
    input_rna = input_rna.replace(".","").replace("-", "").replace("_", "")

    # Convert from WUSS to simple notation
    ref_ss = convert_from_wuss(rfam_ss)

    # Get a list of suboptimal secondary structures
    subopts = get_suboptimals(input_rna, energy_range=2)

    for i, struct in enumerate(subopts, 1):
        print str(i) + ":", struct
        print "Distance:", bp_distance(struct, ref_ss)



if __name__ == '__main__':
    main()