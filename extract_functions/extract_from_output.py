from util.exceptions import SCQPSolutionNotFound

"""
functions extract something from the out cp2k files using re
Example of 7.1 gw cp2k output in test
"""
import re

def main():
    num_o_fun = extract_number_of_independent_orbital_function('test/out.out')
    print(f'number of independent orbital function: {num_o_fun}')

    try:
        my_gw_energies = return_gw_energies('test/out_gw_no_solution_found.out')
    except SCQPSolutionNotFound:
        print("I got the SCQPSolutionNotFoundError. This is a correct result")


"""
    CRASH COURSE ON REGULAR EXPRESSION:
    
    ^       # Match start of string
    [-+]?   # Match a leading + or - (optional)
    [0-9]+  # Match one or more digit
    \.?     # Match a literal . (optional, escaped)
    [0-9]*  # Match zero or more digits
    $       # Match the end of the string

    Character classes
    . 	any character except newline
    \w \d \s 	word, digit, whitespace
    \W \D \S 	not word, digit, whitespace
    [abc] 	any of a, b, or c
    [^abc] 	not a, b, or c
    [a-g] 	character between a & g
    Anchors
    ^abc$ 	start / end of the string
    \b 	word boundary
    Escaped characters
    \. \* \\ 	escaped special characters
    \t \n \r 	tab, linefeed, carriage return
    \u00A9 	unicode escaped ©
    Groups & Lookaround
    (abc) 	capture group
    \1 	backreference to group #1
    (?:abc) 	non-capturing group
    (?=abc) 	positive lookahead
    (?!abc) 	negative lookahead
    Quantifiers & Alternation
    a* a+ a? 	0 or more, 1 or more, 0 or 1
    a{5} a{2,} 	exactly five, two or more
    a{1,3} 	between one & three
    a+? a{2,}? 	match as few as possible
    ab|cd 	match ab or cd
"""

def return_homo_lumo(path_to_file):

    """
    mo_cubes is supposed to output homos and lumos.
    this script is able to read them out from the path_to_file cp2k output file
    
    """


    with open(path_to_file,"r") as fin:
       all_file = fin.read()

    myiter = iter(all_file.splitlines())

    regex = re.compile("^\s*Eigenvalues of the occupied subspace spin\s+(.+)")
    regex_lumos = re.compile("^\s*Lowest eigenvalues of the unoccupied subspace spin\s+(.+)")

    f_begin = '(^\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'  # f = float. expr: any float
    f_middle = '(\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)*'
    f_end = '(\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$)'
    any_floating_point_numbers = re.compile(f_begin+f_middle+f_end)

    homos = []
    lumos = []

    while True:
        line = next(myiter)
        if regex.match(line):
            # print(line)
            line = next(myiter)
            # print(line)
            while True:
                line = next(myiter)
                if any_floating_point_numbers.match(line):
                    current_homos = line.split()
                    current_homos = [float(homo) for homo in current_homos]
                    homos.extend(current_homos)
                else:
                    break
        if regex_lumos.match(line):
            line = next(myiter)
            while True:
                line = next(myiter)
                if any_floating_point_numbers.match(line):
                    current_lumos = line.split()
                    current_lumos = [float(homo) for homo in current_lumos]
                    lumos.extend(current_lumos)
                else:
                    break
            break

    # print(f"All homos: {homos}")
    # print(f"All lumos: {lumos}")

    return homos, lumos


def return_gw_energies(path_to_file):
    with open(path_to_file, "r") as fin:
        all_file = fin.read()

    myiter = iter(all_file.splitlines())
    reason = None

    regex = re.compile("^\s*GW quasiparticle energies")
    regex1 = re.compile("^\s*[0-9]+ \( occ \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")
    regex2 = re.compile("^\s*[0-9]+ \( vir \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")
    # regex1 = re.compile("^\s*[0-9]+\s+\(")

    line_occ = None
    line_vir = None

    while True:
        line = next(myiter, -1)
        if line == -1:
            break
        if regex.match(line):
            for i in range(0,10):
                line = next(myiter)
            if regex1.match(line):
                # print("\nI have detected the first occ! The next line is potentailly vir!")
                line_vir = line
                while not regex2.match(line_vir):
                    # print("This line is not yet virt")
                    line_occ = line_vir
                    line_vir = next(myiter)
                # print("Finally! This is vit line!")
            # print("Occupied:", line_occ),
            # print("Virtual:", line_vir)

    if isinstance(line_occ, str) and isinstance(line_vir, str) and \
            line_occ.split().__len__() == 10 and line_vir.split().__len__() == 10:


        homo = float(line_occ.split()[4])
        lumo = float(line_vir.split()[4])

        occ = float(line_occ.split()[9])
        vir = float(line_vir.split()[9])

        # print(f"GW HOMO: {occ} eV")
        # print(f"GW LUMO: {vir} eV")
        return occ, vir, homo, lumo
    else:
        print("gw energies not extracted. I will check, if the reason is: Self-consistent quasi-particle solution not found")
        str_to_find = "^.*Self-consistent quasi-particle solution not found.*$"
        regex = re.compile(str_to_find)
        # with open(path_to_file, "r") as fin:
        #     all_file = fin.read()
        myiter = iter(all_file.splitlines())
        while True:
            line = next(myiter, -1)
            if line == -1:
                break
            if regex.match(line):
                print('The run crashes because the self-consistent quasi-particle solution not found')
                raise SCQPSolutionNotFound 

        # if True:  # here it checks if the scqp solution is not found. if yes, return the exception that will be processed accordingly in the calling script.
        #     raise SCQPSolutionNotFound
        # else:
        print("gs energies were not found")
        return None, None, None, None

def extract_total_energy(path_to_file):
    import re
    """
    returns total energy from file path_to_file
    """
    with open(path_to_file, "r") as fin:
        regex = re.compile(" ENERGY\| Total FORCE_EVAL \( QS \) energy \(a\.u\.\):\s+(.+)\n")
        for line in fin:
            match = regex.match(line)
            if match:
                return match.groups()[0]


def extract_number_of_independent_orbital_function(path_to_file):
    import re
    """
    returns total energy from file path_to_file
    """
    with open(path_to_file, "r") as fin:

        regex = re.compile(" Number of independent orbital functions:\s+([0-9]+)\n")
        for line in fin:
            match = regex.match(line)
            if match:
                return int(match.groups()[0])


if __name__ == '__main__':
    main()