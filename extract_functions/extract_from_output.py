from util.exceptions import SCQPSolutionNotFound, SCFNotConvergedNotPossibleToRunMP2

"""
functions extract something from the out cp2k files using re
Example of 7.1 gw cp2k output in test
"""
import re

def main():
    num_o_fun = extract_number_of_independent_orbital_function('test/out.out')
    print(f'number of independent orbital function: {num_o_fun}')

    # try:
    #     my_gw_energies = return_gw_energies('test/out_gw_no_solution_found.out')
    # except SCQPSolutionNotFound:
    #     print("I got the SCQPSolutionNotFoundError. This is a correct result")

    # TEST: capture the exeption: reason of failed run: "SCF Not Converged. Not Possible to Run MP2"
    print('\nTEST: I test that I capture the exception: "SCF Not Converged. Not Possible to Run MP2"')
    try:
        my_gw_energies = return_gw_energies('test/scf_not_converged_not_possible_mp2.out')
        print('TEST IS FAILED. Everythingg was extracted\n')
    except SCFNotConvergedNotPossibleToRunMP2:
        print("I got the SCFNotConvergedNotPossibleToRunMP2. This is a correct result \nTEST IS PASSED \n")
    except FileNotFoundError:
        print('TEST NOT PASSED because of FileNotFoundError')
    except:
        print('TEST NOT PASSED')

    # TEST: scf DFT
    print('this should work -->')
    my_gw_energies = return_gw_energies_advanced('test/out_scf_gw.out')

    print('this should fail -->')
    my_gw_energies = return_gw_energies_advanced('test/scf_not_converged_not_possible_mp2.out')

    # print('This is the last line')


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
    """
    at present, it goes through the all occurences of the GW energies that looks like that -->

  GW quasiparticle energies
  -------------------------

  The GW quasiparticle energies are calculated according to:
  E_GW = E_SCF + Z * ( Sigc(E_SCF) + Sigx - vxc )

  The energy unit of the following table is eV.  Sigc_fit is a very conservative
  estimate of the statistical error of the fitting.

            MO      E_SCF       Sigc   Sigc_fit   Sigx-vxc          Z       E_GW
    20 ( occ )    -12.173      2.843      0.000     -6.007      1.000    -15.338
    21 ( occ )    -11.506      3.444      0.000     -5.781      1.000    -13.844
    22 ( occ )    -11.117      3.766      0.000     -6.847      1.000    -14.198
    23 ( occ )    -10.526      1.933      0.000     -4.335      1.000    -12.928
    24 ( occ )     -9.587      3.129      0.000     -6.135      1.000    -12.594
    25 ( occ )     -8.926      2.603      0.000     -5.521      1.000    -11.844
    26 ( occ )     -8.763      2.142      0.000     -4.872      1.000    -11.494
    27 ( occ )     -7.568      2.576      0.000     -5.395      1.000    -10.388
    28 ( occ )     -6.978      0.832      0.000     -3.496      1.000     -9.641
    29 ( occ )     -6.765      1.888      0.000     -4.853      1.000     -9.730
    30 ( vir )     -4.257     -2.675      0.000      5.399      1.000     -1.533
    31 ( vir )     -1.921     -3.047      0.000      5.737      1.000      0.768
    32 ( vir )     -1.062     -1.351      0.000      3.130      1.000      0.717
    33 ( vir )     -0.209     -1.072      0.000      2.539      1.000      1.258
    34 ( vir )      0.196     -3.106      0.000      5.618      1.000      2.708
    35 ( vir )      0.527     -1.207      0.000      2.765      1.000      2.085
    36 ( vir )      0.608     -2.120      0.000      3.816      1.000      2.304
    37 ( vir )      0.726     -2.840      0.000      4.714      1.000      2.599
    38 ( vir )      0.740     -0.645      0.000      1.805      1.000      1.900
    39 ( vir )      1.372     -0.587      0.000      1.596      1.000      2.381

  GW HOMO-LUMO gap (eV)                                                     8.20

  PERFORMANCE| PDGEMM flop rate (Gflops / MPI rank):                       26.64

    <--
    and focus on the last last occurence of that, more specifically on the HOMO and LUMO in their initial order (DFT).
    In this case: -->
    29 ( occ )     -6.765      1.888      0.000     -4.853      1.000     -9.730
    30 ( vir )     -4.257     -2.675      0.000      5.399      1.000     -1.533
    <--
    and it takes HOMO as -6.765, LUMO as -4.257, occ as -9.73 and vir as -1.533.

    This is only what we need if HOMO/LUMO order is retained in vir/occ order and if there is only one GW iteration!
    Otherwise, one may need: HOMO LUMO from the first iteration (initial DFT value) and:
    occ/vir from the first (G0W0) or the last (ev-GW) iteration.
    """
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
        print('GW energies were not extracted.')
        _through_an_exception_if_the_string_is_found(string='Self-consistent quasi-particle solution not found',
                                                     exception=SCQPSolutionNotFound,
                                                     all_file=all_file)
        _through_an_exception_if_the_string_is_found(string='SCF not converged: not possible to run MP2',
                                                     exception=SCFNotConvergedNotPossibleToRunMP2,
                                                     all_file=all_file)
        print("gs energies were not found")  # if neither of those
        return None, None, None, None


def return_gw_energies_advanced(path_to_file):
    """


  GW quasiparticle energies
  -------------------------

  The GW quasiparticle energies are calculated according to:
  E_GW = E_SCF + Z * ( Sigc(E_SCF) + Sigx - vxc )

  The energy unit of the following table is eV.  Sigc_fit is a very conservative
  estimate of the statistical error of the fitting.

            MO      E_SCF       Sigc   Sigc_fit   Sigx-vxc          Z       E_GW
    20 ( occ )    -12.173      2.843      0.000     -6.007      1.000    -15.338
    21 ( occ )    -11.506      3.444      0.000     -5.781      1.000    -13.844
    22 ( occ )    -11.117      3.766      0.000     -6.847      1.000    -14.198
    23 ( occ )    -10.526      1.933      0.000     -4.335      1.000    -12.928
    24 ( occ )     -9.587      3.129      0.000     -6.135      1.000    -12.594
    25 ( occ )     -8.926      2.603      0.000     -5.521      1.000    -11.844
    26 ( occ )     -8.763      2.142      0.000     -4.872      1.000    -11.494
    27 ( occ )     -7.568      2.576      0.000     -5.395      1.000    -10.388
    28 ( occ )     -6.978      0.832      0.000     -3.496      1.000     -9.641
    29 ( occ )     -6.765      1.888      0.000     -4.853      1.000     -9.730
    30 ( vir )     -4.257     -2.675      0.000      5.399      1.000     -1.533
    31 ( vir )     -1.921     -3.047      0.000      5.737      1.000      0.768
    32 ( vir )     -1.062     -1.351      0.000      3.130      1.000      0.717
    33 ( vir )     -0.209     -1.072      0.000      2.539      1.000      1.258
    34 ( vir )      0.196     -3.106      0.000      5.618      1.000      2.708
    35 ( vir )      0.527     -1.207      0.000      2.765      1.000      2.085
    36 ( vir )      0.608     -2.120      0.000      3.816      1.000      2.304
    37 ( vir )      0.726     -2.840      0.000      4.714      1.000      2.599
    38 ( vir )      0.740     -0.645      0.000      1.805      1.000      1.900
    39 ( vir )      1.372     -0.587      0.000      1.596      1.000      2.381

  GW HOMO-LUMO gap (eV)                                                     8.20

  PERFORMANCE| PDGEMM flop rate (Gflops / MPI rank):                       26.64

    and it takes HOMO as -6.765, LUMO as -4.257, occ as -9.73 and vir as -1.533.

    return:


    """
    with open(path_to_file, "r") as fin:
        all_file = fin.read()

    myiter = iter(all_file.splitlines())
    reason = None
    header_list = None
    num_gw_iter = None

    occ, vir, homo, lumo, occ_0, vir_0, occ_scf, vir_scf = None, None, None, None, None, None, None, None

    regex = re.compile("^\s*GW quasiparticle energies")
    regex1 = re.compile("^\s*[0-9]+ \( occ \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")  # 12 (occ) <one number>
    reg_occ_groups = re.compile("^\s*([0-9]+) \( (occ) \) ((?:\s*[-+]?[0-9]*\.[0-9]*\s*)*)")
    reg_vir_groups = re.compile("^\s*([0-9]+) \( (vir) \) ((?:\s*[-+]?[0-9]*\.[0-9]*\s*)*)")

    regex2 = re.compile("^\s*[0-9]+ \( vir \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")
    # regex1 = re.compile("^\s*[0-9]+\s+\(")

    line_occ = None
    line_vir = None

    num_gw_iter = 0
    # from collections import namedtuple
    import pandas as pd

    header_is_extracted = False
    list_of_gw_output = []
    while True:
        line = next(myiter, -1)
        if line == -1:
            break
        if regex.match(line):  # GW header
            for i in range(0, 8):  # 10 is the number of lines to the actual data
                line = next(myiter)
            if not header_is_extracted:
                header_txt = next(myiter)  # header of the table
                header_list = header_txt.split()
                header_list.insert(1, *['occ_or_vir'])  # type of the orbital: virtual or occupied
                header_list.append(*['num_gw_iter'])  # type of the orbital: virtual or occupied
                # print(f'this is the header: {header_list}')
                header_is_extracted = True
            else:
                next(myiter)  # all headers are the same
            line = next(myiter)
            occ_line = []
            vir_line = []
            while regex1.match(line):  # first group (occ)
                tmp = re.search(reg_occ_groups, line).groups()
                list_of_str = ' '.join(tmp).split()
                single_line = [int(list_of_str[0]), list_of_str[1], *[float(x) for x in list_of_str[2:]], num_gw_iter]
                occ_line.append(single_line)
                line = next(myiter)
                # now we check if this line is not vir
                while regex2.match(line):  # match virs
                    tmp = re.search(reg_vir_groups, line).groups()
                    list_of_str = ' '.join(tmp).split()
                    single_line = [int(list_of_str[0]), list_of_str[1], *[float(x) for x in list_of_str[2:]], num_gw_iter]
                    vir_line.append(single_line)
                    line = next(myiter)

            # list_of_gw_output.append(occ_line)
            list_of_gw_output.extend([*occ_line, *vir_line])
            num_gw_iter += 1
    my_e_df = pd.DataFrame.from_records(data=list_of_gw_output, columns=header_list)
    print(my_e_df)

    # why there are two G0W0: because orbitals sometimes change their order
    # DFT part
    if header_is_extracted:
        homo = max(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'occ')]['E_SCF'])
        lumo = min(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'vir')]['E_SCF'])

        # G0W0 part
        occ_0 = max(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'occ')]['E_GW'])
        vir_0 = min(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'vir')]['E_GW'])

        # GW part
        occ_scf = max(my_e_df.loc[(my_e_df['num_gw_iter'] == num_gw_iter-1) & (my_e_df['occ_or_vir'] == 'occ')]['E_GW'])
        vir_scf = min(my_e_df.loc[(my_e_df['num_gw_iter'] == num_gw_iter-1) & (my_e_df['occ_or_vir'] == 'vir')]['E_GW'])

        # wrong G0W0 part (correct if order of HOMO/LUMO and occ/vir are same)
        occ = list(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'occ')]['E_GW'])[-1]
        vir = list(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'vir')]['E_GW'])[0]
    else:
        pass

    if isinstance(occ, float) and isinstance(vir, float) and isinstance(homo, float) and isinstance(lumo, float) :
        print(f"G0W0 HOMO wrong: {occ} eV")
        print(f"G0W0 LUMO wrong: {vir} eV")
        print(f"GW HOMO scf: {occ_scf} eV")
        print(f"GW LUMO scf: {vir_scf} eV")
        print(f"G0W0 HOMO: {occ_0} eV")
        print(f"G0W0 LUMO: {vir_0} eV")
        print(f"HOMO: {homo} eV")
        print(f"LUMO: {lumo} eV")
        return occ, vir, homo, lumo, occ_scf, vir_scf, occ_0, vir_0
    else:
        print('GW energies were not extracted.')
        _through_an_exception_if_the_string_is_found(string='Self-consistent quasi-particle solution not found',
                                                     exception=SCQPSolutionNotFound,
                                                     all_file=all_file)
        _through_an_exception_if_the_string_is_found(string='SCF not converged: not possible to run MP2',
                                                     exception=SCFNotConvergedNotPossibleToRunMP2,
                                                     all_file=all_file)


def _through_an_exception_if_the_string_is_found(string, exception, all_file):
    print(f"I will check, if the reason is: {string}")
    str_to_find = f"^.*{string}.*$"
    regex = re.compile(str_to_find)
    # with open(path_to_file, "r") as fin:
    #     all_file = fin.read()
    myiter = iter(all_file.splitlines())
    while True:
        line = next(myiter, -1)
        if line == -1:
            break
        if regex.match(line):
            print(f'The run crashes because indeed: {string}')
            raise exception
    print("No, this is not the reason indicated above. ")


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