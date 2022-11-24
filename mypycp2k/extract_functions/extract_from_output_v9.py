from mypycp2k.util.exceptions import SCQPSolutionNotFound, SCFNotConvergedNotPossibleToRunMP2, NaNInGW, LargeSigc, \
    IterationLimit
import pandas as pd
import numpy as np
import re
from scipy.constants import physical_constants
Ha_to_eV = physical_constants['Hartree energy in eV'][0]

"""
functions extract something from the out cp2k files using re
Example of 7.1 gw cp2k output in test


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


def main():
    num_o_fun = extract_number_of_independent_orbital_function('test/out.out')
    my_print(f'number of independent orbital function: {num_o_fun}')

    # TEST: capture the exeption: reason of failed run: "SCF Not Converged. Not Possible to Run MP2"
    my_print('\nTEST: I test that I capture the exception: "SCF Not Converged. Not Possible to Run MP2"')
    try:
        my_gw_energies = return_gw_energies('test/scf_not_converged_not_possible_mp2.out')
        my_print('TEST IS FAILED. Everythingg was extracted\n')
    except SCFNotConvergedNotPossibleToRunMP2:
        my_print("I got the SCFNotConvergedNotPossibleToRunMP2. This is a correct result \nTEST IS PASSED \n")
    except FileNotFoundError:
        my_print('TEST NOT PASSED because of FileNotFoundError')
    except:
        my_print('TEST NOT PASSED')

    # TEST: scf DFT
    my_print('this should work -->')
    my_gw_energies = return_gw_energies_advanced('test/out_scf_gw.out')

    my_print('this should fail -->')
    try:
        my_gw_energies = return_gw_energies_advanced('test/scf_not_converged_not_possible_mp2.out')
    except SCFNotConvergedNotPossibleToRunMP2:
        my_print('the exception is captured. test passed')

    # print('This is the last line')
    # TEST: scf DFT
    my_print('this should crash because of nan->')
    try:
        my_gw_energies = return_gw_energies_advanced('test/pot_en_abnormal.out')
    except NaNInGW:
        my_print('this has crashed')
    my_print("<-- this should crash")

    my_print('this should crash because of large Sigc ->')
    try:
        my_gw_energies = return_gw_energies_advanced('test/false_converge.out')
    except LargeSigc:
        my_print('the exception of large Sigc is captured. The test is passed')
    my_print("<-- this should crash because of large Sigc")

    # return_gw_energies_advanced('test/false_converge.out')  # this will raise an exception!

    my_print('this should crash because of 20 ITERATIONS->')
    try:
        my_gw_energies = return_gw_energies_advanced('test/20_iterations.out')
    except IterationLimit:
        my_print('the exception of IterationLimit is captured. The test is passed')
    my_print("<-- this should crash because of 20 ITERATIONS")

    # my_gw_energies = return_gw_energies_advanced('test/false_converge.out')


def return_homo_lumo(path_to_file):
    """
    reads HOMOs and LUMOs from the output file.
    mo_cubes function has to be activated for homos and lumos to be output
    Units: Hartree.

    """

    with open(path_to_file, "r") as fin:
        all_file = fin.read()

    myiter = iter(all_file.splitlines())

    regex = re.compile(r"^\s*Eigenvalues of the occupied subspace spin\s+(.+)")
    regex_lumos = re.compile(r"^\s*Lowest Eigenvalues of the unoccupied subspace spin\s+(.+)")  # version 9.1

    f_begin = r'(^\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'  # f = float. expr: any float
    f_middle = r'(\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)*'
    f_end =    r'(\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$)'
    any_floating_point_numbers = re.compile(f_begin + f_middle + f_end)

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
            line = next(myiter)  # cp2k 9.1, OT
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

    regex = re.compile(r"^\s*GW quasiparticle energies")
    regex1 = re.compile(r"^\s*[0-9]+ \( occ \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")
    regex2 = re.compile(r"^\s*[0-9]+ \( vir \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")
    # regex1 = re.compile("^\s*[0-9]+\s+\(")

    line_occ = None
    line_vir = None

    while True:
        line = next(myiter, -1)
        if line == -1:
            break
        if regex.match(line):
            for i in range(0, 10):
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
        my_print('GW energies were not extracted.')
        _through_an_exception_if_the_string_is_found(string='Self-consistent quasi-particle solution not found',
                                                     exception=SCQPSolutionNotFound,
                                                     all_file=all_file)
        _through_an_exception_if_the_string_is_found(string='SCF not converged: not possible to run MP2',
                                                     exception=SCFNotConvergedNotPossibleToRunMP2,
                                                     all_file=all_file)
        my_print("gs energies were not found")  # if neither of those
        return None, None, None, None


def return_gw_energies_advanced(path_to_file, silent=False):
    """
    Returns HOMO / LUMO / GW HOMO / GW LUMO / G0W0 HOMO / G0W0 LUMO / "wrong" G0W0 HOMO and LUMO
    vocabulary:


  The GW quasiparticle energies are calculated according to:

                      E_GW = E_SCF + Sigc(E_GW) + Sigx - vxc

  Upper equation is solved self-consistently for E_GW, see Eq. (12) in J. Phys.
  Chem. Lett. 9, 306 (2018), doi: 10.1021/acs.jpclett.7b02740


  ------------
  G0W0 results
  ------------

  Molecular orbital   E_SCF (eV)       Sigc (eV)   Sigx-vxc (eV)       E_GW (eV)
     1 ( occ )          -510.160          -0.237         -44.087        -554.483
     2 ( occ )          -269.800          -4.974         -32.182        -306.956
     3 ( occ )           -25.231           4.214         -11.413         -32.430
     4 ( occ )           -16.731           2.500          -8.086         -22.317
     5 ( occ )           -12.336           2.173          -6.269         -16.432
     6 ( occ )           -10.427           2.071          -6.032         -14.387
     7 ( occ )           -10.127           1.711          -5.789         -14.205
     8 ( occ )            -7.372           2.089          -6.159         -11.442
     9 ( occ )            -5.597           2.519          -6.505          -9.584
    10 ( vir )             1.617          -1.040           4.960           5.537
    11 ( vir )             3.152          -1.099           4.879           6.931
    12 ( vir )             4.019          -1.224           5.258           8.053
    13 ( vir )             4.223          -1.205           5.163           8.181
    14 ( vir )             4.945          -1.723           6.136           9.357
    15 ( vir )            13.605          -3.512           8.156          18.249
    16 ( vir )            14.499          -1.890           6.662          19.271
    17 ( vir )            14.749          -2.079           6.581          19.251
    18 ( vir )            20.404          -2.065           8.375          26.713
    19 ( vir )            22.145          -4.760           8.775          26.161

  GW HOMO-LUMO gap (eV)                                                    15.12

  INTEG_INFO| Started with Integration point   1 of  16
  INTEG_INFO| Started with Integration point   2 of  16
  INTEG_INFO| Started with Integration point   3 of  16
  INTEG_INFO| Started with Integration point   4 of  16
  INTEG_INFO| Started with Integration point   5 of  16
  INTEG_INFO| Started with Integration point   6 of  16
  INTEG_INFO| Started with Integration point   7 of  16
  INTEG_INFO| Started with Integration point   8 of  16
  INTEG_INFO| Started with Integration point   9 of  16
  INTEG_INFO| Started with Integration point  10 of  16
  INTEG_INFO| Started with Integration point  11 of  16
  INTEG_INFO| Started with Integration point  12 of  16
  INTEG_INFO| Started with Integration point  13 of  16
  INTEG_INFO| Started with Integration point  14 of  16
  INTEG_INFO| Started with Integration point  15 of  16
  INTEG_INFO| Started with Integration point  16 of  16
  PERFORMANCE| DGEMM flop rate (Gflops / MPI rank):                        16.01


  ---------------------------------------
  Eigenvalue-selfconsistency cycle:    2
  ---------------------------------------

  Molecular orbital   E_SCF (eV)       Sigc (eV)   Sigx-vxc (eV)       E_GW (eV)
     1 ( occ )          -554.483          -0.371         -44.087        -554.617
     2 ( occ )          -306.956          -0.197         -32.182        -302.180
     3 ( occ )           -32.430           3.461         -11.413         -33.182
     4 ( occ )           -22.317           2.343          -8.086         -22.474
     5 ( occ )           -16.432           1.253          -6.269         -17.352
     6 ( occ )           -14.387           1.103          -6.032         -15.356
     7 ( occ )           -14.205           1.008          -5.789         -14.908
     8 ( occ )           -11.442           1.381          -6.159         -12.149
     9 ( occ )            -9.584           1.743          -6.505         -10.359
    10 ( vir )             5.537          -0.740           4.960           5.836
    11 ( vir )             6.931          -0.758           4.879           7.273
    12 ( vir )             8.053          -0.862           5.258           8.415
    13 ( vir )             8.181          -0.829           5.163           8.557
    14 ( vir )             9.357          -1.230           6.136           9.851
    15 ( vir )            18.249          -2.649           8.156          19.112
    16 ( vir )            19.271          -1.304           6.662          19.858
    17 ( vir )            19.251          -1.249           6.581          20.080
    18 ( vir )            26.713          -1.793           8.375          26.985
    19 ( vir )            26.161          -2.826           8.775          28.094

  GW HOMO-LUMO gap (eV)                                                    16.20
    @param silent:
    @type path_to_file: str
    @return  occ, vir, homo, lumo, occ_scf, vir_scf, occ_0, vir_0
    occ: G0W0 HOMO "wrong"
    vir: G0W0 LUMO "wrong"
    homo: dft HOMO
    lumo: dft LUMO
    occ_scf: ev-GW HOMO
    vir_scf: ev-GW LUMO
    occ_0: G0W0 HOMO
    vir_0: G0W0 LUMO
    """
    with open(path_to_file, "r") as fin:
        all_file = fin.read()

    myiter = iter(all_file.splitlines())
    reason = None
    header_list = None
    num_gw_iter = None

    occ, vir, homo, lumo, occ_0, vir_0, occ_scf, vir_scf = None, None, None, None, None, None, None, None

    regex_g0w0 = re.compile(r"^\s*G0W0 results")
    regex_gw = re.compile(r"^\s*Eigenvalue-selfconsistency cycle")
    regex_occ = re.compile(r"^\s*[0-9]+ \( occ \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")  # 12 (occ) <one number>
    regex_vir = re.compile(r"^\s*[0-9]+ \( vir \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")
    reg_occ_groups = re.compile(r"^\s*([0-9]+) \( (occ) \) ((?:\s*[-+]?[0-9]*\.[0-9]*\s*)*)")
    reg_vir_groups = re.compile(r"^\s*([0-9]+) \( (vir) \) ((?:\s*[-+]?[0-9]*\.[0-9]*\s*)*)")

    header_list_real = ['Molecular orbital', 'E_SCF', 'Sigc',   'Sigx-vxc', 'E_GW']  # this header will be used.

    # regex_occ = re.compile("^\s*[0-9]+\s+\(")

    num_gw_iter = 0

    header_is_extracted = False
    list_of_gw_output = []
    while True:
        line = next(myiter, -1)
        if line == -1:
            break
        # g0w0 comes first, but I will do it differently
        # if regex_g0w0.match(line):
        #     pass
        if regex_gw.match(line) or regex_g0w0.match(line):  # GW header or G0W0 header
            for i in range(0, 2):  # from the match line to the header of the GW table
                line = next(myiter)
            if not header_is_extracted:  # we need a header
                header_txt = next(myiter)  # header of the table
                # header_list = header_txt.split()  # this does not work in 9.1. I use predefined list
                header_list = header_list_real
                header_list.insert(1, *['occ_or_vir'])  # type of the orbital: virtual or occupied
                header_list.append(*['num_gw_iter'])  # type of the orbital: virtual or occupied
                # print(f'this is the header: {header_list}')
                if header_txt:
                    header_is_extracted = True
            else:
                next(myiter)  # all headers are the same
            line = next(myiter)
            occ_line = []
            vir_line = []
            while regex_occ.match(line):  # first group (occ)
                tmp = re.search(reg_occ_groups, line).groups()
                list_of_str = ' '.join(tmp).split()
                # if you check it here, you loose time for "if"
                single_line = [int(list_of_str[0]), list_of_str[1], *[float(x) for x in list_of_str[2:]], num_gw_iter]
                occ_line.append(single_line)
                line = next(myiter)
                # now we check if this line is not vir
                while regex_vir.match(line):  # match virs
                    tmp = re.search(reg_vir_groups, line).groups()
                    list_of_str = ' '.join(tmp).split()
                    single_line = [int(list_of_str[0]), list_of_str[1], *[float(x) for x in list_of_str[2:]],
                                   num_gw_iter]
                    vir_line.append(single_line)
                    line = next(myiter)

            # list_of_gw_output.append(occ_line)
            list_of_gw_output.extend([*occ_line, *vir_line])
            num_gw_iter += 1
    my_e_df = pd.DataFrame.from_records(data=list_of_gw_output, columns=header_list)
    # print(my_e_df)

    # why there are two G0W0: because orbitals sometimes change their order
    # DFT part
    if header_is_extracted:
        homo = max(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'occ')]['E_SCF'])
        lumo = min(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'vir')]['E_SCF'])

        # G0W0 part
        occ_0 = max(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'occ')]['E_GW'])
        vir_0 = min(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'vir')]['E_GW'])

        # GW part
        occ_scf = max(
            my_e_df.loc[(my_e_df['num_gw_iter'] == num_gw_iter - 1) & (my_e_df['occ_or_vir'] == 'occ')]['E_GW'],
            default="not extracted")
        vir_scf = min(
            my_e_df.loc[(my_e_df['num_gw_iter'] == num_gw_iter - 1) & (my_e_df['occ_or_vir'] == 'vir')]['E_GW'],
            default="not extracted")

        # wrong G0W0 part (only correct if order of HOMO/LUMO and occ/vir are same)
        occ = list(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'occ')]['E_GW'])[-1]
        vir = list(my_e_df.loc[(my_e_df['num_gw_iter'] == 0) & (my_e_df['occ_or_vir'] == 'vir')]['E_GW'])[0]

        if occ_scf == "not extracted":
            my_print(f'The run crashes because there is a NaN in the last scf iteration')
            raise NaNInGW
    else:
        pass

    # check for large numbers in the xc energy
    if header_is_extracted:
        max_Sigc = max(np.abs(my_e_df.loc[(my_e_df['num_gw_iter'] == num_gw_iter - 1)]['Sigc']))
        if not silent:
            my_print(f'maximum Sigc = {max_Sigc}')  # debug
        else:
            pass
        if max_Sigc > 20.0:  # todo: hard-coded!
            my_print('Unphysically large |Sigc| (>20.0 eV). I return LargeSigc expection')
            raise LargeSigc

    # check if 20 iter limit is reached
    if header_is_extracted:
        if max(my_e_df['num_gw_iter']) == 19:  # todo: hard-coded. assume, 20 iteration is the limit!
            my_print('GW scf not converged: 20 iterations reached! I return IterationLimit exception')
            raise IterationLimit

    if isinstance(occ, float) and isinstance(vir, float) and isinstance(homo, float) and isinstance(lumo, float):
        my_print(f"G0W0 HOMO wrong: {occ} eV")
        my_print(f"G0W0 LUMO wrong: {vir} eV")
        my_print(f"GW HOMO scf: {occ_scf} eV")
        my_print(f"GW LUMO scf: {vir_scf} eV")
        my_print(f"G0W0 HOMO: {occ_0} eV")
        my_print(f"G0W0 LUMO: {vir_0} eV")
        my_print(f"HOMO: {homo} eV")
        my_print(f"LUMO: {lumo} eV")
        return occ, vir, homo, lumo, occ_scf, vir_scf, occ_0, vir_0
    else:
        my_print('GW energies were not extracted.')
        _through_an_exception_if_the_string_is_found(string='Self-consistent quasi-particle solution not found',
                                                     exception=SCQPSolutionNotFound,
                                                     all_file=all_file)
        _through_an_exception_if_the_string_is_found(string='SCF not converged: not possible to run MP2',
                                                     exception=SCFNotConvergedNotPossibleToRunMP2,
                                                     all_file=all_file)


def _through_an_exception_if_the_string_is_found(string, exception, all_file):
    my_print(f"I will check, if the following bad string not in the output: {string}")
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
            my_print(f'The run crashes because: {string}')
            raise exception
    my_print(f"No, this is no string {string} in the output. ")


def _extract_total_energy(path_to_file):
    import re
    """
    returns total energy from file path_to_file [eV]
    """
    with open(path_to_file, "r") as fin:
        regex = re.compile(r" ENERGY\| Total FORCE_EVAL \( \w+ \) energy \[a\.u\.\]:\s+(.+)\n")
        for line in fin:
            match = regex.match(line)
            if match:
                return float(match.groups()[0]) * Ha_to_eV


def extract_esp_charges_accuracy(path_to_file):
    import re
    """
    returns accuracy of eps fit
    """

    f_begin = r'(^\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)'  # f = float. expr: any float
    f_middle = r'(\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)*'
    f_end =    r'(\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$)'

    regex1 = re.compile(r"^\s*[0-9]+ \( occ \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")

    with open(path_to_file, "r") as fin:
        regex_rms = re.compile(r"^\s*Root-mean-square \(RMS\) error of RESP fit:\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*")
        regex_rrms = re.compile(r"^\s*Relative root-mean-square \(RRMS\) error of RESP fit:\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*")
        for line in fin:
            match_rms = regex_rms.match(line)
            if match_rms:
                rms = float(match_rms.groups()[0])  # units are not clear * Ha_to_eV
            match_rrms = regex_rrms.match(line)
            if match_rrms:
                rrms = float(match_rrms.groups()[0])
    return rms, rrms

def extract_number_of_independent_orbital_function(path_to_file):
    import re
    """
    returns total energy from file path_to_file
    """
    with open(path_to_file, "r") as fin:

        regex = re.compile(r" Number of independent orbital functions:\s+([0-9]+)\n")
        for line in fin:
            match = regex.match(line)
            if match:
                return int(match.groups()[0])


# def my_print(x, *args, **kwargs):
#     my_print(x, *args, **kwargs)


# todo: this is hard coded. One has to change here silen manually. Set it to False to see the output in out files.

def my_print(x, silent=True, *args, **kwargs):
    if not silent:
        print(x, *args, **kwargs)
    else:
        pass


if __name__ == '__main__':
    main()
