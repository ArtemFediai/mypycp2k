def return_homo_lumo(path_to_file):

    import re
    """
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


    with open('out.out',"r") as fin:
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
            print(line)
            line = next(myiter)
            print(line)
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


    print(f"All homos: {homos}")
    print(f"All lumos: {lumos}")

    return homos, lumos