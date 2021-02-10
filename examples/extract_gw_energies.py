import re
from copy import deepcopy
from itertools import repeat


path_to_file = "my_run_folder/out2.out"

with open(path_to_file, "r") as fin:
    all_file = fin.read()

myiter = iter(all_file.splitlines())


regex = re.compile("^\s*GW quasiparticle energies")
regex1 = re.compile("^\s*[0-9]+ \( occ \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")
regex2 = re.compile("^\s*[0-9]+ \( vir \)(\s+[-+]?[0-9]*\.[0-9]*\s*)")
# regex1 = re.compile("^\s*[0-9]+\s+\(")

while True:
    line = next(myiter, -1)
    if line == -1:
        break
    # print(line)
    if regex.match(line):
        # print(line)
        for i in range(0,10):
            line = next(myiter)
        # print(line)
        if regex1.match(line):
                print("\nI have detected the first occ! The next line is potentailly vir!")
                # line_occ = deepcopy(line)
                line_vir = line
                while not regex2.match(line_vir):
                    print("This line is not yet virt")
                    line_occ = line_vir
                    line_vir = next(myiter)
                print("Finally! This is vit line!")
        print("Occupied:", line_occ),
        print("Virtual:", line_vir)
        # break
        # continue
