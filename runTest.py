import subprocess
import difflib

files = ['g1_2', 'g1_4', 'g1_20', 'g1_4096', 'g2_1', 'g2_2', 'g2_3', 'g2_4', 'g2_4096', 'g3_27_w1_sl', 'g3_703_w1_sl', 'g3_4089_wd_nosl', 'g3_4394_wd_nosl_az']

for fileName in files:
    file = open(r"./graphs/run/" + fileName + ".run", "r")
    command = file.read()
    command = command[:-1]
    print("test " + command)
    output = subprocess.check_output(r"./gal ./graphs/dot/" + command, shell=True)
    output = str(output, 'utf-8')
    wantedOutputFile = open(r"./graphs/run/" + fileName + ".out", "r")
    wantedOutput = wantedOutputFile.read()
    for line in difflib.unified_diff(output, wantedOutput, fromfile='testOutput', tofile='wantedOutput', lineterm=''):
        print(line)
