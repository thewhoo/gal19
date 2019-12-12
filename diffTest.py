import subprocess
import difflib

files = ['g1_2', 'g1_4', 'g1_20', 'g1_4096', 'g2_1', 'g2_2', 'g2_3', 'g2_4', 'g2_4096', 'g3_27_w1_sl', 'g3_703_w1_sl', 'g3_4089_wd_nosl', 'g3_4394_wd_nosl_az']

tests = {
    'g3_dense_w1_nosl.dot' : ['a x 12000 8', 'a x 12000 0'],
}

graph_dir = 'graphs/dot'
bin = './cmake-build-debug/gal19'
for graph, cmds in tests.items():
    outputs = []
    for c in cmds:
        outputs.append(str(subprocess.check_output(f"{bin} {graph_dir}/{graph} {c}", shell=True), 'utf-8'))
    for line in difflib.unified_diff(outputs[0], outputs[1], fromfile=cmds[0], tofile=cmds[1], lineterm=''):
        print(line)