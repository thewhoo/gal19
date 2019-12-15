import subprocess
import difflib
from collections import defaultdict

files = ['g1_2', 'g1_4', 'g1_20', 'g1_4096', 'g2_1', 'g2_2', 'g2_3', 'g2_4', 'g2_4096', 'g3_27_w1_sl', 'g3_703_w1_sl', 'g3_4089_wd_nosl', 'g3_4394_wd_nosl_az']

tests = {
    #1 thread
    'g1.dot100' : ['a d 2 0', ' a d 2 1'],
    'g1.dot101' : ['a d 4 0', 'a d 4 1'],
    'g1.dot102' : ['a d 20 0', 'a d 20 1'],
    'g1.dot103' : ['a d 4096 0', 'a d 4096 1'],
    'g2_noloop.dot104' : ['a d 1 0', 'a d 1 1'],
    'g2_noloop.dot105': ['a d 3 0', 'a d 3 1'],
    'g2_noloop.dot106' : ['a d 6 0', 'a d 6 1'],
    'g2_noloop.dot107' : ['a d 4096 0', 'a d 4096 1'],
    'g3_dense_w1.dot108' : ['a d 27 0', 'a d 27 1'],
    'g3_dense_w1.dot109' : ['a d 703 0', 'a d 703 1'],
    'g3_dense_wd_nosl.dot110' : ['a d 4089 0', 'a d 4089 1'],
    'g3_dense_wd_nosl.dot111' : ['a z 4394 0', 'a z 4394 1'], #max for 5 threads - bad_allc otherwise
    'g3_dense_w1_nosl.dot112' : ['a x 15650 0', 'a x 15650 1'],
    #2 threads
    'g1.dot201' : ['a d 2 0', ' a d 2 2'],
    'g1.dot202' : ['a d 4 0', 'a d 4 2'],
    'g1.dot203' : ['a d 20 0', 'a d 20 2'],
    'g1.dot204' : ['a d 4096 0', 'a d 4096 2'],
    'g2_noloop.dot205' : ['a d 1 0', 'a d 1 2'],
    'g2_noloop.dot206': ['a d 3 0', 'a d 3 2'],
    'g2_noloop.dot207' : ['a d 6 0', 'a d 6 2'],
    'g2_noloop.dot208' : ['a d 4096 0', 'a d 4096 2'],
    'g3_dense_w1.dot209' : ['a d 27 0', 'a d 27 2'],
    'g3_dense_w1.dot210' : ['a d 703 0', 'a d 703 2'],
    'g3_dense_wd_nosl.dot211' : ['a d 4089 0', 'a d 4089 2'],
    'g3_dense_wd_nosl.dot212' : ['a z 4394 0', 'a z 4394 2'], #max for 5 threads - bad_allc otherwise
    'g3_dense_w1_nosl.dot213' : ['a x 15650 0', 'a x 15650 2'],
    #3 threads
    'g1.dot301' : ['a d 2 0', ' a d 2 3'],
    'g1.dot302' : ['a d 4 0', 'a d 4 3'],
    'g1.dot303' : ['a d 20 0', 'a d 20 3'],
    'g1.dot304' : ['a d 4096 0', 'a d 4096 3'],
    'g2_noloop.dot305' : ['a d 1 0', 'a d 1 3'],
    'g2_noloop.dot306': ['a d 3 0', 'a d 3 3'],
    'g2_noloop.dot307' : ['a d 6 0', 'a d 6 3'],
    'g2_noloop.dot308' : ['a d 4096 0', 'a d 4096 3'],
    'g3_dense_w1.dot309' : ['a d 27 0', 'a d 27 3'],
    'g3_dense_w1.dot310' : ['a d 703 0', 'a d 703 3'],
    'g3_dense_wd_nosl.dot311' : ['a d 4089 0', 'a d 4089 3'],
    'g3_dense_wd_nosl.dot312' : ['a z 4394 0', 'a z 4394 3'], #max for 5 threads - bad_allc otherwise
    'g3_dense_w1_nosl.dot313' : ['a x 15650 0', 'a x 15650 3'],
    #4 threads
    'g1.dot400' : ['a d 2 0', ' a d 2 4'],
    'g1.dot401' : ['a d 4 0', 'a d 4 4'],
    'g1.dot402' : ['a d 20 0', 'a d 20 4'],
    'g1.dot403' : ['a d 4096 0', 'a d 4096 4'],
    'g2_noloop.dot404' : ['a d 1 0', 'a d 1 4'],
    'g2_noloop.dot405': ['a d 3 0', 'a d 3 4'],
    'g2_noloop.dot406' : ['a d 6 0', 'a d 6 4'],
    'g2_noloop.dot407' : ['a d 4096 0', 'a d 4096 4'],
    'g3_dense_w1.dot408' : ['a d 27 0', 'a d 27 4'],
    'g3_dense_w1.dot409' : ['a d 703 0', 'a d 703 4'],
    'g3_dense_wd_nosl.dot410' : ['a d 4089 0', 'a d 4089 4'],
    'g3_dense_wd_nosl.dot411' : ['a z 4394 0', 'a z 4394 4'], #max for 5 threads - bad_allc otherwise
    'g3_dense_w1_nosl.dot412' : ['a x 15650 0', 'a x 15650 4'],
    #5 threads
    'g1.dot500' : ['a d 2 0', ' a d 2 5'],
    'g1.dot501' : ['a d 4 0', 'a d 4 5'],
    'g1.dot502' : ['a d 20 0', 'a d 20 5'],
    'g1.dot503' : ['a d 4096 0', 'a d 4096 5'],
    'g2_noloop.dot504' : ['a d 1 0', 'a d 1 5'],
    'g2_noloop.dot505': ['a d 3 0', 'a d 3 5'],
    'g2_noloop.dot506' : ['a d 6 0', 'a d 6 5'],
    'g2_noloop.dot507' : ['a d 4096 0', 'a d 4096 5'],
    'g3_dense_w1.dot508' : ['a d 27 0', 'a d 27 5'],
    'g3_dense_w1.dot509' : ['a d 703 0', 'a d 703 5'],
    'g3_dense_wd_nosl.dot510' : ['a d 4089 0', 'a d 4089 5'],
    'g3_dense_wd_nosl.dot511' : ['a z 4394 0', 'a z 4394 5'], #max for 5 threads - bad_allc otherwise
    'g3_dense_w1_nosl.dot512' : ['a x 15650 0', 'a x 15650 5'],
    #6 threads
    'g1.dot600' : ['a d 2 0', ' a d 2 6'],
    'g1.dot601' : ['a d 4 0', 'a d 4 6'],
    'g1.dot602' : ['a d 20 0', 'a d 20 6'],
    'g1.dot603' : ['a d 4096 0', 'a d 4096 6'],
    'g2_noloop.dot604' : ['a d 1 0', 'a d 1 6'],
    'g2_noloop.dot605': ['a d 3 0', 'a d 3 6'],
    'g2_noloop.dot606' : ['a d 6 0', 'a d 6 6'],
    'g2_noloop.dot607' : ['a d 4096 0', 'a d 4096 6'],
    'g3_dense_w1.dot608' : ['a d 27 0', 'a d 27 6'],
    'g3_dense_w1.dot609' : ['a d 703 0', 'a d 703 6'],
    'g3_dense_wd_nosl.dot610' : ['a d 4089 0', 'a d 4089 6'],
    'g3_dense_w1_nosl.dot611' : ['a x 15650 0', 'a x 15650 6'],
    #7 threads
    'g1.dot700' : ['a d 2 0', ' a d 2 7'],
    'g1.dot701' : ['a d 4 0', 'a d 4 7'],
    'g1.dot702' : ['a d 20 0', 'a d 20 7'],
    'g1.dot703' : ['a d 4096 0', 'a d 4096 7'],
    'g2_noloop.dot704' : ['a d 1 0', 'a d 1 7'],
    'g2_noloop.dot705': ['a d 3 0', 'a d 3 7'],
    'g2_noloop.dot706' : ['a d 6 0', 'a d 6 7'],
    'g2_noloop.dot707' : ['a d 4096 0', 'a d 4096 7'],
    'g3_dense_w1.dot708' : ['a d 27 0', 'a d 27 7'],
    'g3_dense_w1.dot709' : ['a d 703 0', 'a d 703 7'],
    'g3_dense_wd_nosl.dot710' : ['a d 4089 0', 'a d 4089 7'],
    'g3_dense_w1_nosl.dot711' : ['a x 15650 0', 'a x 15650 7'],
    #8 threads
    'g1.dot800' : ['a d 2 0', ' a d 2 8'],
    'g1.dot801' : ['a d 4 0', 'a d 4 8'],
    'g1.dot802' : ['a d 20 0', 'a d 20 8'],
    'g1.dot803' : ['a d 4096 0', 'a d 4096 8'],
    'g2_noloop.dot804' : ['a d 1 0', 'a d 1 8'],
    'g2_noloop.dot805': ['a d 3 0', 'a d 3 8'],
    'g2_noloop.dot806' : ['a d 6 0', 'a d 6 8'],
    'g2_noloop.dot807' : ['a d 4096 0', 'a d 4096 8'],
    'g3_dense_w1.dot808' : ['a d 27 0', 'a d 27 8'],
    'g3_dense_w1.dot809' : ['a d 703 0', 'a d 703 8'],
    'g3_dense_wd_nosl.dot810' : ['a d 4089 0', 'a d 4089 8'],
    'g3_dense_w1_nosl.dot811' : ['a x 15650 0', 'a x 15650 8'],
}

graph_dir = 'graphs/dot'
bin = './gal'
times = defaultdict(list)
for i in range(12):
    for graph, cmds in tests.items():
        outputs = []
        for c in cmds:
            outputs.append((str(subprocess.check_output(f"{bin} {graph_dir}/{graph[:-3]} {c}", shell=True), 'utf-8')).splitlines())
        #print name
        print(graph, cmds[0], '/', cmds[1])
        #print("\tseq: " + outputs[0][-1])
        #print("\tpar: " + outputs[1][-1])
        times[graph[:-3] + "seq " + cmds[0][:-2]].append(outputs[0][-1])
        times[graph[:-3] + "par(" + graph[-3]  + ") " + cmds[0][:-2]].append(outputs[0][-1])

for key, value in times.items():
    print(key)
    for time in value:
        print("\t" + time)