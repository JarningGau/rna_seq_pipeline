import re
import os
import argparse
import threading

## Parameters
parser = argparse.ArgumentParser()

parser.add_argument('-n','--threads', action='store', dest='threads', type=int, help='threads')
parser.add_argument('-i', action='store', dest='input_path', help='fastq path')
parser.add_argument('-t','--target', action='store', dest='target_path', help='target path')

paras = parser.parse_args()


class runParallel(threading.Thread):
    def __init__(self, cmds):
        super(runParallel, self).__init__()
        self.cmds = cmds

    def run(self):
        for cmd in self.cmds:
            os.system(cmd)


def make_parallel(cmds, threads):
    '''
    Divide tasks into blocks for parallel running.
    Put the cmd in parallel into the same bundle.
    The bundle size equals the threads.
    '''
    threads = int(threads)
    cmd_list = []
    i,j = 0,0
    for cmd in cmds:
        if j == 0:
            cmd_list.append(list())
            i += 1
        cmd_list[i-1].append(cmd)
        j = (j+1) % threads
    return cmd_list

def exe_parallel(cmds, threads):
    cmds_list = make_parallel(cmds, threads)
    for cmd_batch in cmds_list:
        for cmd in cmd_batch:
            t = runParallel(cmd)
            t.start()
        t.join()

def load_fq(path):
    fq1 = sorted([f for f in os.listdir(path) if f.endswith("_1.fastq.gz") and 'unassigned' not in f])
    fq2 = sorted([f for f in os.listdir(path) if f.endswith("_2.fastq.gz") and 'unassigned' not in f])
    return zip(fq1, fq2)


def umi_tools_extract(in1, in2, out1, out2, log, err):
    '''
    If the barcode is at the 5'end of fq1, then -I fq1 --read2-in fq2.
    Otherwise, -I fq2 --read2-in fq1
    '''
    cmd = "umi_tools extract -I %s --bc-pattern=NNNNNNNN --read2-in %s -S %s --read2-out %s -L %s -E %s" % (in1, in2, out1, out2, log, err)
    return [cmd]


if __name__ == '__main__':
    fqs = load_fq(paras.input_path)
    os.path.exists(paras.target_path) or os.makedirs(paras.target_path)
    cell_id_re = re.compile(r"cell-\d+")
    cmds = []
    for fq1, fq2 in fqs:
        in1 = os.path.join(paras.input_path, fq1)
        in2 = os.path.join(paras.input_path, fq2)
        cell_id = cell_id_re.findall(fq1)[0]
        out1 = os.path.join(paras.target_path, cell_id+"_2.fastq.gz")
        out2 = os.path.join(paras.target_path, cell_id+"_1.fastq.gz")
        log = os.path.join(paras.target_path, cell_id+".log")
        err = os.path.join(paras.target_path, cell_id+".err")
        cmds.append(umi_tools_extract(in1, in2, out1, out2, log, err))
    exe_parallel(cmds, paras.threads)
