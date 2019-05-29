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

def exe_parallel(cmd, threads):
    cmds_list = make_parallel(cmd, threads)
    for cmd_batch in cmds_list:
        for cmd in cmd_batch:
            t = runParallel(cmd)
            t.start()
        t.join()

def dmp(fq1, fq2, target_path):
    # flexbar 3.4.0
    cmd = ["flexbar --threads 4 -z GZ -reads %s -reads2 %s -b barcode.fa -bt LTAIL -be 0.25 -bu -t %s" % (fq2, fq1, target_path)]
    return(cmd)

def load_fq(path):
    fq1 = sorted([f for f in os.listdir(path) if f.endswith("_1.fastq.gz")])
    fq2 = sorted([f for f in os.listdir(path) if f.endswith("_2.fastq.gz")])
    return zip(fq1, fq2)


if __name__ == '__main__':
    fqs = load_fq(paras.input_path)
    cmds = []
    for fq1, fq2 in fqs:
        srr_id = fq1.split("_")[0]
        _fq1 = os.path.join(paras.input_path, fq1)
        _fq2 = os.path.join(paras.input_path, fq2)
        _target_path = os.path.join(paras.target_path, srr_id) + "/"
        os.path.exists(_target_path) or os.makedirs(_target_path)
        cmds.append(dmp(_fq1, _fq2, _target_path))
    exe_parallel(cmds, paras.threads)
