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
        if type(self.cmds) == str:
            os.system(self.cmds)
        else:
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

def exe_tandom(cmds):
    for cmd in cmds:
        os.system(cmd)

def flexbar(threads, fq1, cell_id):
    cmd1 = "flexbar --adapter-min-overlap 7 -ae LEFT -as AAGCAGTGGTATCAACGCAGAGTACATGGG --min-read-length 16 --threads %s --reads %s --target %s.trimOligo" % (threads, fq1, cell_id)
    cmd2 = "flexbar --adapter-min-overlap 7 -ae RIGHT -a adapters_polyAnT.fa --min-read-length 16 --threads %s --reads %s.trimOligo.fastq --target %s.trimOligo.trimPolyA" % (threads, cell_id, cell_id)
    cmd3 = "rm %s.trimOligo.fastq" % cell_id
    return [cmd1, cmd2, cmd3]

def get_fq1(path):
    return sorted([fq for fq in os.listdir(path) if fq.endswith("_1.fastq.gz")])

def gzip(path):
    cmds = []
    for fq in os.listdir(path):
        cmds.append("gzip %s" %os.path.join(path, fq))
    return cmds


if __name__ == '__main__':
    # step1. trim switch oligo(left mode), polyA & adapter (right mode)
    os.path.exists(paras.target_path) or os.makedirs(paras.target_path)
    for fq in get_fq1(paras.input_path):
        fq1 = os.path.join(paras.input_path, fq)
        cell_id = os.path.join(paras.target_path, fq.split("_")[0])
        exe_tandom(flexbar(paras.threads, fq1, cell_id))
    # step2. gzip compress
    exe_parallel(gzip(paras.target_path), paras.threads)