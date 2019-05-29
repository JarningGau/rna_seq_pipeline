import os
import argparse
import threading

## Parameters
parser = argparse.ArgumentParser()
parser.add_argument('-s','--run2sample', action='store', dest='run2sample',
        help='a file map run id (SRR) to sample id (GSM)')
parser.add_argument('-i', action='store', dest='input_path',
        help='path contains fastqs to merge')
parser.add_argument('-t', action='store', dest='target_path',
        help='path where merged fastq write in')
parser.add_argument('-n', action='store', dest='threads', type=int, 
        help='threads')

paras = parser.parse_args()


class runParallel(threading.Thread):
    def __init__(self, cmds):
        super(runParallel, self).__init__()
        self.cmds = cmds

    def run(self):
        if type(self.cmds) == str:
            #print(self.cmds)
            os.system(self.cmds)
        else:
            for cmd in self.cmds:
                #print(cmd)
                os.system(cmd)

def make_parallel(cmds, threads):
    '''
    Divide tasks into blocks for parallel running.
    Put the cmd in parallel into the same bundle.
    The bundle size equals the threads.
    '''
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
    cmd_list = make_parallel(cmds, threads)
    for cmd_batch in cmd_list:
        for cmd in cmd_batch:
            t = runParallel(cmd)
            t.start()
        t.join()

def merge_fqs(fqs, merged_fq):
    cmd = "zcat %s | gzip > %s" % (fqs, merged_fq)
    return cmd

def load_run2sample_hash():
    run2sample = {}
    for line in open(paras.run2sample):
        k,v = line.strip().split("\t")
        run2sample[k] = v
    return run2sample

def load_fq(path):
    fq1s = sorted([fq for fq in os.listdir(path) if fq.endswith("_1.fastq.gz")])
    fq2s = sorted([fq for fq in os.listdir(path) if fq.endswith("_2.fastq.gz")])
    return zip(fq1s,fq2s)

def sample2run_hash():
    run2sample = load_run2sample_hash()
    fqs = load_fq(paras.input_path)
    sample2run = {}
    for fq1, fq2 in fqs:
        run_id = fq1.split("_")[0]
        sample = run2sample[run_id]
        if sample not in sample2run:
            sample2run[sample] = {"fq1":[], "fq2":[]}
        sample2run[sample]["fq1"].append(os.path.join(paras.input_path, fq1))
        sample2run[sample]["fq2"].append(os.path.join(paras.input_path, fq2))
    return sample2run


if __name__ == '__main__':
    cmds = []
    os.path.exists(paras.target_path) or os.makedirs(paras.target_path)
    sample2run = sample2run_hash()
    for sample in sample2run.keys():
        merged_fq1 = os.path.join(paras.target_path, sample+"_1.fastq.gz")
        fq1s = " ".join(sample2run[sample]["fq1"])
        merged_fq2 = os.path.join(paras.target_path, sample+"_2.fastq.gz")
        fq2s = " ".join(sample2run[sample]["fq2"])
        cmds.append(merge_fqs(fq1s, merged_fq1))
        cmds.append(merge_fqs(fq2s, merged_fq2))
    exe_parallel(cmds, paras.threads)
