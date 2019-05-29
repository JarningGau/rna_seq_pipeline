import os
import threading
import argparse

## Parameters
parser = argparse.ArgumentParser()

parser.add_argument('--whitelist', action='store_const', dest='mode',
        const='whitelist',
        help='prepare the whitelist of cell barcodes')
parser.add_argument('-i', action='store', dest='input_path',
        help='path for where fastqs loacate')
parser.add_argument('--cell_num', action='store', dest='cell_num_file', default=None,
        help='file records cell counts per sample: <cell_num><tab><sample_id>')
parser.add_argument('-n','--threads', action='store', dest='threads', type=int,
        help='threads')

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

def whitelist(sample_id, cell_num):
    # CCCCCCCCCCCCNNNNNNNN for drop seqs (12+8)
    if cell_num is None:
        cmd = "umi_tools whitelist \
        --stdin %s_2.fastq.gz \
        --bc-pattern=CCCCCCCCCCCCNNNNNNNN \
        --log2stderr >%s.whitelist" % (sample_id, sample_id)
    else:
        cmd = "umi_tools whitelist \
        --stdin %s_2.fastq.gz \
        --bc-pattern=CCCCCCCCCCCCNNNNNNNN \
        --set-cell-number %s \
        --log2stderr >%s.whitelist" % (sample_id, cell_num, sample_id)
    return cmd

def load_cell_count_hash():
    if paras.cell_num_file is None:
        return {}
    cell_count_h = {}
    with open(paras.cell_num_file) as fi:
        for line in fi:
            v,k = line.strip().split('\t')
            cell_count_h[k] = v
    return cell_count_h

def extract(sample_id):
    cmd = "umi_tools extract \
    --bc-pattern=CCCCCCCCCCCCNNNNNNNN \
    --stdin %s_2.fastq.gz \
    --stdout %s_2.extracted.fastq.gz \
    --read2-in %s_1.fastq.gz \
    --read2-out=%s_1.extracted.fastq.gz \
    --filter-cell-barcode \
    --error-correct-cell \
    --whitelist=%s.whitelist" % (sample_id, sample_id, sample_id, sample_id, sample_id)
    return cmd

def load_fqs():
    fq1s = sorted([f for f in os.listdir(paras.input_path) if f.endswith("_1.fastq.gz")])
    fq2s = sorted([f for f in os.listdir(paras.input_path) if f.endswith("_2.fastq.gz")])
    return zip(fq1s, fq2s)


if __name__ == '__main__':
    if paras.mode == "whitelist":
        cmds = []
        cell_count_h = load_cell_count_hash()
        for fq1, fq2 in load_fqs():
            sample_id = fq1.split("_")[0]
            cell_num = cell_count_h.get(sample_id, None)
            cmds.append(whitelist(os.path.join(paras.input_path, sample_id), cell_num))
        exe_parallel(cmds, paras.threads)
    if paras.mode == "extract":
        cmds = []
        for fq1, fq2 in load_fqs():
            sample_id = os.path.join(paras.input_path, fq1.split("_")[0])
            cmds.append(extract(sample_id))
        exe_parallel(cmds, paras.threads)