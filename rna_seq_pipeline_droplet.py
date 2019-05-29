import os
import sys
import threading
import shutil
import argparse


## Parameters
parser = argparse.ArgumentParser()

parser.add_argument('-c','--config', action='store', dest='config',
        help='config file')
parser.add_argument('--version', action='version', version='rna_seq_pipeline 0.1')

paras = parser.parse_args()


class runParallel(threading.Thread):
    def __init__(self, cmds):
        super(runParallel, self).__init__()
        self.cmds = cmds

    def run(self):
        if type(self.cmds) == str:
            os.system(self.cmds)
            #print(self.cmds)
        else:
            for cmd in self.cmds:
                os.system(cmd)
                #print(cmd)


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

def exe_tandem(cmds):
    for cmd in cmds: os.system(cmd)

def load_config():
    config_dict = {}
    config = [l.rstrip().split("=") for l in open(paras.config, 'r').readlines() if not l.startswith("#")]
    for k,v in config:
        config_dict[k] = v
    return config_dict

def mkdirs(path):
    os.path.exists(path) or os.makedirs(path)

def check_config(config_dict):
    if not os.path.exists(config_dict["fq_path"]):
        print("No such a fastq path: %s" % config_dict["fq_path"])
        sys.exit(-1)
    mkdirs(config_dict["log_path"])
    mkdirs(config_dict["sam_path"])
    mkdirs(config_dict["bam_path"])
    mkdirs(config_dict["count_path"])

def write_config(config_dict):
    shutil.copy(paras.config, config_dict["current_path"])

def get_files(path, appendix):
    return sorted([f for f in os.listdir(path) if f.endswith(appendix)])

def load_fastq(config_dict):
    if config_dict["fq2"] == "NA":
        fqs = get_files(config_dict["fq_path"], config_dict["fq1"])
    else:
        fq1s = get_files(config_dict["fq_path"], config_dict["fq1"])
        fq2s = get_files(config_dict["fq_path"], config_dict["fq2"])
        fqs = zip(fq1s, fq2s)
    return fqs

def hisat2(config_dict, sample_id, fq1, fq2, sam, bam):
    # for unpaired reads | -U <r> |
    if fq2 == "NA":
        cmd1 = "hisat2 -p %s --rg-id=%s --rg PL:ILLUMINA -x %s --dta -U %s -S %s" %(
            config_dict["cpu"], sample_id, config_dict["hisat2index"], fq1, sam)
    # for paired reads | -1 <m1> -2 <m2> |
    else:
        cmd1 = "hisat2 -p %s --rg-id=%s --rg PL:ILLUMINA -x %s --dta -1 %s -2 %s -S %s" %(
            config_dict["cpu"], sample_id, config_dict["hisat2index"], fq1, fq2, sam)
    cmd2 = "samtools sort -@ %s %s > %s" % (config_dict["cpu"], sam, bam)
    cmd3 = "rm %s" % sam
    return [cmd1, cmd2, cmd3]

def aligner(config_dict):
    fqs = load_fastq(config_dict)
    cmds = []
    if config_dict["aligner"] == "hisat2":
        for fq in fqs:
            if type(fq) == tuple:
                sample_id = fq[0].split(config_dict["sep"])[0]
                fq1 = os.path.join(config_dict["fq_path"], fq[0])
                fq2 = os.path.join(config_dict["fq_path"], fq[1])
            else:
                sample_id = fq.split(config_dict["sep"])[0]
                fq1 = os.path.join(config_dict["fq_path"], fq)
                fq2 = "NA"
            sam = os.path.join(config_dict["sam_path"], sample_id+".sam")
            bam = os.path.join(config_dict["bam_path"], sample_id+".sorted.bam")
            cmds.append(hisat2(config_dict, sample_id, fq1, fq2, sam, bam))
    return cmds

def assign(config_dict):
    cmds = []
    for bam in [f for f in os.listdir(config_dict["bam_path"]) if f.endswith("sorted.bam")]:
        sample_id = os.path.join(config_dict["bam_path"], bam.split(".")[0])
        cmd1 = "featureCounts -a %s -o %s.gene_assigned -R BAM %s.sorted.bam -T %s" % (config_dict["gene_gtf"], sample_id, sample_id, config_dict["cpu"])
        cmd2 = "samtools sort -@ %s %s.sorted.bam.featureCounts.bam > %s.sorted.assigned.bam" % (config_dict["cpu"], sample_id, sample_id)
        cmd3 = "samtools index %s.sorted.assigned.bam" % sample_id
        cmds.append([cmd1, cmd2, cmd3])
    return cmds

def count(config_dict):
    cmds = []
    for bam in [f for f in os.listdir(config_dict["bam_path"]) if f.endswith("sorted.assigned.bam")]:
        sample_id = bam.split(".")[0]
        cmd = "umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I %s.sorted.assigned.bam -S %s.counts.tsv.gz" % (os.path.join(config_dict["bam_path"], sample_id), os.path.join(config_dict["count_path"], sample_id))
        cmds.append(cmd)
    return cmds


if __name__ == '__main__':
    config_dict = load_config()
    check_config(config_dict)
    write_config(config_dict)
    # step1. alignment
    for cmds in aligner(config_dict): exe_tandem(cmds)
    # step2. assignment
    for cmds in assign(config_dict): exe_tandem(cmds)
    # step3. count
    exe_parallel(count(config_dict), config_dict["cpu"])