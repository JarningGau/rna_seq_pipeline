import shutil
import os
import sys
import re
import threading
import argparse
import pandas as pd


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
    cmds_list = make_parallel(cmd, config_dict["cpu"])
    for cmd_batch in cmds_list:
        for cmd in cmd_batch:
            t = runParallel(cmd)
            t.start()
        t.join()

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

def load_sam(config_dict):
    sams = get_files(config_dict["sam_path"], "sam")
    return sams

def load_bam(config_dict):
    bams = get_files(config_dict["bam_path"], "bam")
    return bams

def load_flagstat(config_dict):
    flagstats = get_files(config_dict["bam_path"], "flagstat")
    return flagstats

def load_count_file(config_dict):
    count_files = get_files(config_dict["count_path"], "count")
    return count_files

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
            cmds.append(hisat2(config_dict, sample_id, fq1, fq2, sam))
    return cmds

def hisat2(config_dict, sample_id, fq1, fq2, sam):
    # for unpaired reads | -U <r> |
    if fq2 == "NA":
        cmd = "hisat2 -p %s --rg-id=%s --rg PL:ILLUMINA -x %s --dta -U %s -S %s" %(
            config_dict["cpu"], sample_id, config_dict["hisat2index"], fq1, sam)
    # for paired reads | -1 <m1> -2 <m2> |
    else:
        cmd = "hisat2 -p %s --rg-id=%s --rg PL:ILLUMINA -x %s --dat -1 %s -2 %s -S %s" %(
            config_dict["cpu"], sample_id, config_dict["hisat2index"], fq1, fq2, sam)
    return cmd

def sam2bam(config_dict):
    cmds = []
    sams = load_sam(config_dict)
    if len(sams)*2 <= config_dict["cpu"]:
        threads = int(config_dict["cpu"]) / len(sams)
    else:
        threads = 1
    for sam in sams:
        sample_id = sam.split(".")[0]
        sam = os.path.join(config_dict["sam_path"], sam)
        bam = os.path.join(config_dict["bam_path"], sample_id+".sorted.bam")
        flagstat = os.path.join(config_dict["bam_path"], sample_id+".flagstat")
        cmd1 = "samtools sort -@ %s %s > %s" % (threads, sam, bam)
        cmd2 = "samtools index %s" % bam
        cmd3 = "samtools flagstat %s > %s" % (bam, flagstat)
        #cmds.append([cmd1, cmd2, cmd3])
        cmds.append([cmd3])
    return cmds

def counter(config_dict):
    bams = load_bam(config_dict)
    cmds = []
    if config_dict["counter"] == "htseq-count":
        for bam in bams:
            sample_id = bam.split(".")[0]
            bam = os.path.join(config_dict["bam_path"], bam)
            count_file = os.path.join(config_dict["count_path"], sample_id+".count")
            cmd = htseq(config_dict, bam, count_file)
            cmds.append(cmd)
    if config_dict["counter"] == "featureCounts":
        for bam in bams:
            sample_id = bam.split(".")[0]
            bam = os.path.join(config_dict["bam_path"], bam)
            count_file = os.path.join(config_dict["count_path"], sample_id+".count")
            cmd = feature_counts(config_dict, bam, count_file)
            cmds.append(cmd)
    return cmds

def htseq(config_dict, bam, count_file):
    # default settings
    # --mode union
    # --minaqual 10
    cmd = ["htseq-count --format bam --order pos --stranded no --type exon --idattr=gene_id --additional-attr gene_name --additional-attr gene_type %s %s > %s" %(bam, config_dict["gene_gtf"], count_file)]
    return cmd

def feature_counts(config_dict, bam, count_file):
    cmd = ["featureCounts --extraAttributes gene_name,gene_type -a %s -o %s %s" % (config_dict["gene_gtf"], count_file, bam)]
    return cmd

def export_matrix(config_dict, mat_file):
    count_files = load_count_file(config_dict)
    col_name,mat = [],[]
    if config_dict["counter"] == "htseq-count":
        for cf in count_files:
            sample_id = cf.split(".")[0]
            cf = os.path.join(config_dict["count_path"], cf)
            mat.append(pd.read_csv(cf, sep='\t', header=None, index_col=0)[3])
            col_name.append(sample_id)
        mat = pd.concat(mat, axis=1)
        mat = mat.drop(mat.index[len(mat)-5:])
    if config_dict["counter"] == "featureCounts":
        for cf in count_files:
            sample_id = cf.split(".")[0]
            cf = os.path.join(config_dict["count_path"], cf)
            selected_col = os.path.join(config_dict["bam_path"], sample_id+".sorted.bam")
            mat.append(pd.read_csv(cf, sep='\t', header=1, index_col=0)[selected_col])
            col_name.append(sample_id)
        mat = pd.concat(mat, axis=1)
    mat.columns = col_name
    mat.to_csv(mat_file, sep='\t')

def export_align_rates(config_dict, output):
    flagstats = load_flagstat(config_dict)
    mr = re.compile(r"([\d\.]+%)")
    fo = open(output, "w")
    if config_dict["fq2"] == "NA":
        for fs in flagstats:
            sample_id = fs.split(".")[0]
            fs = os.path.join(config_dict["bam_path"], fs)
            align_rates = mr.findall("".join(open(fs).readlines()))[0]
            align_rates = "%.4lf" % (float(align_rates.rstrip("%")) / 100)
            fo.write("%s\t%s\n" % (sample_id, align_rates))
    fo.close()

def export_assign_rates(config_dict, output):
    if config_dict["counter"] == "htseq-count":
        pass
    if config_dict["counter"] == "featureCounts":
        fo = open(output, "w")
        summaries = get_files(config_dict["count_path"], "summary")
        for summ in summaries:
            sample_id = summ.split(".")[0]
            summ = os.path.join(config_dict["count_path"], summ)
            pdfm = pd.read_csv(summ, sep="\t", header=None, skiprows=1, index_col=0)
            ass = int(pdfm.loc["Assigned"])
            total = sum(pdfm[1])
            ass_r = "%.4lf" % (float(ass) / total)
            fo.write("%s\t%s\n" % (sample_id, ass_r))


if __name__ == '__main__':
    config_dict = load_config()
    check_config(config_dict)
    write_config(config_dict)
    # step1. alignment
    for cmd in aligner(config_dict): os.system(cmd)
    # step2. sam to bam
    cmds = sam2bam(config_dict)
    exe_parallel(cmds, config_dict["cpu"])
    # step3. counting
    cmds = counter(config_dict)
    exe_parallel(cmds, config_dict["cpu"])
    # step4. export expression matrix
    mat_file = os.path.join(config_dict["current_path"], config_dict["matrix_file"]+".mtx")
    align_meta = os.path.join(config_dict["current_path"], config_dict["matrix_file"]+".alignrates")
    assign_meta = os.path.join(config_dict["current_path"], config_dict["matrix_file"]+".assignrates")
    export_matrix(config_dict, mat_file)
    export_align_rates(config_dict, align_meta)
    export_assign_rates(config_dict, assign_meta)
