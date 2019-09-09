#!/bin/bash

#Author: Yunjiang Qiu <serein927@gmail.com>
#File: run.sh
#Create Date: 2015-09-22 14:19:36

snakemake -p --rerun-incomplete -k -j 128 --cluster "qsub -l {cluster.ppn} -l {cluster.time} -N {params.jobname} -q {cluster.queue} -o pbslog/{params.jobname}.pbs.out -e pbslog/{params.jobname}.pbs.err" --jobscript jobscript.pbs --jobname "{rulename}.{jobid}.pbs" --cluster-config cluster.json 2> run.log
