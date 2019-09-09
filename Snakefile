# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

shell.prefix("set -o pipefail; ")
shell.prefix("set -e; ")
shell.prefix("set -u; ")
configfile: "config.json"

include_prefix = "rules/"

include:
    include_prefix + "bwa_mem_hic_wasp.rules"
include:    
    include_prefix + "haplotype.rules"

rule all:
    input:
        expand("mapping/out/{unit}.sorted.haplotyped.bam", unit=config["units"])
