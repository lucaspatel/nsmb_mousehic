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
include:
    include_prefix + "statistics.rules",

wildcard_constraints:
   unit = '\w+',
   bin = '(?!all)\w+'

rule all:
    input: 
        expand("mapping/out/{unit}.all.mcool", unit=config["units"]),
        expand("mapping/out/{unit}.{bin}.mcool", unit=config["units"], bin=config["bins"]),
        expand("mapping/out/{unit}.all.pairing.bedgraph", unit=config["units"]),
        expand("mapping/out/{unit}.{bin}.pairing.bedgraph", unit=config["units"], bin=config["bins"]),
        expand("mapping/out/{unit}.all.counts.zip", unit=config["units"]),
        expand("mapping/out/{unit}.{bin}.counts.zip", unit=config["units"], bin=config["bins"])
