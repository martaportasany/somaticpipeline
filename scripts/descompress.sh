#!/bin/bash
mkdir -p {output.dir}
cp {input.vcf} {output.dir}
cp {input.tsv} {output.dir}
cd {output.dir}
gzip -d {input.vcf} {input.tsv}
