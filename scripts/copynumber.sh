#!/bin/bash
cd {output.dir}
samtools mpileup -E -f {input.ref} {input.tumor} {input.normal} | java -jar {input.VarScanpath} copynumber varScan --mpileup 1

