The unknown passed variants could be found:

```
awk ' BEGIN {FS="\t"} {if ($7 == "PASS" && $3 == ".") {print $0}}' variants/NA12878.rmdup.realign.hc.filter.snpeff.dbsnp.vcf
```

