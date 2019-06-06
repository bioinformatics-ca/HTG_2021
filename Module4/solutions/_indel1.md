The number of indel can be computed using this command:

```
grep -v "^#" variants/NA12878.rmdup.realign.hc.vcf | awk '{ if(length($4) != length($5)) { print $0 } }' | wc -l 
```

We have 78 INDELs in the realigned vcf



