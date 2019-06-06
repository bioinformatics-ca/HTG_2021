The precentage could be estimate directly by this command:

```
awk ' BEGIN {FS="\t" ; p=0; d=0 ; OFS="\t"} \
{if ($7 == "PASS") {p++; if ($3 != ".") {d++}}} \
END {print "PASS:", p , "dbSNP:" , d , "%:", d/p} ' variants/NA12878.rmdup.realign.hc.filter.snpeff.dbsnp.vcf
```

> PASS: &nbsp;&nbsp;&nbsp;&nbsp; 403 &nbsp;&nbsp;&nbsp;&nbsp; dbSNP: &nbsp;&nbsp;&nbsp;&nbsp; 380 &nbsp;&nbsp;&nbsp;&nbsp; %: &nbsp;&nbsp;&nbsp;&nbsp; 0.942928   

For those unfamiliar with `awk commands`, they can use these commands and then compute manually the percentage.


```
grep -v ^# variants/NA12878.rmdup.realign.hc.filter.snpeff.dbsnp.vcf | grep PASS | grep -c rs
grep -v ^# variants/NA12878.rmdup.realign.hc.filter.snpeff.dbsnp.vcf | grep -c PASS
```
> 379   
> 403   
