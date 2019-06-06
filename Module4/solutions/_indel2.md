Largest indels size can be computed by extending the previous command while printing the size of the indel:

```
grep -v "^#" variants/NA12878.rmdup.realign.hc.vcf | \
awk '{ if(length($4) != length($5)) { print sqrt((length($4) - length($5))^2) "\t"$0 } }' | \
sort -k1,1nr | head 
```

The indel sizes are computed as the absolute value of the difference in charater size between the 2 alleles.


The largest ones are 15bp long `chr1 17798774` and `chr1 17817280`. 




