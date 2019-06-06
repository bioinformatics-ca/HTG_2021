```
head -10  alignment/NA12878/NA12878.sorted.dup.recal.metric.alignment.tsv | tail -4 | cut -f7
``` 

The alignment rate is 99.7%


Usually, we consider: 
 
   - A good alignment if > 90%
   - Reference assembly issues if [75-90]%
   - Probably a mismatch between sample and reference if < 75 %
