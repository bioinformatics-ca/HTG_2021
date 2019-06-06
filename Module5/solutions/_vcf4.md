You can estimate the raw number of variant by type in each sample using the following command:

```
for i in SVvariants/*bcf ; do \
 echo $i ; \
 bcftools view $i | grep -v "^#" | awk ' {print $5} ' | sort | uniq -c  ;  \
done
```

You should get:

|Sample|Deletion|Duplication|Inversion|
|--|--|--|--|
|NA12878|92|15|16|
|NA12891|134|12|9|
|NA12892|113|11|8|
