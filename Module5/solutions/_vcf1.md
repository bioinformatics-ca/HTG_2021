You can estimate the raw number of variant in each sample using the following command:

```
for i in SVvariants/*bcf ; do \
 echo $i ; \
 bcftools view $i | grep -v "^#" | wc -l  ;  \
done
```

You should get:

|NA12878|NA12891|NA12892|
|--|--|--|
|123|155|132|

