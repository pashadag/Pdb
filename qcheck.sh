#/bin/sh
REF=$1
CONTIGS=$2
OUT=$CONTIGS

~/blat/blat -out=psl -fastMap -noHead $REF $CONTIGS $OUT.psl
#~/blat/blat -fastMap -minIdentity=99 $REF $CONTIGS $!.psl
cat $OUT.psl | item 10 2 11 16 17 | awk '{ if (($2 == 0) && (($5 - $4 + 0) == $3)) print $1 }' > $OUT.match
cat $OUT.match | sort -n | uniq -c | wc -l

