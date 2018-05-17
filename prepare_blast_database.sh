#!/bin/bash

PATH_TO_BACTNOG_DB=/mirror/eggnog/eggnog_4.5/data/bactNOG/bactNOG.raw_algs.tar.gz

TMPDIR=/tmp/bactNOG
ENOG_LIST=$(pwd)/data/curated_34_enogs.txt
TAXID_LIST=$(pwd)/data/tax_ids_used.txt

TARGET_DIR=$(pwd)/data/databases
mkdir -p ${TARGET_DIR} $TMPDIR

cd $TMPDIR

DATABASE_EXISTING=1
for ENOG in $(less ${ENOG_LIST}); do
	if ! [ -f ${TARGET_DIR}/$ENOG".fa" ]; then
		DATABASE_EXISTING=0
	fi
done


if [ ${DATABASE_EXISTING} -eq 0 ]; then
	tar -xf ${PATH_TO_BACTNOG_DB} -C .
fi

module load ncbiblastplus

for ENOG in $(less ${ENOG_LIST}); do
	if [ ${DATABASE_EXISTING} -eq 0 ]; then
		cp bactNOG_raw_algs/*"."$ENOG"."* $ENOG".fa"
		sed 's/\..*$/\t/;s/^>/\t>/;' $ENOG".fa" | tr -d '\n\-' | tr "\t" "\n" > ${TARGET_DIR}/$ENOG".fa"
		#some taxids changed (merged) after eggnog release. sed to replace the headers (so they fit the taxonomy tree)
		sed -i 's/>245018$/>649756/;s/>439483$/>929558/;s/>469595$/>1639133/;s/>469596$/>100884/;s/>525366$/>1423814/;' ${TARGET_DIR}/$ENOG".fa"
	fi
	makeblastdb -in ${TARGET_DIR}/$ENOG".fa" -dbtype prot 
done

if ! [ -f ${TAXID_LIST} ]; then
	grep -o "[0-9][0-9]*" ${TARGET_DIR}/*".fa" | cut -f2 -d":" | sort | uniq > ${TAXID_LIST}
fi

rm -r $TMPDIR
