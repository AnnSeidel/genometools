#!/bin/sh

# author: Stefan Kurtz, Nov. 2015

set -e -x

if test $# -eq 3
then
  minlen=$1
  minidentity=$2
  referencefile=$3
  queryfile=""
else
  if test $# -eq 4
  then
    minlen=$1
    minidentity=$2
    referencefile=$3
    queryfile=$4
  else
    echo "Usage: $0 <minlen> <minidentity> <referencefile> [queryfile]"
    exit 1
  fi
fi

reference=`basename $referencefile`
seedlength=14
maxfreq=21
gtbin=${GTINSTALL}/bin/gt

if test "${queryfile}" != ""
then
  query=`basename $queryfile`
  outputprefix="${reference}-${query}"
else
  query=${reference}
  outputprefix="${reference}"
fi

${gtbin} encseq encode -sds no -md5 no -des no -indexname ${reference} \
                                           ${referencefile}
if test "${queryfile}" != ""
then
  ${gtbin} encseq encode -sds no -md5 no -des no -indexname ${query} \
                                           ${query}
  qiioption="-qii ${query}"
else
  qiioption=""
fi

${gtbin} seed_extend -parts 4 -ii ${reference} ${qiioption} -t ${maxfreq} -l ${minlen} \
                    -seedlength 14 -diagbandwidth no -minidentity ${minidentity} \
                    -v -overlappingseeds -bias-parameters -no-reverse \
                    -history 60 > ${outputprefix}-se.matches
rm -f ${query}.ssp ${query}.fasta ${query}.esq
rm -f ${reference}.ssp ${reference}.fasta ${reference}.esq

# note that the following script randomly replaces wildcards by
# characters over the base alphabet. So for the same sequence and parameters
# the set of matches often varies over different runs.
scripts/convert2myersformat.rb $referencefile > ${reference}.fasta

if test $PACKAGES != ""
then
  MYERSPROG="$PACKAGES"/myers
else
  MYERSPROG="../myers"
fi

rm -f ${reference}.db
${MYERSPROG}/DAZZ_DB/fasta2DB ${reference}.db ${reference}.fasta
if test "${queryfile}" != ""
then
  scripts/convert2myersformat.rb ${queryfile} > ${query}.fasta
  ${MYERSPROG}/DAZZ_DB/fasta2DB ${query}.db ${query}.fasta
fi

echo "# Fields: s. len, s. seqnum, s. start, strand, q. len, q. seqnum, q. start, score, editdist, % identity" > ${outputprefix}-da.matches
scripts/rdj-spacepeak.sh -hashmark ${MYERSPROG}/DALIGNER/daligner -t${maxfreq} -I -A -Y -e0.${minidentity} \
                   -k${seedlength} -l${minlen} \
                   ${query}.db ${reference}.db >> ${outputprefix}-da.matches
rm -f ${reference}.db ${query}.db .${reference}.idx .${reference}.bps
rm -f ${query}.${reference}*.las .${query}.idx .${query}.bps
