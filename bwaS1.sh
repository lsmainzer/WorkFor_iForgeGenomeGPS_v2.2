#!/bin/sh
#redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# != 11 ]
then
        MSG="parameter mismatch"
        echo -e "jobid:${PBS_JOBID}\nprogram=$0 stopped at line=$LINENO.\nReason=$MSG" | ssh iforge "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine""
        exit 1;
else
	set -x
	echo `date`
        scriptfile=$0
        aligndir=$1
        parms=$2
        ref=$3
        outputdir=$4
	outputfile=$5
        R=$6
        scriptdir=$7
        elog=$8
        olog=$9
        email=${10}
        qsubfile=${11}
        LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"
        parameters=$( echo $parms | tr "_" " " )

        ## checking quality scores to gather additional params
        qscores=$scriptdir/checkFastqQualityScores.pl
        ill2sanger=`perl $qscores $R 10000`
        exitcode=$?
        if [ exitcode -ne 0 ]
        then
            MSG="QC scores for $R not generated.  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
        if [ $ill2sanger -gt 65 ]
        then
           qual="-I"
        else
           qual=" "
        fi

        cd $outputdir
        $aligndir/bwa aln $parameters $qual $ref $R > $outputfile
        exitcode=$?
        if [ exitcode -ne 0 ]
        then
            MSG="bwa aln command failed on $R.  exitcode=$exitcode. alignment failed"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit $exitcode;
        fi
        if [ ! -s $outputfile ]
        then
            MSG="$outputfile aligned file not created. alignment failed"
	   echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | ssh iforge  "mailx -s '[Support #200] Mayo variant identification pipeline' "$redmine,$email""
            exit 1;
        fi
        echo `date`
fi