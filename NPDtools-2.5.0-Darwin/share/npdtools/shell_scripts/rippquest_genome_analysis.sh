#!/bin/sh

sequence_fpath="$1"
compound_type=$2
mod_file="$3"

#pfam_file=$PWD/pfam_anal/${compound_type}_all.hmm
#pfam_file=pfam_anal/Pfam-A
interval_siz=10000
#min_contig_siz=10000
min_contig_siz=1000
offset=20000 # switch to 0 resulted in griseomycin loosing CypM and CypL 
max_orf_siz=100

BASEDIR=$(dirname $(dirname "$0"))
PYTHONBASE=$BASEDIR/python_libs/
SHELLBASE=$BASEDIR/shell_scripts/

hmmsearch_path=""
uname_str=`uname`
if [ "$uname_str" = "Linux" ]; then
  os=linux
elif [ "$uname_str" = "Darwin" ]; then
  os=osx
fi
hmmsearch_path="$BASEDIR/external_tools/hmmer/$os/hmmer/binaries/hmmsearch"

sed_cmd=sed

if [ -s "${sequence_fpath}".6_frame.format ]; then
	echo "#INFO : translated sequence ${sequence_fpath}.6_frame.format exists"
else
	echo "#INFO : translating the sequence"
	echo "#INFO : python $PYTHONBASE/transeq.py \"${sequence_fpath}\".format -o \"${sequence_fpath}\".6_frame.format --frame 6 --wide"
	#cat "${sequence_fpath}" | grep -v -e ">" | tr -d '\n' | ${sed_cmd} -e '1i>seq\' | ${sed_cmd} -e '$a\'  > "${sequence_fpath}".format
	cat "${sequence_fpath}" | ${sed_cmd} 's/>.*/	/' | tr -d '\n' | tr '\t' '\n' | tail -n +2 | ${sed_cmd} '$a\' | awk 'length>'$min_contig_siz | cat -n | awk '{print $1-1,"\t",$2}' | ${sed_cmd} 's/^/>ctg_/g' | tr -d ' ' | tr '\t' '\n' > "${sequence_fpath}".format
	if [ ! -s "${sequence_fpath}".format ]; then
	    echo "#ERROR: all contigs are too short to proceed!"
	    exit 1
	fi
	python $PYTHONBASE/transeq.py "${sequence_fpath}".format -o "${sequence_fpath}".6_frame.format --frame 6 --wide
	#transeq "${sequence_fpath}".format "${sequence_fpath}".6_frame -frame 6
	#cat "${sequence_fpath}".6_frame | fasta_formatter  > "${sequence_fpath}".6_frame.format
fi


list_enzyme=`cat $mod_file | grep "^enabler" | awk '{print $2}'`
list_single=`cat $mod_file | grep "^single" | awk '{print $7}'`
list_double=`cat $mod_file | grep "^double" | awk '{print $10}'`

ripp_hmm_search_dirname="${sequence_fpath}"_ripp_hmm_search
mkdir -p "${ripp_hmm_search_dirname}"
echo "#INFO : running hmmsearch"
echo $list_enzyme $list_single $list_double | tr ' ' '\n' | grep -v "^all$" | sort | uniq | while read pattern
do
	echo "#INFO : $hmmsearch_path --max -o \"${ripp_hmm_search_dirname}\"/${pattern}.log --domtblout \"${ripp_hmm_search_dirname}\"/${pattern}.txt $BASEDIR/pfam_anal/ripp_hmm/${pattern}.hmm \"${sequence_fpath}\".6_frame.format"
	$hmmsearch_path --max -o "${ripp_hmm_search_dirname}"/${pattern}.log --domtblout "${ripp_hmm_search_dirname}"/${pattern}.txt $BASEDIR/pfam_anal/ripp_hmm/${pattern}.hmm "${sequence_fpath}".6_frame.format
done

echo $list_enzyme | tr ' ' '\n' | ${sed_cmd} 's|^|"'"${ripp_hmm_search_dirname}"'"/|g' | ${sed_cmd} 's|$|.txt|g' | xargs cat | grep -v "#" | awk '$7<1e-05' | awk '{print "enabler",$1,$18,$19}' > "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals_aa

echo $list_single $list_double | tr ' ' '\n' | grep -v "^all$" | sort | uniq | while read pattern
do
	if [ ! -z $pattern ]; then
		echo $pattern | ${sed_cmd} 's|^|"'"${ripp_hmm_search_dirname}"'"/|g' | ${sed_cmd} 's|$|.txt|g' | xargs cat | grep -v "#" | awk '$7<1e-05' | awk '{print "'$pattern'",$1,$18,$19}' >> "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals_aa
	fi
done

cat "${sequence_fpath}".format | grep ">" > "${sequence_fpath}".format.contig_name
cat "${sequence_fpath}".format | grep -v ">" | awk '{print length}'  > "${sequence_fpath}".format.contig_siz
paste "${sequence_fpath}".format.contig_name "${sequence_fpath}".format.contig_siz > "${sequence_fpath}".format.contig_info

cat /dev/null > "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals
#genome_size=$(cat "${sequence_fpath}".format | tail -n 1 | wc -c)
cat "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals_aa | while read enzyme protein_id start_pos stop_pos
do
	frame=$(echo $protein_id | awk -F'_' '{print $NF}')
	contig=$(echo $protein_id | awk -F'_' '{print $1"_"$2}')
	#genome_size=$(grep -A1 "^>"$contig"$" ${sequence_fpath}.format | tail -n 1 | awk '{print length}')
	genome_size=$(grep "^>"$contig "${sequence_fpath}".format.contig_info | awk '{print $2}')
	start_nuc=$(echo $frame $start_pos $stop_pos $genome_size | awk '{if($1>3) print $4-$3*3; else print $2*3;}')
	stop_nuc=$(echo $frame $start_pos $stop_pos $genome_size | awk '{if($1>3) print $4-$2*3; else print $3*3;}')
	start_nuc_minus=$(echo $start_nuc $interval_siz | awk '{if($1-$2<=0) print 1; else print $1-$2;}')
	stop_nuc_plus=$(echo $stop_nuc $interval_siz $genome_size | awk '{if($1+$2>3*$3) print 3*$3; else print $1+$2;}')
	echo $start_nuc_minus $stop_nuc_plus | awk '{print "'$enzyme'","'$contig'",$1,$2}' >> "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals
done

contig_list=`cat ${sequence_fpath}.6_frame.format.${compound_type}.txt.intervals | awk '{print $2}' | sort | uniq`

cat /dev/null > "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged

echo $contig_list | tr ' ' '\n' | while read contig
do
	grep "^enabler" "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals | grep " "$contig" " | ${sed_cmd} 's/^enabler '$contig' //g' > "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.raw
	python $PYTHONBASE/merge_overlap_interval.py "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.raw | ${sed_cmd} 's/^/enabler '$contig' /g' >> "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged
done
grep -v "^enabler" "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals  >> "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged


cat /dev/null > "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.enzyme
cat "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged | grep -v "^enabler" | while read enzyme contig start_pos stop_pos
do
        med_pos=`echo $start_pos $stop_pos | awk '{print int(($1+$2)/2)}'`
        echo $enzyme" "$contig" "$med_pos >> "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.enzyme
done

cat "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged | grep "^enabler" | while read enzyme contig start_pos stop_pos
do
        printf '%b' "$enzyme $contig $start_pos $stop_pos "
        list=""
        while read enzyme_add ctg med_pos
        do
                #echo $enzyme_add $med_pos
                start_pos_off=$(($start_pos-$offset))
                stop_pos_off=$(($stop_pos+$offset))
                if [ "$med_pos" -ge "$start_pos_off" -a "$stop_pos_off" -ge "$med_pos" -a "$contig" = "$ctg" ]; then
                        #echo -n " "$enzyme_add" "$med_pos
                        list=$list" "$enzyme_add
                fi
        done < "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.enzyme
        list_uniq=`echo $list | tr ' ' '\n' | sort | uniq | tr '\n' ' '`
        echo $list_uniq | awk '{print NF}' | tr '\n' ' '
        echo $list_uniq
done > "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged.with_mod

cat "${sequence_fpath}".format "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged | grep "^enabler" > "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged.enabler
sh $SHELLBASE/extract_intervals.sh "${sequence_fpath}".format "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged.enabler > "${sequence_fpath}".format.${compound_type}.intervals.fasta

if [ -s "${sequence_fpath}".format.${compound_type}.intervals.fasta ]
then
        echo 'orf found'
        python $PYTHONBASE/getorf.py "${sequence_fpath}".format.${compound_type}.intervals.fasta -o "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format -w # between stop codons
        #prodigal -i "${sequence_fpath}".format.${compound_type}.intervals.fasta -a "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf
	#cat "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf | ${sed_cmd} 's/*$//g' | awk '{print $1}' | fasta_formatter > "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format
else
        echo 'no orf found'
        #cat /dev/null > "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format
        exit 1
fi

cat "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format  | tr '\n' '\t' | ${sed_cmd} 's/	>/ >/g' | ${sed_cmd} '$a\' | tr ' ' '\n' | awk '{print $1,substr($2, length($2)-'$max_orf_siz',length($2))}' | grep " .*M" | awk '{print $1}' > "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.header

cat "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format  | tr '\n' '\t' | ${sed_cmd} 's/	>/ >/g' | ${sed_cmd} '$a\' | tr ' ' '\n' | awk '{print $1,substr($2, length($2)-'$max_orf_siz',length($2))}' | grep " .*M" | awk '{print $2}' | ${sed_cmd} 's/^[^M]*M/M/g' > "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.seq

paste "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.header "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.seq | tr '\t' '\n' > "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.revised

cat "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.revised | awk 'NR % 2 == 1' | awk -F_ '{print $3}' > "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.idx

cat -n "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged.with_mod | awk '{printf $1-1" "; out="";for(i=6;i<=NF;i++){out=out"_"$i}; print out}' > "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged.with_mod.edit

awk 'FNR==NR{a[$1]=$2;next} {print a[$1]"\n"}' "${sequence_fpath}".6_frame.format.${compound_type}.txt.intervals.merged.with_mod.edit "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.idx > "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.add

paste "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.revised "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.add | tr -d '\t' | tr -d ' ' > "${sequence_fpath}".format.${compound_type}.intervals.fasta.orf.format.edit

echo "#DONE"
