I=0
bas=`basename "$1" | awk -F. '{print $1}' | sed 's/_//g'`
#bas=`basename $1`
cat "$2" | while read enzyme contig start stop
do
	#echo ">"$bas"_"$I"_"$start"_"$stop
	echo ">"$contig"_"$I"_"$start"_"$stop
	grep -A1 "^>$contig$" "$1" | tail -n 1 | cut -c $start-$stop
	I=$(($I+1))
done
