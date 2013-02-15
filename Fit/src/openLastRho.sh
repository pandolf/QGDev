#!/bin/bash
echo "Parsing file $1" >&2

PTLAST=0
RHOLAST=0

COUNT=0;
PT1=0;
PT2=0;
RHO1=0;
RHO2=0;

TMPFILE=/tmp/amarini/$RANDOM.txt
#emulate tac: perl -e 'print reverse <>'
tac $1 |while read line; do
		(( COUNT++ ))
		(( m=COUNT % 4 ))
		[ $m == 1 ] && echo -n -e "/" >&2
		[ $m == 2 ] && echo -n -e "\b|" >&2
		[ $m == 3 ] && echo -n -e "\b\\" >&2
		[ $m == 0 ] && echo -n -e "\b." >&2
		[ $COUNT -gt 40 ] && echo -n -e "\r" >&2 && COUNT=0 && echo -n -e "-----------\r" >&2
	echo ${line} | grep '^\[' >/dev/null && echo $line && continue; 
	#echo b >&2
	echo ${line} | grep '^{' >/dev/null && echo $line && continue; 
	#echo c >&2
	PT1=$(echo ${line} | awk '{print $1}')
	PT2=$(echo ${line} | awk '{print $2}')
	RHO1=$(echo ${line} | awk '{print $3}')
	RHO2=$(echo ${line} | awk '{print $4}')
	#echo d >&2
	if	[ "${PT1}" != "${PTLAST}" ] && [ "${RHO1}" != "$RHOLAST" ] ; then
		{ echo $line | sed "s:^$PT1 $PT2 $RHO1 ${RHO2} :$PT1 $PT2 $RHO1 100 :";} 
		#{ echo $line | sed "s:\ ${RHO2}\ : 100 :";} 
	#	echo a $PT1 $PT2 $RHO1 $RHO2 -- $PTLAST $RHOLAST >&2
	else
		 echo ${line}
	#	echo b $PT1 $PT2 $RHO1 $RHO2 -- $PTLAST $RHOLAST >&2
	fi;
	
	#echo e >&2
	PTLAST=${PT1}
	RHOLAST=${RHO1}
done > $TMPFILE 
tac $TMPFILE 

exit 0;

