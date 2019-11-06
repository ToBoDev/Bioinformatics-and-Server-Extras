#!/bin/bash
PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
# check and send ip address to email

#run this first to match correct grep pattern for your system (i.e. 'iet addr:<IP first 3>', 'eno1', eth1', etc)
MYIP=`ifconfig | grep 'inet addr:' | awk '{print $2}' | cut -d ":" -f2`;
TIME=`date`;

LASTIPFILE='/.last_ip_addr';
LASTIP=`cat ${LASTIPFILE}`

ipadd=$(cat <<EOF 
{"text":"New IP address for $HOSTNAME is $MYIP on $TIME Have a great day!"} 
EOF
)

if [[ ${MYIP} != ${LASTIP} ]]
then
	curl -X POST -H 'Content-type: application/json' --data "$ipadd" <INSERT SLACK WEBHOOK URL HERE - details:https://api.slack.com/web>;
	echo ${MYIP} > ${LASTIPFILE}
fi

#save this script somewhere on your server and add it as a cronjob! It should send a ping to your specified slack channel whenever you IP address changes (such as during power outages)
