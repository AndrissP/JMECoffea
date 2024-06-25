# set the cvmfs environment
source /cvmfs/sft-nightlies.cern.ch/lcg/views/dev4/latest/x86_64-centos7-gcc11-opt/setup.sh
export SSL_CERT_DIR='/etc/pki/tls/certs:/etc/grid-security/certificates'
in_file=$1
out_file=$2
LINES=$(cat $in_file)

for FILE in $LINES
do
   site_alias=$(dasgoclient -query="site file=$FILE" | grep T2 | head -n 1)
	if [[ -z "$site_alias" ]] ; then
     echo "Discarded file = " $FILE
	else
		echo file = ${FILE}, site = ${site_alias}
		#echo full name =  $full_name
		echo ${FILE} ${site_alias} >> $out_file
	fi
done

