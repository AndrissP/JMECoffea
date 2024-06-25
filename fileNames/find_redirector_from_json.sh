#!/bin/bash

# Input: List of data files and JSON file path
in_file=$1
out_file=$2
JSON_FILE="site_endpoints.json"
JSON_URL="https://cmssst.web.cern.ch/cmssst/site_info/site_endpoints.json"
wget $JSON_URL -O $JSON_FILE
LINES=$(cat $in_file)

for FILE in $LINES
do
   # Find the site alias
   site_alias=$(dasgoclient -query="site file=$FILE" | grep T2 | head -n 1)
	echo site_alias: $site_alias
   if [[ -z "$site_alias" ]] ; then
     echo "Discarded file = " $FILE
   else
     echo searching for redirector

     # Find the redirector from the JSON file
     redirector=$(jq -r --arg site_alias "$site_alias" '.data[] | select(.site == $site_alias and .type == "XROOTD") | .endpoint' "$JSON_FILE")

     # Output the results
     echo "File: $FILE, Site Alias: $site_alias, Redirector: $redirector"
	  full_name=root://${redirector}/${FILE}
	  echo full_name: $full_name
	  echo $full_name >> $out_file
   fi
done

