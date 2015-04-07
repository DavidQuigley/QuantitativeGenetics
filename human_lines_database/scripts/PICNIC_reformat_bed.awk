#!/usr/bin/awk -f
BEGIN {print "IDENTIFIER\tcell_line\tminor_allele\tmajor_allele"} 
{ 
gsub("id:","",$2); 
gsub(",minor:","\t",$2); 
gsub(",major:","\t",$2); 
print $1"\t"$2
} 
END {}  
