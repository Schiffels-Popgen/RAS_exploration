BEGIN {FS=""; OFS=""}
{
  for(i = 1; i <= NF; i++)
  if($i == 0) {$i = 2} else {if($i == 2) $i = 0}
  print
}