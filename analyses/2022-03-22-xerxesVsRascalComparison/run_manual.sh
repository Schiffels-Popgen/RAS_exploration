cat all_vars_m1.0_chr1.freqsum | awk '
BEGIN {
  sumRAS00 = 0.0
  norm = 0
}
NR > 1 {norm += 1}
NR > 1 && $5 + $6 +$7 +$8 >= 2 && $5 + $6 +$7 +$8 <= 5 {
  add = $5 * $5 / 4
  sumRAS00 += add
  # if (add > 0) print $1, $2, $3, $4, $5
}
END {
  print norm, sumRAS00 / norm
}'
