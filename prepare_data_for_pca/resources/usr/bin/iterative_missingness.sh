#!/bin/bash

# Arguments:
# $1 = starting threshold (e.g. 90)
# $2 = ending threshold (e.g. 98)
# $3 = increment (e.g. 2)

root="ukb_maf"
plink="plink"

# Initial missingness threshold as proportion
aspercent=$(echo "$1 / 100" | bc -l)
genomind_1=$(echo "1 - $aspercent" | bc -l)

# Step 1: SNP filtering
$plink \
  --bfile ${root} \
  --geno $genomind_1 \
  --make-bed \
  --out ${root}_SNP$1

# Step 2: Sample filtering
$plink \
  --bfile ${root}_SNP$1 \
  --mind $genomind_1 \
  --make-bed \
  --out ${root}_sample$1.SNP$1

# Begin iterative loop
newstep=$(($1 + $3))

for i in $(seq $newstep $3 $2); do
  aspercent=$(echo "$i / 100" | bc -l)
  genomind=$(echo "1 - $aspercent" | bc -l)
  prefix=$(($i - $3))

  # SNP filtering at new threshold
  $plink \
    --bfile ${root}_sample${prefix}.SNP${prefix} \
    --geno $genomind \
    --make-bed \
    --out ${root}_sample${prefix}.SNP${i}

  # Sample filtering at new threshold
  $plink \
    --bfile ${root}_sample${prefix}.SNP${i} \
    --mind $genomind \
    --make-bed \
    --out ${root}_sample${i}.SNP${i}
done

# Final output after all filtering steps
$plink \
  --bfile ${root}_sample$2.SNP$2 \
  --make-bed \
  --out ukb_cleaned
