#!/bin/bash

# Read input script
read -p "input script: " script0
script=${script0}.py
echo 'input script: '$script

## Generation Trial.py
# Extract lines we need
head_num=$(cat $script | awk '/Initialize\ parser/ {print FNR}')
tail_num=$(cat $script | awk '/UserInput\ =\ parser.parse_args\(\)/ {print FNR}')
echo 'first line: '$head_num
echo 'last line: '$tail_num

# Extract input letter
cat $script | grep inputs.add_argument > tmp1
sed -i '1d' tmp1
sed "s/',/='\ +/" tmp1 > tmp2
sed "s/inputs.add_argument('-/              ',/" tmp2 > tmp1

# Extract dest name
cat $script | grep dest\= > tmp2
sed "s/                    dest='/UserInput./" tmp2 > tmp3
sed "s/',/\ +/" tmp3 > tmp2
paste -d\  tmp1 tmp2 > tmp3

# Replacing with the input file name
sed "s/XXX/$script0/" /home/wud18/python/Trial_maker/Trial.tail1 > tmp1
sed "s/XXX/$script0/" /home/wud18/python/Trial_maker/Trial.tail2 > tmp2

# Combining all lines
head -n $tail_num $script > tmp0
for((i=1;i<$head_num;i++)); do sed -i '1d' tmp0; done
cat /home/wud18/python/Trial_maker/Trial.head tmp0 tmp1 tmp3 tmp2 > ${script0}Trial.py
mv ${script0}Trial.py /home/wud18/python

## Generating the corresponding .slurm file
# Read input parameters
read -p "input partation: " partation
read -p "input tasks_per_node: " tasks_per_node
read -p "input memory: " memory
read -p "input time(days): " days

# Replacing with input parameters
sed "s/PARTITION/$partation/" /home/wud18/python/Trial_maker/Trial.slurm > tmp6
sed "s/TASKS_PER_NODE/$tasks_per_node/" tmp6 > tmp7
sed "s/MEM/$memory/" tmp7 > tmp6
sed "s/DAYS/$days/" tmp6 > tmp7

# Extract input letters
cat $script | grep inputs.add_argument > tmp4
sed -i '1d' tmp4
sed "s/',/+/" tmp4 > tmp5
sed "s/inputs.add_argument('//" tmp5 > tmp4
paste -d\  tmp4 tmp4 > tmp5
sed 's/+\ -/\ "${/' tmp5 > tmp4
sed 's/+/}"/' tmp4 > tmp5
paste -s -d\  tmp5 > tmp4
sed "s/-s/python\ \$sc\ -s/" tmp4 > tmp5

# Combining all lines
cat tmp7 tmp5 > ${script0}Submit.slurm
mv ${script0}Submit.slurm /home/wud18/bash

# Removing temporary files
rm tmp0 tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7
