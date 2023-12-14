#!/usr/bin/bash

echo hi

input_dir=$1

echo input dir $input_dir

if [ -d "$input_dir" ];
then

   echo "found input directory $input_dir"

   if [ -d $input_dir/fit-results ];
   then
      echo "fit-results subdirectory already exists"
   else
      echo "making fit-results subdirectory"
      mkdir -p $input_dir/fit-results
   fi

   for wf in $input_dir/boot*.npy
   do
      fileonly=${wf##*/}
      fitlog="${fileonly/.npy/.fitlog}"
      fitlogwithpath="${input_dir}/fit-results/${fitlog}"
      echo fileonly $fileonly    fitlog $fitlog    fitlogwithpath $fitlogwithpath
      echo
      echo command: python3 run-ndfit-1a.py  ${input_dir} train-and-true-samples.npy ${fileonly} ${fitlogwithpath}
      python3 run-ndfit-1a.py  ${input_dir} train-and-true-samples.npy ${fileonly} |& tee ${fitlogwithpath}
   done

else
   echo "directory $input_dir does not exist."
   exit -1
fi

