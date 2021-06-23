###### PLEASE CONFIGURE THESE PARAMETERS #####
SETUPFILE=/home/pyadmin/software/setup_scripts/setup-pythia8.sh #contains env. variables
COMPILER=compile_pythia.sh
export CASEFILE=case.sh
source $CASEFILE
NJOBS=50 #WARNING: BE AWARE THAT THE FILES PRODUCED BY EACH JOB WILL HAVE 
	 #1/NJOBS AS NORMALIZATION, TO PRODUCE MERGED FILES PROPERLY NORMALIZED.
echo "----------------------------------"
echo "----------------------------------"
echo "----------------------------------"
echo $CASE
echo "----------------------------------"
echo "----------------------------------"
echo "----------------------------------"

###### 
rm *.root *.exe
source $SETUPFILE
./$COMPILER splittinggluon

rm -rf $OUTPUTFOLDER
mkdir $OUTPUTFOLDER
cd $OUTPUTFOLDER
for i in $( eval echo {1..$NJOBS} )
do
   mkdir file_$i
   cd file_$i
   cp ../../splittinggluon.exe .
   cp ../../config_splitting.yaml .
   echo $RANDOM
   ./splittinggluon.exe $CASE $RANDOM $NJOBS &
   cd ..
done

