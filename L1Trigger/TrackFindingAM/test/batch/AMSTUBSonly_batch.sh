#!/bin/bash

###########################################
#
# Main script for parallel CMSSW AM-based pattern 
# recognition on lxbatch
# 
# !!!Working on EDM files!!!
#
# The jobs themselves are launched by PR_processor_parallel.sh
#
# source AMSTUBSonly_batch.sh p1 p2 p3
# with:
# p1 : The directory containing the data file you want to analyze
# p2 : Name of the SE subdirectory where you will store the data
# p3 : The global tag name
#
# For more details, and examples, have a look at:
# 
# http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCTuto (STEP V.2)
#
# Author: S.Viret (viret@in2p3.fr)
# Date  : 01/03/2016
#
# Script tested with release CMSSW_6_2_0_SLHC27
#
###########################################


# Here we retrieve the main parameters for the job 

INPUTDIR=${1}   # Directory where the input root files are
OUTPUTDIR=${2}  # Name of the output directory in the SE
GTAG=${3}       # Global tag

###################################
#
# The list of parameters you can modify is below
#
###################################

export LFC_HOST=lyogrid06.in2p3.fr
export LFC_INPUT=xrootd.unl.edu

# I/O Directory info on the grid
BAREPATH=/dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/PR/Synchro/Input
SRMPATH=srm://$LFC_HOST/$BAREPATH
XRDPATH=root://$LFC_HOST/$BAREPATH

# The queue over which you want to send the job
BQUEUE=2nd

###########################################################
###########################################################
# You are not supposed to touch the rest of the script !!!!
###########################################################
###########################################################

# The SE directory containing the output EDM file with the PR output
OUTDIR_GRID=$SRMPATH/${INPUTDIR}_$OUTPUTDIR
OUTDIR_XROOT=$XRDPATH/${INPUTDIR}_$OUTPUTDIR

MATTER=`basename $OUTDIR_GRID`

OUTDIRTMP=TMPDAT_$MATTER

echo 'The data will be read from directory: '$INDIR_GRID
echo 'The final pattern reco output files will be written in: '$OUTDIR_GRID

# Following lines suppose that you have a certificate installed on lxplus. To do that follow the instructions given here:
#
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideLcgAccess
#

voms-proxy-init --voms cms --valid 200:00 -out $HOME/.globus/gridproxy.cert
export X509_USER_PROXY=${HOME}/.globus/gridproxy.cert

rm temp,files.txt

python /afs/cern.ch/work/s/sviret/testarea/SLHC_sample/das_client.py --query="dataset dataset=/RelValTTbar_14TeV/"$INPUTDIR"/GEN-SIM-DIGI-RAW release=CMSSW_6_2_0_SLHC26_patch2 | sort dataset.creation_time-" --limit=30 > temp


samplelist=`cat temp`
for f in $samplelist
do
python /afs/cern.ch/work/s/sviret/testarea/SLHC_sample/das_client.py --query="file dataset="$f --limit=3000 | grep store > files.txt
#sed "s'/store'root://xrootd.unl.edu//store'" -i files.txt

filelist=`cat files.txt`
done 

gfal-mkdir $OUTDIR_GRID

# Then we setup some environment variables

PACKDIR=$PWD            # This is where the scripts are created 
RELEASEDIR=$CMSSW_BASE  # This is where the release is installed

# We loop over the data directory in order to find all the files to process

nsj=0

for ll in $filelist 
do   
    l=`basename $ll`
    echo $l
    i=0
    j=$NPFILE

    nsj=$(( $nsj + 1))

    echo "#\!/bin/bash" > run_${nsj}_${MATTER}.sh
    echo "$PACKDIR/run_STUBS_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
    echo "$PACKDIR/run_RM_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
    
    echo "#\!/bin/bash" > run_STUBS_${nsj}_${MATTER}.sh
    echo "#\!/bin/bash" > run_RM_${nsj}_${MATTER}.sh
    
    chmod 755 run_${nsj}_${MATTER}.sh
    chmod 755 run_STUBS_${nsj}_${MATTER}.sh
    chmod 755 run_RM_${nsj}_${MATTER}.sh

    echo 'Working with file '$l

    # First look if the file has been processed

    OUTM=`echo $l | cut -d. -f1`

    OUTF=${OUTM}"_with_STUBS.root"

    AMPR_FILE=${OUTDIR_GRID}/$OUTF

    dealF=`gfal-ls $AMPR_FILE | wc -l`

    if [ $dealF != "0" ]; then
	
	echo "File "$l" has already processed, skip it..."
	continue;
    fi

    echo "$PACKDIR/PR_processor_parallel.sh STUBS $ll $AMPR_FILE $RELEASEDIR $GTAG $OUTDIRSTO" >> run_STUBS_${nsj}_${MATTER}.sh

    echo "cd $OUTDIRSTO" >> run_RM_${nsj}_${MATTER}.sh
    echo "rm *${OUTM}_*root" >> run_RM_${nsj}_${MATTER}.sh

    if [ ${4} == "BATCH" ]; then	
	bsub -R "pool>5000" -q $BQUEUE -e /dev/null -o /tmp/${LOGNAME}_out.txt run_${nsj}_${MATTER}.sh
    fi
   
done

