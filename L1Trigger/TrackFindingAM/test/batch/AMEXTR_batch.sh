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
# source AMPR_batch.sh p1 p2 p3
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

# I/O Directory info on the grid
BAREPATH=/dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/PR/ReSynchro
SRMPATH=srm://$LFC_HOST$BAREPATH
XRDPATH=root://$LFC_HOST/$BAREPATH

# The queue over which you want to send the job
BQUEUE=1nd

###########################################################
###########################################################
# You are not supposed to touch the rest of the script !!!!
###########################################################
###########################################################


# The SE directory containing the input EDM files

INDIR_GRID=$SRMPATH/$INPUTDIR
INDIR_XROOT=$XRDPATH/$INPUTDIR

# The SE directory containing the output EDM file with the PR output
OUTDIR_GRID=$SRMPATH/$OUTPUTDIR
OUTDIR_XROOT=$XRDPATH/$OUTPUTDIR

MATTER=`basename $OUTPUTDIR`

OUTDIRTMP=TMPDAT_$MATTER
OUTDIRSTO=TMP_$MATTER

echo 'The data will be read from directory: '$INDIR_GRID
echo 'The final pattern reco output files will be written in: '$OUTDIR_GRID

# Following lines suppose that you have a certificate installed on lxplus. To do that follow the instructions given here:
#
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideLcgAccess
#

voms-proxy-init --voms cms --valid 100:00 -out $HOME/.globus/gridproxy.cert
export X509_USER_PROXY=${HOME}/.globus/gridproxy.cert

gfal-mkdir $SRMPATH/$OUTPUTDIR

# Then we setup some environment variables

PACKDIR=$PWD            # This is where the scripts are created 
RELEASEDIR=$CMSSW_BASE  # This is where the release is installed

# We loop over the data directory in order to find all the files to process

nsj=0

for ll in `gfal-ls $INDIR_GRID | grep AMPR` 
do   
    l=`basename $ll`
    echo $l
    i=0
    j=$NPFILE

    nsj=$(( $nsj + 1))

    echo "#\!/bin/bash" > run_${nsj}_${MATTER}.sh
    echo "$PACKDIR/run_EXTR_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
    echo "$PACKDIR/run_RM_${nsj}_${MATTER}.sh" >> run_${nsj}_${MATTER}.sh
    
    echo "#\!/bin/bash" > run_EXTR_${nsj}_${MATTER}.sh
    echo "#\!/bin/bash" > run_RM_${nsj}_${MATTER}.sh
    
    chmod 755 run_${nsj}_${MATTER}.sh
    chmod 755 run_EXTR_${nsj}_${MATTER}.sh
    chmod 755 run_RM_${nsj}_${MATTER}.sh

    echo 'Working with file '$l

    # First look if the file has been processed

    OUTM=`echo $l | cut -d. -f1`

    OUTE=${OUTM}"_extr.root"

    EXTR_FILE=${OUTDIR_GRID}/$OUTE

    dealI=`gfal-ls $EXTR_FILE | wc -l`

    if [ $dealI != "0" ]; then	
	echo "File "$l" has already processed, skip it..."
	continue;
    fi
	
    #
    # Last step, all the merged files for the given input
    # file have been processed. Then launch the final merging 
    # 

    echo "$PACKDIR/PR_processor_parallel.sh  EXTR $INDIR_XROOT/$l $OUTE -1 $OUTDIR_GRID $RELEASEDIR $GTAG $OUTDIRSTO $OUTB" >> run_EXTR_${nsj}_${MATTER}.sh
    echo "cd $OUTDIRTMP" >> run_RM_${nsj}_${MATTER}.sh
    echo "rm *${OUTM}_*root" >> run_RM_${nsj}_${MATTER}.sh

    if [ ${4} == "BATCH" ]; then	
	bsub -R "pool>7000" -q $BQUEUE -e /dev/null -o /tmp/${LOGNAME}_out.txt run_${nsj}_${MATTER}.sh
	#bsub -R "pool>20000" -q $BQUEUE run_${nsj}_${MATTER}.sh
    fi
   
done

