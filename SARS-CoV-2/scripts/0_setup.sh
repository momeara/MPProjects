#!/bin/bash


############################
# Setup access to raw data #
############################

if [ ! -d ~/bucket ]
then
    mkdir ~/bucket
fi
if [ $(mountpoint bucket) == "bucket is a mountpoint" ]
then
    s3fs \
        sextoncov19 \
        -o use_cache=/home/ubuntu/tmp \
        -o uid=1001 \
        -o mp_umask=002 \
        -o multireq_max=5 \
        -o iam_role="SextonS3" \
        -o allow_other \
        ~/bucket
fi

if [ ! -d ~/bucket ]
then
    mkdir ~/bucket-insitro
fi
if [ $(mountpoint bucket-insitro) == "bucket-instro is a mountpoint" ]
then
    s3fs \
        umich-insitro \
        -o use_cache=/home/ubuntu/tmp \
        -o uid=1001 \
        -o mp_umask=002 \
        -o multireq_max=5 \
        -o iam_role="SextonS3" \
        -o allow_other ~/bucket_umich-insitro
fi


###################################################
# sync study metadata into the raw_data directory #
###################################################

cp ~/bucket/screen_metadata/* raw_data/
