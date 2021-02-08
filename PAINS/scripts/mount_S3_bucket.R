
mount_S3_bucket <- function(){
    system("
sudo s3fs \\
     cellprofilerdata \\
     -o use_cache=/home/ubuntu/tmp \\
     -o uid=1001 \\
     -o mp_umask=002 \\
     -o multireq_max=5 \\
     -o iam_role=\"SextonS3\" \\
     -o allow_other \\
     ~/bucket_cellprofilerdata")
}
