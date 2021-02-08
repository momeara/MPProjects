


system("
s3fs \\
    sextoncov19 \\
    -o use_cache=/home/ubuntu/tmp \\
    -o uid=1001 \\
    -o mp_umask=002 \\
    -o multireq_max=5 \\
    -o iam_role=\"SextonS3\" \\
    -o allow_other ~/bucket")


system("
s3fs \\
    umich-insitro \\
    -o use_cache=/home/ubuntu/tmp \\
    -o uid=1001 \\
    -o mp_umask=002 \\
    -o multireq_max=5 \\
    -o iam_role=\"SextonS3\" \\
    -o allow_other ~/bucket_umich-insitro")

system("mkdir product/figures/TS2_infected")
