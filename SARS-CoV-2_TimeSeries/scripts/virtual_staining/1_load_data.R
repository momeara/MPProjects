
# load a single well from the 48-hour infection condition
# DPC is digital phase contrast, where through-focus is used
# drive partially coherent phase constast.
system("
destination=raw_data/virtual_staining/20200531T160447_48_hour_DPC
mkdir -p $destination
for W in {0001..0384}; do
for F in {0001..0012}; do
for T in 0001; do
for Z in {001..030}; do
for C in 1; do
   aws s3 cp s3://umich-insitro/CQ1/20200531T160447_48_hour_DPC/Image/W${W}F${F}T${T}Z${Z}C${C}.tif ${destination}/
done
done
done
done
done
")

# the stains for the same well
system("
destination=raw_data/virtual_staining/20200530T160342_48_hour
mkdir -p $destination
for W in {0001..0384}; do
for F in {0001..0012}; do
for T in 0001; do
for Z in {001..010}; do
for C in 1; do
   aws s3 cp s3://umich-insitro/CQ1/20200530T160342_48_hour/Image/W${W}F${F}T${T}Z${Z}C${C}.tif ${destination}/
done
done
done
done
done
")
