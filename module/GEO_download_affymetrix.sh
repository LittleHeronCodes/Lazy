# GSE_ID=
# DIR=

echo $GSE_ID
mkdir -p $DIR/$GSE_ID

URL='https://www.ncbi.nlm.nih.gov/geo/download/?acc='${GSE_ID}'&format=file'
wget $URL -O $DIR/$GSE_ID/${GSE_ID}_RAW.tar

tar -xf $DIR/$GSE_ID/${GSE_ID}_RAW.tar -C $DIR/$GSE_ID/
rename -v 's/\.cel./\.CEL./' *.cel.*
gunzip -f $GSE_ID/*.CEL.gz


