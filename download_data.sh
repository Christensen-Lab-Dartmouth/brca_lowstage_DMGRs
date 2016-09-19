# All data is publicly available and was retrieved from The Cancer Genome Atlas
# Downloaded March 2015

url='https://zenodo.org/record/153930/files/TCGA_BRCA_IDAT.tar.gz'
wget --directory-prefix 'I.Data_Processing/Data' $url

# md5sum 'I.Data_Processing/Data/TCGA_BRCA_IDAT.tar.gz'
# 8f45c465a1077e1f52f4d6b2a24c5986

base_fh='I.Data_Processing/Data/'
tar --gzip --extract --file $base_fh'TCGA_BRCA_IDAT.tar.gz' --directory=$base_fh
