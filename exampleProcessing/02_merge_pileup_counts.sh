# This script takes all the A,C,G,T pileups for single cells and merges them
directory=$1
samplename=$2

# Combine all files
cat ${directory}/*.A.txt | gzip > "${directory}/${samplename}_all.A.txt.gz"
cat ${directory}/*.C.txt | gzip > "${directory}/${samplename}_all.C.txt.gz"
cat ${directory}/*.G.txt | gzip > "${directory}/${samplename}_all.G.txt.gz"
cat ${directory}/*.T.txt | gzip > "${directory}/${samplename}_all.T.txt.gz"
cat ${directory}/*.coverage.txt | gzip > "${directory}/${samplename}_all.coverage.txt.gz"
