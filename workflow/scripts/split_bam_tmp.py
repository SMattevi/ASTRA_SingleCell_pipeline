import argparse
import os
import fileinput
import pysam

commandLineArgumentParser = argparse.ArgumentParser()
commandLineArgumentParser.add_argument("-cl", "--cluster_file",  help="cluster file")
commandLineArgumentParser.add_argument("-b","--bam", help="bam file")
commandLineArgumentParser.add_argument("-o","--out", help="output directory")

commandLineArguments = commandLineArgumentParser.parse_args()

clfile = commandLineArguments.cluster_file
bam = commandLineArguments.bam
outdir = commandLineArguments.out

if not os.path.exists(outdir):
    os.makedirs(outdir)

try:
    with open(clfile) as CLUSTER:
        myclcell = {}
        myclfiles = {}
        myclfilehandles = {}

        myheader = CLUSTER
        print("I am reading cluster information file: ", clfile)

        mypaircounter = 0 
        
        for line in myheader:
            line = line.replace('\n', '') 
            mypaircounter=mypaircounter+1
            
            myarray = line.split("\t")
            mybarcode= myarray[0]
            mycluster= myarray[1]
            myclcell[mybarcode] = mycluster
            myclfile = outdir+"/"+"cluster_"+mycluster+".sam"
            myclfiles[mycluster] = myclfile
        print("I have successfully read",mypaircounter,"cell-cluster pairs from", clfile)

        print("I am writing head to sam files")

        cl_counter = 0
        for mycluster in myclfiles:
            print(mycluster)
            cl_counter=cl_counter+1
            myclfile=myclfiles[mycluster]
            cluster_file_handles = {}
            cluster_file_handles[mycluster] = open(myclfile, 'w')
            try:
                with open("samtools view -H "+bam) as BAMHEAD:
                    myhh=BAMHEAD
                    cluster_file_handles[mycluster].write(myhh)
                BAMHEAD.close()
            except IOError:
                print("Cannot read file ", "samtools view -H ",bam)
        print("I have found",cl_counter,"clusters in your data and initialized one sam file for each of them.")
except IOError:
    print("Cannot read file ", clfile)



