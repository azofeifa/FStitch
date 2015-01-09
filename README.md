#FStitch
##Usage
Scientists looking to classify regions of high read density in Globab Run On sequencing will find Fast Read Stitcher (FStitch) most useful as it identifies putative nascent transcripts _de novo_<sup>1</sup>. However, users may also find this package useful as a ChIP-seq peak caller.
##Output

![alt text](https://github.com/azofeifa/FStitch/master/images/IGV_SNAP.png “FStitch Output in IGV“)

##System Requirements
FStitch is written in the C++ programming language and uses OpenMP to parallelize portions of the program.  With this in mind, users will need to have a GCC compilers later than version 4.2 to compile and run FStitch. For mac users, downloading the latest Xcode will update the GCC compiler need be. To check you compiler version, 

$gcc —-version

or 

$g++ —-version

Note, for those running FStitch on a compute cluster, commonly you will need to perform a ‘module load gcc<version>’ to compile FStitch. Please ask your sys admins for questions on module load behavior. 
##Setup
If your compiler is up to date, you can compile FStitch by moving into the FastReadStitcher/ directory and running 

$sh setup.sh

This runs “make” in the src/ directory. If everything compiles, you should see at the end of the compilation:

$=========================================

$Sucessfully Compiled

Importantly, will you now see the executable “FStitch” in the src directory. This will be the command used for the following computations. 
##Bedgraph Files
The fast read stitcher attempts to classify and identify contiguous regions of read coverage that are showing strong signal over background mapping noise. With this in mind, FStitch requires a BedGraph file. Where for each genomic position, the number of reads mapping to that position are provided. This commonly known as a BedGraph file<sup>2</sup>.   





##FStitch train
FStitch uses two probabilistic models to classify regions of high read density that may be indicative of nascent transcription (GRO-seq) or a read coverage peak (ChIP-seq): Logistic Regression and a Hidden Markov Model. The logistic regression coefficients are estimated via a user defined label training file.  Sense we are classifying regions as signal or noise, FStitch requires regions of the genome that show characteristic transcription or high read dense profiles and regions of the genome that display noise or not a profile of nascent transcription or a read dense region. With this information, FStitch trains a logistic regression classifier and then couples it to a hidden markov model. The transition parameters for the HMM are learned via the Baum Welch algorithm and thus do not require user label training data.  


This file comprise four columns that are separated by tabs.  



###parameters
##FStitch segment
###parameters
## Understanding and Interpreting Output
##References
1. 

2.


