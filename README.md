#FStitch
##Usage and Output
Scientists looking to classify regions of high read density in Globab Run On sequencing will find Fast Read Stitcher (FStitch) most useful as it 
identifies putative nascent transcripts _de novo_<sup>1</sup>. However, users may also find this package useful as a ChIP-seq peak caller.


![Alt text](https://github.com/azofeifa/FStitch/blob/master/images/IGV_SNAP.png)


The above IGV Snap shot displays the classifications given by FStitch. Color ‘green’ is read data considered inactive noise. Color ‘blue’ is a putative nascent transcript on the forward strand and ‘red’ on the reverse strand. We identify both genes undergoing nascent transcription and many regions that are unannotated characteristic of enhancer elements. 


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
The fast read stitcher program attempts to classify and identify contiguous regions of read coverage that are showing strong signal over background mapping noise. With this in mind, FStitch requires a BedGraph file. Where for each genomic position, the number of reads mapping to that position are provided. This commonly known as a BedGraph file<sup>2</sup>. Briefly a BedGraph file consists of four columns: chromosome, start genomic coordinate, stop genomic coordinate, coverage. Below is an example:

  
![Alt text](https://github.com/azofeifa/FStitch/blob/master/images/BedGraphScreenShot.png)

Note: FStitch does not accept bed graph files where 0 coverage values are reported. In short, you can convert your bam files to a bed graph file format using the _bedtools_<sup>3</sup> command:

$bedtools genomecov -ibam <bamfile> -g <genome_file> -bg

We note that specifying five prime (-5) in the “genomecov” may allow for cleaner annotations however unspecified five prime bed works just fine as well. 
 


##FStitch train
FStitch uses two probabilistic models to classify regions of high read density that may be indicative of nascent transcription (GRO-seq) or a read coverage peak (ChIP-seq): Logistic Regression and a Hidden Markov Model. The logistic regression coefficients are estimated via a user defined label training file.  Sense we are classifying regions as signal or noise, FStitch requires regions of the genome that show characteristic transcription or high read dense profiles and regions of the genome that display noise or not a profile of nascent transcription or a read dense region. With this information, FStitch trains a logistic regression classifier and then couples it to a hidden markov model. The transition parameters for the HMM are learned via the Baum Welch algorithm and thus do not require user label training data.  

In short, FStitch requires regions the user considers active transcription (or a peak) and regions considered inactive (simply noise). We note that the more regions provided to FStitch the more accurate the classifications however we have in Cross Validation<sup>1</sup> analysis that roughly 5-10 regions of active and inactive regions will yield highly accurate classifications. These regions are provided to FStitch using a specific file format with four columns separated by tabs: chromosome, genomic coordinate start, genomic coordinate stop, (0 if “noise” or 1 “inactive”). An example is given below:

![Alt text](https://github.com/azofeifa/FStitch/blob/master/images/TrainingFileImage.png)
    


###parameters
##FStitch segment
###parameters
## Understanding and Interpreting Output
##References
1. 

2.


