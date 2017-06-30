# FStitch
Below is the README for compiling, running and understanding Fast Read Stitcher (FStitch). This software is written primarily for scientists looking to identify putative nascent transcripts _de novo_ in Global Run-On sequencing<sup>1</sup>. However, users may also find this package useful as a ChIP-seq peak caller.

![Alt text](https://github.com/azofeifa/FStitch/blob/master/images/IGV_SNAP.png)


Above: IGV snap shot displays the classifications given by FStitch. Color ‘green’ is read data considered inactive noise. Color ‘blue’ is a putative nascent transcript on the forward strand and ‘red’ on the reverse strand. We identify both genes undergoing nascent transcription and many regions that are unannotated, characteristic of enhancer elements. 

## Brief Overview of FStitch Commands
Here are the minimal commands needed to run FStitch from start to finish; for greater detail on usage and file types see below. 
```
$ FStitch train -i </path/to/BedGraphFile> -j </path/to/TrainingFile>  -o </path/to/Parameters.out>

$ FStitch segment -i </path/to/forward/BedGraphFile> -j </path/to/reverse/BedGraphFile> -k </path/to/Parameters.out> -o </path/to/Classifications.bed>
```

## System Requirements
FStitch is written in the C++ programming language, with C++11 support and uses OpenMP<sup>4</sup> to parallelize portions of the program.  With this in mind, users will need to have a GCC compilers later than version 4.7 to compile and run FStitch. For mac users, downloading the latest Xcode will update the GCC compiler need be. To check you compiler version, 
```
$ gcc —-version
$ g++ —-version
```
Note, for those running FStitch on a compute cluster, commonly you will need to perform a ‘module load gcc<version>’ to compile FStitch. Please ask your sys admins for questions on module load behavior. 
##Setup
Download the FastReadStitcher/ directory from this url or clone to your local machine. If your compiler is up to date, you can compile FStitch by moving into the FastReadStitcher/ directory and running 
```
$ sh setup.sh
=========================================
Sucessfully Compiled
```
In short, the setup.sh just runs “make” in the src/ directory. If everything compiles, you should see "Sucessfully Compiled" at the end.

Importantly, you will now see the executable “FStitch” in the src directory. This will be the first command used for all the following computations. 
##Input File
The fast read stitcher program attempts to classify and identify contiguous regions of read coverage that are showing strong signal over background mapping noise. With this in mind, FStitch requires a file where for each genomic position, the number of reads mapping to that position are provided. This file is commonly known as a BedGraph file<sup>2</sup>. Briefly a BedGraph file consists of four columns: chromosome, start genomic coordinate, stop genomic coordinate, coverage. Below is an example:
  
![Alt text](https://github.com/azofeifa/FStitch/blob/master/images/BedGraphScreenShot.png)

Note: FStitch does not accept bedgraph files where 0 coverage values are reported and if the data is stranded (like GRO-seq) then there should be one BedGraph file corresponding to the positive strand and the negative strand . In short, you can convert your bam files to a bed graph file format using the _bedtools_<sup>3</sup> command:
```
$ bedtools genomecov -ibam <bamfile> -g <genome_file> -bg -s “+/-“
```
We note that specifying five prime (-5) in the “genomecov” may allow for cleaner annotations however unspecified five prime bed works just fine as well. 

##Running FStitch
The Fast Read Stitcher program is divided into two main commands: “train” and “segment”. “train” estimates the necessary probabilistic model parameters and “segment” pulls the output from “train” and classifies the entire genome into _active_ and _inactive_ regions of regions of high density read coverage. 

Note that a quick reference to the below parameters and software usage can be supplied by 
```
$/src/FStitch -h 

$/src/FStitch --help
```


## FStitch train
FStitch uses two probabilistic models to classify regions of high read density that may be indicative of nascent transcription (GRO-seq) or a read coverage peak (ChIP-seq): Logistic Regression and a Hidden Markov Model. The logistic regression coefficients are estimated via a user defined label training file.  Sense we are classifying regions as signal or noise, FStitch requires regions of the genome that show characteristic transcription or high read dense profiles and regions of the genome that display noise or not a profile of nascent transcription or a read dense region. With this information, FStitch trains a logistic regression classifier and then couples it to a Markov model. The transition parameters for the Markov model are learned via the Baum Welch algorithm and thus do not require user label training data.  

In short, FStitch requires regions the user considers active transcription (or a peak) and regions considered inactive (simply noise). We note that the more regions provided to FStitch the more accurate the classifications however we have in Cross Validation<sup>1</sup> analysis that roughly 10-15 regions of active and inactive regions will yield accurate classifications. These regions are provided to FStitch using a specific file format with four columns separated by tabs: chromosome, genomic coordinate start, genomic coordinate stop, (0 if “noise” or 1 “signal”). An example is given below:

![Alt text](https://github.com/azofeifa/FStitch/blob/master/images/TrainingFileImage2.png)

The segments do not need to be in any order and can be from any chromosome, however each region must not overlap any other segment as this will cause confusion in the learning algorithms for the logistic regression classifier. 

Running FStitch train is simple once you have your data in the correct format and have created the training file above. A description of the parameters for FStitch train are given below

|Flag|Type|Desription|
|----|----|----------|
|-i	 | \</path/to/BedGraphFile> | BedGraph File from above
| -j | \</path/to/TrainingFile> | Training File from above
| -o | \</path/to/anyName> | TrainingParameterOutFile
| -np| number | number of processors, default 8
| -al| number |learning rate for newtons method, default 1
| -cm| number | max number of iterations for Baum-Welch, default 100
| -ct| number | convergence threshold for Baum Welch, default 0.01
| -reg| number | regularization parameter for logistic regression classifier, default 1
Putting this together
```
$ /src/FStitch train -i \</path/to/BedGraphFile\> -j \</path/to/TrainingFile> -o \</path/to/anyName.out>
```
This will produce the a fie called anyName.out that will store the learned parameters for the logistic regression and HMM transition parameters need in “FStitch segment”. Below is one such output for anyName.out

![Alt text](https://github.com/azofeifa/FStitch/blob/master/images/ParameterOutFile.png)

Very important: If FStitch is being used on stranded data, the BedGraph file used in the “FStitch train” command must correspond to the strand indicated in the TrainingFile. For example, if the strand in the training file comes from the forward strand but the user supplies a BedGraph file that is on the reverse strand, then learned parameters will not be accurate. 


## FStitch segment
FStitch segment follows from FStitch train and takes as input the TrainingParameterOutFile (from above, \</path/to/anyName.out>) as input, along with the original BedGraph file. A description of the parameters for FStitch segment are given below

|Flag|Type|Desription|
|----|----|----------|
| -i	| \</path/to/BedGraphFile> |BedGraph File Format from above
| -k 	| \</path/to/anyName.out> |Training Parameter Out File from FStitch train call
| -o	| \</path/to/anyName.bed> |A bed file that gives the regions considered active nascent transcription (or ChIP-seq peak) and noise
| -np 	| number |number of processors, default 8

Putting this together
```
$ /src/FStitch segment -i \</path/to/BedGraphFile\> -k \</path/to/anyName.out> -o \</path/to/anyName.bed> 
```
This will produce a file called anyName.bed, and can be imported into any genome browser)

Note you can use your parameter out file from FStitch train (i.e. anyName.out) to segment other datasets. In fact, using the same parameter out file will gurantee consistency and comparibility across datasets, so this is encouraged.

## Cite
If you find the Fast Read Stitcher program useful for your research please cite:

Joseph Azofeifa, Mary A. Allen, Manuel Lladser, and Robin Dowell. 2014. __FStitch: a fast and simple algorithm for detecting nascent RNA transcripts.__ In Proceedings of the 5th ACM Conference on Bioinformatics, Computational Biology, and Health Informatics (BCB '14). ACM, New York, NY, USA, 174-183. 

##References 
1. Joseph Azofeifa, Mary A. Allen, Manuel Lladser, and Robin Dowell. 2014. __FStitch: a fast and simple algorithm for detecting nascent RNA transcripts.__ In Proceedings of the 5th ACM Conference on Bioinformatics, Computational Biology, and Health Informatics (BCB '14). ACM, New York, NY, USA, 174-183. DOI=10.1145/2649387.2649427 http://doi.acm.org/10.1145/2649387.2649427   
2. http://genome.ucsc.edu/goldenpath/help/bedgraph.html
3. http://bedtools.readthedocs.org/en/latest/
4. http://openmp.org/wp/

