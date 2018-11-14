# SubmitSVFit
```

cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
git clone https://github.com/SVfit/ClassicSVfit TauAnalysis/ClassicSVfit -b fastMTT_21_06_2018
git clone https://github.com/SVfit/SVfitTF TauAnalysis/SVfitTF
cd $CMSSW_BASE/src
scram b -j 4
```

There are a number of current versions of the stand alone svFitter included here.
You can find them in ROOT/bin/SVFitStandAlone... with the names of their executables
defined ROOT/bin/BuildFile.xml.

To run in interactive mode for example:
```
SVFitStandAloneFSA inputFile=coolInputFile.root newOutputFile=1 newFile=tmpOut.root doES=1
```

 - inputFile = obvious
 - newOutputFile = 0/1
   - 0 = update input file with svFit vars
   - 1 = output new file with original TTree and new svFit vars
 - newFile = name of output file, default is newFile.root if none specified
 - doES = apply energy scale adjustments providing nominal, shift UP and shift DOWN values of svFit
   - 0 = default, no shift
   - 1 = apply shifts

To submit jobs to condor:
```
cd test
python svFitSubmitter.py -dr -sd /hdfs/store/user/truggles/svFitTestSept03 -es=1 --jobName svFitSept03Test
```

 - -dr = dryRun and outputs a command for you to run
 - -sd = select directory, the input directory, this will produce a list of all files in that directory to run over<BR>
       you must have your files in a /hdfs directory.
 - -es = apply energy scale, see above
 - --jobName = applys a job name to append to your file names and output directory structure

To get your files from elsewhere to /hdfs do something like this:
```
mkdir /hdfs/store/user/truggles/mySubmitDir
rsync -ahP /nfs_scratch/truggles/httSept04skimMerged/*.root /hdfs/store/user/truggles/httSept04skimMerge/
```

