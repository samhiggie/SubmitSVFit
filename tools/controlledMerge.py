# From truggles
# https://github.com/truggles/Z_to_TauTau_13TeV/blob/master/util/svFitMerger.py

import ROOT
import os, glob, subprocess


# Check if directory exists, make it if not
def checkDir( dirName ) :
    if not os.path.exists( dirName ) : os.makedirs( dirName )




def mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir ) :
    #files = glob.glob(originalDir+'/%s%s_*_%s.root' % (jobId, sample, channel) )
    print originalDir+'/%s/*.root' % (sample)
    files = glob.glob(originalDir+'/%s/*.root' % (sample) )
    checkDir( targetDir+'/%s' % (sample) )

    rep = 0
    runningSize = 0
    toMerge = []
    ints = []
    for file_ in files :

        # Merge to ~ 1000 events per file
        f = ROOT.TFile(file_,'r')
        t = f.Get(ttreePath)
	print file_,t.GetEntries()
        size = t.GetEntries()
        print size,"   ",file_
        runningSize += size
        print "running size: ",runningSize
        toMerge.append( file_ )
        if runningSize > 15000 :
            runningSize = 0
            mergeList = ["hadd", "-f", targetDir+"/%s/ntuple_%i_%s.root" % (sample, rep, channel)]
            for f in toMerge :
                mergeList.append( f )
            subprocess.call( mergeList )
            ints = []
            toMerge = []
            rep += 1
    mergeList = ["hadd", "-f", targetDir+"/%s/ntuple_%i_%s.root" % (sample, rep, channel)]
    for f in toMerge :
        mergeList.append( f )
    if len( mergeList ) > 3 : # greater than 3 means we actually have a file to merge (not empty)
        subprocess.call( mergeList )



if __name__ == '__main__' :


    samples=["Out_DataA","Out_DY2_v1","Out_DY_low","Out_ST_t_antitop","Out_TTToHadronic","Out_W3_v1","Out_W_v1","Out_ZZ_v1","Out_DataB","Out_DY2_v2","Out_EWKWminus","Out_ST_t_top","Out_TTToSemiLeptonic","Out_W3_v2","Out_WW_v1","Out_DataC","Out_DY3_v1","Out_EWKWplus","Out_ST_tW_antitop","Out_VBFHTT","Out_W4_v1","Out_WW_v2","Out_DataD","Out_DY3_v2","Out_EWKZLL","Out_ST_tW_top","Out_W1_v1","Out_W4_v2","Out_WZ_v1","Out_DY1_v1","Out_DY4_v1","Out_EWKZNuNu","Out_ttHTT","Out_W2_v1","Out_WminusHTT","Out_WZ_v2","Out_DY1_v2","Out_DY_v1","Out_GGHTT","Out_TTTo2L2Nu","Out_W2_v2","Out_WplusHTT","Out_ZHTT","Out_embeddedA","Out_embeddedB","Out_embeddedC","Out_embeddedD","Out_WG","Out_GGZHNNTT","Out_GGZHQQTT","Out_GGZHLLTT","Out_GGHWW","Out_GGZHWW","Out_VBFHWW","Out_WminusHWW","Out_WplusHWW","Out_ZHWW"] #mutau 2018

    #samples=["Out_DY1_v1","Out_DY_v2","Out_DataF","Out_EWKZLL_v1","Out_ST_t_antitop","Out_W2_v1","Out_WZ2L2Q","Out_ZZ2L2Q","Out_DY1_v2","Out_DYlow","Out_DataG","Out_EWKZLL_v2","Out_ST_t_top","Out_W3_v1","Out_WZJLLLNu","Out_ZZ4L","Out_DY2_v1","Out_DataBv1","Out_DataH","Out_EWKZNuNu_v1","Out_TT","Out_W4_v1","Out_W_v1","Out_ttHnonBB","Out_DY2_v2","Out_DataBv2","Out_EWKWminus_v1","Out_EWKZNuNu_v2","Out_VBFHTT","Out_W4_v2","Out_W_v2","Out_DY3_v1","Out_DataC","Out_EWKWminus_v2","Out_GGHTT","Out_VV2L2Nu_v1","Out_WW1L1Nu2Q","Out_WminusHTT","Out_DY4_v1","Out_DataD","Out_EWKWplus_v1","Out_ST_tW_antitop","Out_VV2L2Nu_v2","Out_WZ1L1Nu2Q","Out_WplusHTT","Out_DY_v1","Out_DataE","Out_EWKWplus_v2","Out_ST_tW_top","Out_W1_v1","Out_WZ1L3Nu","Out_ZHTT"] #mutau 2016

    #samples=["Out_DY1","Out_DY_v2","Out_DataE","Out_EWKZNuNu","Out_ST_t_top","Out_W1_v1","Out_W4","Out_WplusHTT","Out_DY2","Out_DYlow","Out_DataF","Out_GGHTT","Out_TTTo2L2Nu","Out_W1_v2","Out_WW","Out_ZHTT","Out_DY3","Out_DataB","Out_EWKWminus","Out_ST_tW_antitop","Out_TTToHadronic","Out_W2_v1","Out_WZ","Out_ZZ","Out_DY4","Out_DataC","Out_EWKWplus","Out_ST_tW_top","Out_TTToSemiLeptonic","Out_W2_v2","Out_W_v1","Out_DY_v1","Out_DataD","Out_EWKZLL","Out_ST_t_antitop","Out_VBFHTT","Out_W3","Out_WminusHTT"] #2017

    #samples=["Out_DY1","Out_DY_v2","Out_EWKZNuNu","Out_ST_t_top","Out_W1_v1","Out_W4","Out_WplusHTT","Out_DY2","Out_DYlow","Out_GGHTT","Out_TTTo2L2Nu","Out_W1_v2","Out_WW","Out_ZHTT","Out_DY3","Out_EWKWminus","Out_ST_tW_antitop","Out_TTToHadronic","Out_W2_v1","Out_WZ","Out_ZZ","Out_DY4","Out_EWKWplus","Out_ST_tW_top","Out_TTToSemiLeptonic","Out_W2_v2","Out_W_v1","Out_DY_v1","Out_EWKZLL","Out_ST_t_antitop","Out_VBFHTT","Out_W3","Out_WminusHTT"] #2017 wo data

    #samples=["Out_embeddedB","Out_embeddedC","Out_embeddedD","Out_embeddedE","Out_embeddedF"] # embedded 2017

    #samples=["Out_DY3_v2"]


    originalDir = '/nfs_scratch/caillol/smhem2018_4aug'
    targetDir = '/nfs_scratch/caillol/smhem2018_4aug_merged'
    jobId = ''
    channel = 'xx'
    ttreePath = 'emu_tree'
    for sample in samples :
	print sample
        mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir )


