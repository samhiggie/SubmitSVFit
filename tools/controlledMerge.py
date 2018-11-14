# From truggles
# https://github.com/truggles/Z_to_TauTau_13TeV/blob/master/util/svFitMerger.py

import ROOT
import os, glob, subprocess


# Check if directory exists, make it if not
def checkDir( dirName ) :
    if not os.path.exists( dirName ) : os.makedirs( dirName )




def mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir ) :
    #files = glob.glob(originalDir+'/%s%s_*_%s.root' % (jobId, sample, channel) )
    print originalDir+'/Out_%s/*.root' % (sample)
    files = glob.glob(originalDir+'/Out_%s/*.root' % (sample) )
    checkDir( targetDir+'/Out_%s' % (sample) )

    rep = 0
    runningSize = 0
    runningNumFiles = 0
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
        runningNumFiles += 1
        print "running size: ",runningSize
        toMerge.append( file_ )
        if runningSize > 500 or runningNumFiles == 500 :
            runningSize = 0
            runningNumFiles = 0
            mergeList = ["hadd", "-f", targetDir+"/Out_%s/ntuple_%i_%s.root" % (sample, rep, channel)]
            for f in toMerge :
                mergeList.append( f )
            subprocess.call( mergeList )
            ints = []
            toMerge = []
            rep += 1
    mergeList = ["hadd", "-f", targetDir+"/Out_%s/ntuple_%i_%s.root" % (sample, rep, channel)]
    for f in toMerge :
        mergeList.append( f )
    if len( mergeList ) > 3 : # greater than 3 means we actually have a file to merge (not empty)
        subprocess.call( mergeList )



if __name__ == '__main__' :
    #samples=["ZZ4L_eeem","ZH120_eeem","ZH125_eeem","ZH130_eeem","GGZZ2e2mu_eeem","GGZZ2e2tau_eeem","GGZZ2mu2tau_eeem","GGZZ4tau_eeem","GGZZ4e_eeem","GGZZ4mu_eeem","DoubleEle_eeem","DoubleMu_eeem","ZZ4L_eeet","ZH120_eeet","ZH125_eeet","ZH130_eeet","GGZZ2e2mu_eeet","GGZZ2e2tau_eeet","GGZZ2mu2tau_eeet","GGZZ4tau_eeet","GGZZ4e_eeet","GGZZ4mu_eeet","DoubleEle_eeet","DoubleMu_eeet","ZZ4L_eemt","ZH120_eemt","ZH125_eemt","ZH130_eemt","GGZZ2e2mu_eemt","GGZZ2e2tau_eemt","GGZZ2mu2tau_eemt","GGZZ4tau_eemt","GGZZ4e_eemt","GGZZ4mu_eemt","DoubleEle_eemt","DoubleMu_eemt","ZZ4L_eett","ZH120_eett","ZH125_eett","ZH130_eett","GGZZ2e2mu_eett","GGZZ2e2tau_eett","GGZZ2mu2tau_eett","GGZZ4tau_eett","GGZZ4e_eett","GGZZ4mu_eett","DoubleEle_eett","DoubleMu_eett","ZZ4L_mmem","ZH120_mmem","ZH125_mmem","ZH130_mmem","GGZZ2e2mu_mmem","GGZZ2e2tau_mmem","GGZZ2mu2tau_mmem","GGZZ4tau_mmem","GGZZ4e_mmem","GGZZ4mu_mmem","DoubleEle_mmem","DoubleMu_mmem","ZZ4L_mmet","ZH120_mmet","ZH125_mmet","ZH130_mmet","GGZZ2e2mu_mmet","GGZZ2e2tau_mmet","GGZZ2mu2tau_mmet","GGZZ4tau_mmet","GGZZ4e_mmet","GGZZ4mu_mmet","DoubleEle_mmet","DoubleMu_mmet","ZZ4L_mmmt","ZH120_mmmt","ZH125_mmmt","ZH130_mmmt","GGZZ2e2mu_mmmt","GGZZ2e2tau_mmmt","GGZZ2mu2tau_mmmt","GGZZ4tau_mmmt","GGZZ4e_mmmt","GGZZ4mu_mmmt","DoubleEle_mmmt","DoubleMu_mmmt","ZZ4L_mmtt","ZH120_mmtt","ZH125_mmtt","ZH130_mmtt","GGZZ2e2mu_mmtt","GGZZ2e2tau_mmtt","GGZZ2mu2tau_mmtt","GGZZ4tau_mmtt","GGZZ4e_mmtt","GGZZ4mu_mmtt","DoubleEle_mmtt","DoubleMu_mmtt","WWZ_eeem","WWZ_eeet","WWZ_eemt","WWZ_eett","WWZ_mmem","WWZ_mmet","WWZ_mmmt","WWZ_mmtt","WZZ_eeem","WZZ_eeet","WZZ_eemt","WZZ_eett","WZZ_mmem","WZZ_mmet","WZZ_mmmt","WZZ_mmtt","ZZZ_eeem","ZZZ_eeet","ZZZ_eemt","ZZZ_eett","ZZZ_mmem","ZZZ_mmet","ZZZ_mmmt","ZZZ_mmtt"]
    #samples=["ZZ4L_eett","ZH120_eett","ZH125_eett","ZH130_eett","GGZZ2e2mu_eett","GGZZ2e2tau_eett","GGZZ2mu2tau_eett","GGZZ4tau_eett","GGZZ4e_eett","GGZZ4mu_eett","DoubleEle_eett"]
    #samples=["ZH110_eeem","ZH110_eeet","ZH110_eemt","ZH110_mmem","ZH110_mmet","ZH110_mmmt","ZH110_mmtt","ZH140_eeem","ZH140_eeet","ZH140_eemt","ZH140_mmem","ZH140_mmet","ZH140_mmmt","ZH140_mmtt","ZHWW125_eeem","ZHWW125110_eeet","ZHWW125_eemt","ZHWW125_mmem","ZHWW125_mmet","ZHWW125_mmmt","ZHWW125_mmtt"]
    #samples=["DoubleEle_eett"]#ZH110_eett","ZH140_eett","ZHWW125_eett"]

    #samples=["ZZ4L_eeem","ZZ4L_v2_eeem","ZH110_eeem","ZH120_eeem","ZH125_eeem","ZH130_eeem","ZH140_eeem","ZHWW125_eeem","DoubleEle_eeem","DoubleMu_eeem","SingleEle_eeem","SingleMu_eeem","WWZ_eeem","WZZ_eeem","ZZZ_eeem","ttZJets_eeem","GGZZ2e2mu_eeem","GGZZ2e2tau_eeem","GGZZ2mu2tau_eeem","GGZZ4tau_eeem","GGZZ4e_eeem","GGZZ4mu_eeem","ZZ4L_eeet","ZZ4L_v2_eeet","ZH110_eeet","ZH120_eeet","ZH125_eeet","ZH130_eeet","ZH140_eeet","ZHWW125_eeet","DoubleEle_eeet","DoubleMu_eeet","SingleEle_eeet","SingleMu_eeet","WWZ_eeet","WZZ_eeet","ZZZ_eeet","ttZJets_eeet","GGZZ2e2mu_eeet","GGZZ2e2tau_eeet","GGZZ2mu2tau_eeet","GGZZ4tau_eeet","GGZZ4e_eeet","GGZZ4mu_eeet","ZZ4L_eemt","ZZ4L_v2_eemt","ZH110_eemt","ZH120_eemt","ZH125_eemt","ZH130_eemt","ZH140_eemt","ZHWW125_eemt","DoubleEle_eemt","DoubleMu_eemt","SingleEle_eemt","SingleMu_eemt","WWZ_eemt","WZZ_eemt","ZZZ_eemt","ttZJets_eemt","GGZZ2e2mu_eemt","GGZZ2e2tau_eemt","GGZZ2mu2tau_eemt","GGZZ4tau_eemt","GGZZ4e_eemt","GGZZ4mu_eemt","ZZ4L_eett","ZZ4L_v2_eett","ZH110_eett","ZH120_eett","ZH125_eett","ZH130_eett","ZH140_eett","ZHWW125_eett","DoubleEle_eett","DoubleMu_eett","SingleEle_eett","SingleMu_eett","WWZ_eett","WZZ_eett","ZZZ_eett","ttZJets_eett","GGZZ2e2mu_eett","GGZZ2e2tau_eett","GGZZ2mu2tau_eett","GGZZ4tau_eett","GGZZ4e_eett","GGZZ4mu_eett","ZZ4L_mmem","ZZ4L_v2_mmem","ZH110_mmem","ZH120_mmem","ZH125_mmem","ZH130_mmem","ZH140_mmem","ZHWW125_mmem","DoubleEle_mmem","DoubleMu_mmem","SingleEle_mmem","SingleMu_mmem","WWZ_mmem","WZZ_mmem","ZZZ_mmem","ttZJets_mmem","GGZZ2e2mu_mmem","GGZZ2e2tau_mmem","GGZZ2mu2tau_mmem","GGZZ4tau_mmem","GGZZ4e_mmem","GGZZ4mu_mmem","ZZ4L_mmet","ZZ4L_v2_mmet","ZH110_mmet","ZH120_mmet","ZH125_mmet","ZH130_mmet","ZH140_mmet","ZHWW125_mmet","DoubleEle_mmet","DoubleMu_mmet","SingleEle_mmet","SingleMu_mmet","WWZ_mmet","WZZ_mmet","ZZZ_mmet","ttZJets_mmet","GGZZ2e2mu_mmet","GGZZ2e2tau_mmet","GGZZ2mu2tau_mmet","GGZZ4tau_mmet","GGZZ4e_mmet","GGZZ4mu_mmet","ZZ4L_mmmt","ZZ4L_v2_mmmt","ZH110_mmmt","ZH120_mmmt","ZH125_mmmt","ZH130_mmmt","ZH140_mmmt","ZHWW125_mmmt","DoubleEle_mmmt","DoubleMu_mmmt","SingleEle_mmmt","SingleMu_mmmt","WWZ_mmmt","WZZ_mmmt","ZZZ_mmmt","ttZJets_mmmt","GGZZ2e2mu_mmmt","GGZZ2e2tau_mmmt","GGZZ2mu2tau_mmmt","GGZZ4tau_mmmt","GGZZ4e_mmmt","GGZZ4mu_mmmt","ZZ4L_mmtt","ZZ4L_v2_mmtt","ZH110_mmtt","ZH120_mmtt","ZH125_mmtt","ZH130_mmtt","ZH140_mmtt","ZHWW125_mmtt","DoubleEle_mmtt","DoubleMu_mmtt","SingleEle_mmtt","SingleMu_mmtt","WWZ_mmtt","WZZ_mmtt","ZZZ_mmtt","ttZJets_mmtt","GGZZ2e2mu_mmtt","GGZZ2e2tau_mmtt","GGZZ2mu2tau_mmtt","GGZZ4tau_mmtt","GGZZ4e_mmtt","GGZZ4mu_mmtt"]

    samples=["DoubleEleB_v1_eeem","DoubleEleB_v2_eeem","DoubleEleC_eeem","DoubleEleD_eeem","DoubleEleE_eeem","DoubleEleF_eeem","DoubleEleG_eeem","DoubleEleH_v2_eeem","DoubleEleH_v3_eeem","SingleEleB_v1_eeem","SingleEleB_v2_eeem","SingleEleC_eeem","SingleEleD_eeem","SingleEleE_eeem","SingleEleF_eeem","SingleEleG_eeem","SingleEleH_v2_eeem","SingleEleH_v3_eeem","DoubleMuB_v1_eeem","DoubleMuB_v2_eeem","DoubleMuC_eeem","DoubleMuD_eeem","DoubleMuE_eeem","DoubleMuF_eeem","DoubleMuG_eeem","DoubleMuH_v2_eeem","DoubleMuH_v3_eeem","SingleMuB_v1_eeem","SingleMuB_v2_eeem","SingleMuC_eeem","SingleMuD_eeem","SingleMuE_eeem","SingleMuF_eeem","SingleMuG_eeem","SingleMuH_v2_eeem","SingleMuH_v3_eeem","DoubleEleB_v1_eeet","DoubleEleB_v2_eeet","DoubleEleC_eeet","DoubleEleD_eeet","DoubleEleE_eeet","DoubleEleF_eeet","DoubleEleG_eeet","DoubleEleH_v2_eeet","DoubleEleH_v3_eeet","SingleEleB_v1_eeet","SingleEleB_v2_eeet","SingleEleC_eeet","SingleEleD_eeet","SingleEleE_eeet","SingleEleF_eeet","SingleEleG_eeet","SingleEleH_v2_eeet","SingleEleH_v3_eeet","DoubleMuB_v1_eeet","DoubleMuB_v2_eeet","DoubleMuC_eeet","DoubleMuD_eeet","DoubleMuE_eeet","DoubleMuF_eeet","DoubleMuG_eeet","DoubleMuH_v2_eeet","DoubleMuH_v3_eeet","SingleMuB_v1_eeet","SingleMuB_v2_eeet","SingleMuC_eeet","SingleMuD_eeet","SingleMuE_eeet","SingleMuF_eeet","SingleMuG_eeet","SingleMuH_v2_eeet","SingleMuH_v3_eeet","DoubleEleB_v1_eemt","DoubleEleB_v2_eemt","DoubleEleC_eemt","DoubleEleD_eemt","DoubleEleE_eemt","DoubleEleF_eemt","DoubleEleG_eemt","DoubleEleH_v2_eemt","DoubleEleH_v3_eemt","SingleEleB_v1_eemt","SingleEleB_v2_eemt","SingleEleC_eemt","SingleEleD_eemt","SingleEleE_eemt","SingleEleF_eemt","SingleEleG_eemt","SingleEleH_v2_eemt","SingleEleH_v3_eemt","DoubleMuB_v1_eemt","DoubleMuB_v2_eemt","DoubleMuC_eemt","DoubleMuD_eemt","DoubleMuE_eemt","DoubleMuF_eemt","DoubleMuG_eemt","DoubleMuH_v2_eemt","DoubleMuH_v3_eemt","SingleMuB_v1_eemt","SingleMuB_v2_eemt","SingleMuC_eemt","SingleMuD_eemt","SingleMuE_eemt","SingleMuF_eemt","SingleMuG_eemt","SingleMuH_v2_eemt","SingleMuH_v3_eemt","DoubleEleB_v1_eett","DoubleEleB_v2_eett","DoubleEleC_eett","DoubleEleD_eett","DoubleEleE_eett","DoubleEleF_eett","DoubleEleG_eett","DoubleEleH_v2_eett","DoubleEleH_v3_eett","SingleEleB_v1_eett","SingleEleB_v2_eett","SingleEleC_eett","SingleEleD_eett","SingleEleE_eett","SingleEleF_eett","SingleEleG_eett","SingleEleH_v2_eett","SingleEleH_v3_eett","DoubleMuB_v1_eett","DoubleMuB_v2_eett","DoubleMuC_eett","DoubleMuD_eett","DoubleMuE_eett","DoubleMuF_eett","DoubleMuG_eett","DoubleMuH_v2_eett","DoubleMuH_v3_eett","SingleMuB_v1_eett","SingleMuB_v2_eett","SingleMuC_eett","SingleMuD_eett","SingleMuE_eett","SingleMuF_eett","SingleMuG_eett","SingleMuH_v2_eett","SingleMuH_v3_eett","DoubleEleB_v1_mmem","DoubleEleB_v2_mmem","DoubleEleC_mmem","DoubleEleD_mmem","DoubleEleE_mmem","DoubleEleF_mmem","DoubleEleG_mmem","DoubleEleH_v2_mmem","DoubleEleH_v3_mmem","SingleEleB_v1_mmem","SingleEleB_v2_mmem","SingleEleC_mmem","SingleEleD_mmem","SingleEleE_mmem","SingleEleF_mmem","SingleEleG_mmem","SingleEleH_v2_mmem","SingleEleH_v3_mmem","DoubleMuB_v1_mmem","DoubleMuB_v2_mmem","DoubleMuC_mmem","DoubleMuD_mmem","DoubleMuE_mmem","DoubleMuF_mmem","DoubleMuG_mmem","DoubleMuH_v2_mmem","DoubleMuH_v3_mmem","SingleMuB_v1_mmem","SingleMuB_v2_mmem","SingleMuC_mmem","SingleMuD_mmem","SingleMuE_mmem","SingleMuF_mmem","SingleMuG_mmem","SingleMuH_v2_mmem","SingleMuH_v3_mmem","DoubleEleB_v1_mmet","DoubleEleB_v2_mmet","DoubleEleC_mmet","DoubleEleD_mmet","DoubleEleE_mmet","DoubleEleF_mmet","DoubleEleG_mmet","DoubleEleH_v2_mmet","DoubleEleH_v3_mmet","SingleEleB_v1_mmet","SingleEleB_v2_mmet","SingleEleC_mmet","SingleEleD_mmet","SingleEleE_mmet","SingleEleF_mmet","SingleEleG_mmet","SingleEleH_v2_mmet","SingleEleH_v3_mmet","DoubleMuB_v1_mmet","DoubleMuB_v2_mmet","DoubleMuC_mmet","DoubleMuD_mmet","DoubleMuE_mmet","DoubleMuF_mmet","DoubleMuG_mmet","DoubleMuH_v2_mmet","DoubleMuH_v3_mmet","SingleMuB_v1_mmet","SingleMuB_v2_mmet","SingleMuC_mmet","SingleMuD_mmet","SingleMuE_mmet","SingleMuF_mmet","SingleMuG_mmet","SingleMuH_v2_mmet","SingleMuH_v3_mmet","DoubleEleB_v1_mmmt","DoubleEleB_v2_mmmt","DoubleEleC_mmmt","DoubleEleD_mmmt","DoubleEleE_mmmt","DoubleEleF_mmmt","DoubleEleG_mmmt","DoubleEleH_v2_mmmt","DoubleEleH_v3_mmmt","SingleEleB_v1_mmmt","SingleEleB_v2_mmmt","SingleEleC_mmmt","SingleEleD_mmmt","SingleEleE_mmmt","SingleEleF_mmmt","SingleEleG_mmmt","SingleEleH_v2_mmmt","SingleEleH_v3_mmmt","DoubleMuB_v1_mmmt","DoubleMuB_v2_mmmt","DoubleMuC_mmmt","DoubleMuD_mmmt","DoubleMuE_mmmt","DoubleMuF_mmmt","DoubleMuG_mmmt","DoubleMuH_v2_mmmt","DoubleMuH_v3_mmmt","SingleMuB_v1_mmmt","SingleMuB_v2_mmmt","SingleMuC_mmmt","SingleMuD_mmmt","SingleMuE_mmmt","SingleMuF_mmmt","SingleMuG_mmmt","SingleMuH_v2_mmmt","SingleMuH_v3_mmmt","DoubleEleB_v1_mmtt","DoubleEleB_v2_mmtt","DoubleEleC_mmtt","DoubleEleD_mmtt","DoubleEleE_mmtt","DoubleEleF_mmtt","DoubleEleG_mmtt","DoubleEleH_v2_mmtt","DoubleEleH_v3_mmtt","SingleEleB_v1_mmtt","SingleEleB_v2_mmtt","SingleEleC_mmtt","SingleEleD_mmtt","SingleEleE_mmtt","SingleEleF_mmtt","SingleEleG_mmtt","SingleEleH_v2_mmtt","SingleEleH_v3_mmtt","DoubleMuB_v1_mmtt","DoubleMuB_v2_mmtt","DoubleMuC_mmtt","DoubleMuD_mmtt","DoubleMuE_mmtt","DoubleMuF_mmtt","DoubleMuG_mmtt","DoubleMuH_v2_mmtt","DoubleMuH_v3_mmtt","SingleMuB_v1_mmtt","SingleMuB_v2_mmtt","SingleMuC_mmtt","SingleMuD_mmtt","SingleMuE_mmtt","SingleMuF_mmtt","SingleMuG_mmtt","SingleMuH_v2_mmtt","SingleMuH_v3_mmtt"]
    #samples = ["DY","DY_v2","DY1","DY2","DY3","DY4","W","W_v2","W1","W2","W2_v2","W3","W3_v2","W4","W4_v2","W4_v3","WW","WZ","ZZ","ST_tW_top","ST_tW_antitop","ST_t_top","ST_t_antitop","DYlow","DY1low","DY2low","WZ1L3Nu","WZ1L1Nu2Q","WZ2L2Q","WW1L1Nu2Q","ZZ2L2Q","WZJLLLNu","ZZ4L","ggH125","VBF125","gghww125","vbfhww125","WGToLNuG_v1","WGToLNuG_v2","WGstarToLNuEE","WGstarToLNuMuMu","Wplus125","Wminus125","ZH125","EWKWminus_v1","EWKWminus_v2","EWKWminus_v3","EWKWplus_v1","EWKWplus_v2","EWKWplus_v3","EWKZ_v1","EWKZ_v2","EWKZ_v3","TT","VV2L2Nu","VV2L2Nu_v2","ggH110","VBF110","Wplus110","Wminus110","ZH110","ggH120","VBF120","Wplus120","Wminus120","ZH120","ggH130","VBF130","Wplus130","Wminus130","ZH130","ggH140","VBF140","Wplus140","Wminus140","ZH140","EWKZNuNu_v1","EWKZNuNu_v2","EWKZNuNu_v3","ttHToTauTau125"] 
    #samples = ["DY","DY_v2","DY1","DY2","DY3","DY4","W","W_v2","W1","W2","W2_v2","W3","W3_v2","W4","W4_v2","W4_v3","ST_tW_top","ST_tW_antitop","ST_t_top","ST_t_antitop","DYlow","DY1low","DY2low","WZ1L3Nu","WZ1L1Nu2Q","WW1L1Nu2Q","ZZ2L2Q","ggH125","VBF125","gghww125","vbfhww125","WGToLNuG_v1","WGToLNuG_v2","WGstarToLNuEE","WGstarToLNuMuMu","Wplus125","Wminus125","ZH125","EWKWminus_v1","EWKWminus_v2","EWKWminus_v3","EWKWplus_v1","EWKWplus_v2","EWKWplus_v3","EWKZ_v1","EWKZ_v2","EWKZ_v3","ggH110","VBF110","Wplus110","Wminus110","ZH110","ggH120","VBF120","Wplus120","Wminus120","ZH120","ggH130","VBF130","Wplus130","Wminus130","ZH130","ggH140","VBF140","Wplus140","Wminus140","ZH140","EWKZNuNu_v1","EWKZNuNu_v2","EWKZNuNu_v3","ttHToTauTau125"]
    #samples = ["ggH110","VBF110","Wplus110","Wminus110","ZH110","ggH120","VBF120","Wplus120","Wminus120","ZH120""ggH130","VBF130","Wplus130","Wminus130","ZH130","ggH140","VBF140","Wplus140","Wminus140","ZH140","EWKZNuNu_v1","EWKZNuNu_v2","EWKZNuNu_v3","ttHToTauTau125"]
    #samples=["DataC","DataD","DataE","DataF","DataG","DataH1","DataH2","DataH3","DataB1","DataB2"]
    #samples=["DataH2"]
    #samples = ["DY","DY_v2","DY1","DY2","DY3","DY4","W","W_v2","W1","W2","W2_v2","W3","W3_v2","W4","W4_v2","W4_v3","WW","WZ","ZZ","ST_tW_top","ST_tW_antitop","ST_t_top","ST_t_antitop","DYlow","DY1low","DY2low","WZ1L3Nu","WZ1L1Nu2Q","ggH125","VBF125","gghww125","vbfhww125","WGToLNuG_v1","WGToLNuG_v2","WGstarToLNuEE","WGstarToLNuMuMu","Wplus125","Wminus125","ZH125","EWKWMinus_v1","EWKWminus_v2","EWKWminus_v3","EWKWplus_v1","EWKWplus_v2","EWKWplus_v3","EWKZ_v1","EWKZ_v2","EWKZ_v3","DataC","DataD","DataE","DataF","DataG","DataH1","DataH2","DataH3","DataB1","DataB2"]#,"TT","VV2L2Nu","VV2L2Nu_v2","WZ2L2Q","WW1L1Nu2Q","ZZ2L2Q","WZJLLLNu","ZZ4L"]
    #samples=["WZ1L1Nu2Q"]

    originalDir = '/nfs_scratch/caillol/ZH_oct30'
    targetDir = '/nfs_scratch/caillol/ZHoct30_merged'
    jobId = ''
    channel = 'xx'
    ttreePath = 'RLE_tree'
    for sample in samples :
	print sample
        mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir )


    #samples = ['WJets', 'WJets1', 'WJets2', 'WJets3', 'WJets4'] # B-F data, add WW, WZ3l1nu, WZ, ZZ4l, ZZ, removed WWZ, WZZ
    #targetDir = '/nfs_scratch/truggles/httOct31svFitPrepMerged/WJets'
    #jobId = ''
    #channel = 'tt'
    #ttreePath = 'tt/final/Ntuple'
    #for sample in samples :
    #    mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir )

    #samples = ['EWKWPlus', 'EWKWMinus', 'EWKZ2l', 'EWKZ2nu', 'WWW', 'ZZZ', 'T-tchan', 'Tbar-tchan', 'TT', 'Tbar-tW', 'T-tW', 'WW1l1nu2q', 'WZ1l1nu2q', 'WZ1l3nu', 'WZ2l2q', 'ZZ2l2q', 'VV', 'dataTT-B', 'dataTT-C', 'dataTT-D', 'dataTT-E', 'dataTT-F', 'WW', 'WZ3l1nu', 'WZ', 'ZZ4l', 'ZZ'] # B-F data, add WW, WZ3l1nu, WZ, ZZ4l, ZZ, removed WWZ, WZZ
    #samples = ['WW', 'WZ3l1nu', 'WZ', 'ZZ4l', 'ZZ'] # B-F data, add WW, WZ3l1nu, WZ, ZZ4l, ZZ, removed WWZ, WZZ
    #targetDir = '/nfs_scratch/truggles/httOct31svFitPrepMerged/AllOthers'
    #jobId = ''
    #channel = 'tt'
    #ttreePath = 'tt/final/Ntuple'
    #for sample in samples :
    #    mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir )


