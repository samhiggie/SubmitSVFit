#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include "TFile.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TKey.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include <math.h> 
#include "TMath.h" 
#include <limits>
#include "TSystem.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"


void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int doES, int doUES, int doRES, int doJES, int year) ;
void CopyFile(const char *fname, optutl::CommandLineParser parser);
void CopyDir(TDirectory *source,optutl::CommandLineParser parser);
void runSVFit(std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons, double measuredMETx, double measuredMETy, TMatrixD &covMET, float num, float &svFitMass, float& svFitPt, float &svFitEta, float &svFitPhi);

int main (int argc, char* argv[]) 
{
    optutl::CommandLineParser parser ("Sets Event Weights in the ntuple");
    parser.addOption("branch",optutl::CommandLineParser::kString,"Branch","__svFit__");
    parser.addOption("newFile",optutl::CommandLineParser::kString,"newFile","newFile.root");
    parser.addOption("inputFile",optutl::CommandLineParser::kString,"input File");
    parser.addOption("newOutputFile",optutl::CommandLineParser::kDouble,"New Output File",0.0);
    parser.addOption("doES",optutl::CommandLineParser::kDouble,"doES",0.0);
    parser.addOption("doUES",optutl::CommandLineParser::kDouble,"doUES",0.0);
    parser.addOption("doRES",optutl::CommandLineParser::kDouble,"doRES",0.0);
    parser.addOption("doJES",optutl::CommandLineParser::kDouble,"doJES",0.0);
    parser.addOption("year",optutl::CommandLineParser::kDouble,"year",0.0);

    parser.parseArguments (argc, argv);

    std::cout << "EXTRA COMMANDS:"
        << "\n --- doES: " << parser.doubleValue("doES") 
        << "\n --- doUES: " << parser.doubleValue("doUES")
        << "\n --- doRES: " << parser.doubleValue("doRES") 
        << "\n --- doJES: " << parser.doubleValue("doJES") << std::endl;

    char TreeToUse[80]="first" ;

    TFile *fProduce;//= new TFile(parser.stringValue("newFile").c_str(),"UPDATE");

    if(parser.doubleValue("newOutputFile")>0){
        TFile *f = new TFile(parser.stringValue("inputFile").c_str(),"READ");
        std::cout<<"Creating new outputfile"<<std::endl;
        std::string newFileName = parser.stringValue("newFile");

        fProduce = new TFile(newFileName.c_str(),"RECREATE");
        copyFiles(parser, f, fProduce);//new TFile(parser.stringValue("inputFile").c_str()+"SVFit","UPDATE");
        fProduce = new TFile(newFileName.c_str(),"UPDATE");
        std::cout<<"listing the directories================="<<std::endl;
        fProduce->ls();
        readdir(fProduce,parser,TreeToUse,parser.doubleValue("doES"),parser.doubleValue("doUES"),parser.doubleValue("doRES"),parser.doubleValue("doJES"),parser.doubleValue("year"));

        fProduce->Close();
        f->Close();
    }
    else{
        TFile *f = new TFile(parser.stringValue("inputFile").c_str(),"UPDATE");
        readdir(f,parser,TreeToUse,parser.doubleValue("doES"),parser.doubleValue("doUES"),parser.doubleValue("doRES"),parser.doubleValue("doJES"),parser.doubleValue("year"));
        f->Close();
    }


} 

void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int doES, int doUES, int doRES, int doJES, int year) 
{

    TDirectory *dirsav = gDirectory;
    TIter next(dir->GetListOfKeys());
    TKey *key;
    char stringA[80]="first";
    dir->cd();      
    int k=0;
    while ((key = (TKey*)next())) {
        printf("Found key=%s \n",key->GetName());

        TObject *obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
            dir->cd(key->GetName());
            TDirectory *subdir = gDirectory;
            sprintf(TreeToUse,"%s",key->GetName());
            if (std::string(key->GetName()).find("_tree")){ readdir(subdir,parser,TreeToUse,parser.doubleValue("doES"),parser.doubleValue("doUES"),parser.doubleValue("doRES"),parser.doubleValue("doJES"),parser.doubleValue("year"));
            }

            dirsav->cd();
        }
        else if(k<1 && obj->IsA()->InheritsFrom(TTree::Class())) {
            k++;
            std::cout<<"ici!!!!"<<std::endl;

            TTree *t = (TTree*)obj;

            std::vector<std::string> uncertNames = {
                "AbsoluteFlavMap",
                "AbsoluteMPFBias",
                "AbsoluteScale",
                "AbsoluteStat",
                "AbsoluteSample",
                "FlavorQCD",
                "Fragmentation",
                "PileUpDataMC",
                "PileUpPtBB",
                "PileUpPtEC1",
                "PileUpPtEC2",
                "PileUpPtHF",
                "PileUpPtRef",
                "RelativeBal",
                "RelativeSample",
                "RelativeFSR",
                "RelativeJEREC1",
                "RelativeJEREC2",
                "RelativeJERHF",
                "RelativePtBB",
                "RelativePtEC1",
                "RelativePtEC2",
                "RelativePtHF",
                "RelativeStatEC",
                "RelativeStatFSR",
                "RelativeStatHF",
                "SinglePionECAL",
                "SinglePionHCAL",
                "SubTotalAbsolute",
                "SubTotalMC",
                "SubTotalPileUp",
                "SubTotalPt",
                "SubTotalRelative",
                "SubTotalScale",
                "TimePtEta",
                "TotalNoFlavorNoTime",
                "TotalNoFlavor",
                "TotalNoTime",
                "Total",
                "Eta3to5",
                "Eta0to5",
                "Eta0to3",
                "EC2",
                "Closure",
            }; // end uncertNames
            // uncertNames={"AbsoluteFlavMap", "AbsoluteMPFBias", "AbsoluteScale", "AbsoluteStat", "AbsoluteSample", "FlavorQCD", "Fragmentation", "PileUpDataMC", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF", "PileUpPtRef", "RelativeBal", "RelativeSample","RelativeFSR", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeStatEC", "RelativeStatFSR", "RelativeStatHF", "SinglePionECAL", "SinglePionHCAL", "SubTotalAbsolute", "SubTotalMC", "SubTotalPileUp", "SubTotalPt", "SubTotalRelative", "SubTotalScale", "TimePtEta", "TotalNoFlavorNoTime","TotalNoFlavor","TotalNoTime","Total", "Closure"};

            //variables for calculations
            unsigned long long evt;
            int run, lumi;
            float pt1;
            float eta1;
            float phi1;
            float gen_match_1;
            float pt2;
            float eta2;
            float phi2;
            float m2;
            float gen_match_2;
            float decayMode=-999.;
            float decayMode2;
            float covMet00;
            float covMet10;
            float covMet01;
            float covMet11;
            float pfCovMatrix00;
            float pfCovMatrix10;
            float pfCovMatrix01;
            float pfCovMatrix11;
            float svFitMass = -10;
            float svFitPt = -10;
            float svFitEta = -10;
            float svFitPhi = -10;
            float svFitMass_UP = -10;
            float svFitPt_UP = -10;
            float svFitEta_UP = -10;
            float svFitPhi_UP = -10;
            float svFitMass_DOWN = -10;
            float svFitPt_DOWN = -10;
            float svFitEta_DOWN = -10;
            float svFitPhi_DOWN = -10;
            float svFitMass_UESUp = -10;
            float svFitPt_UESUp = -10;
            float svFitEta_UESUp = -10;
            float svFitPhi_UESUp = -10;
            float svFitMass_UESDown = -10;
            float svFitPt_UESDown = -10;
            float svFitEta_UESDown = -10;
            float svFitPhi_UESDown = -10;
            float svFitMass_MESUp = -10;
            float svFitPt_MESUp = -10;
            float svFitEta_MESUp = -10;
            float svFitPhi_MESUp = -10;
            float svFitMass_MESDown = -10;
            float svFitPt_MESDown = -10;
            float svFitEta_MESDown = -10;
            float svFitPhi_MESDown = -10;

            //float mvamet_ex, // uncorrected mva met px (float)
            //  mvamet_ey, // uncorrected mva met py (float)
            float  genPx=-999.    , // generator Z/W/Higgs px (float)
                   genPy =-999.   , // generator Z/W/Higgs py (float)
                   visPx =-999.   , // generator visible Z/W/Higgs px (float)
                   visPy =-999.   ; // generator visible Z/W/Higgs py (float)

            int njets =-999.   ;  // number of jets (hadronic jet multiplicity) (int)

            // define MET
            double measuredMETx = 0.;
            double measuredMETy = 0.;
            float mvamet;
            float mvametphi;
            float pfmet;
            float pfmetphi;
            TLorentzVector TMet(0,0,0,0);
            // define MET covariance
            TMatrixD covMET(2, 2);

            // MET Uncertainties
            float uncMetPtUp;
            float uncMetPtDown;
            float uncMetPhiUp;
            float uncMetPhiDown;
            float clusteredMetPtUp;
            float clusteredMetPtDown;
            float clusteredMetPhiUp;
            float clusteredMetPhiDown;
            double uncMetUpMETx = 0.;
            double uncMetUpMETy = 0.;
            double uncMetDownMETx = 0.;
            double uncMetDownMETy = 0.;
            double clusteredMetUpMETx = 0.;
            double clusteredMetUpMETy = 0.;
            double clusteredMetDownMETx = 0.;
            double clusteredMetDownMETy = 0.;
            //Branches that we need for the calculations

            //new vector of TBranches for all the uncertainty shapes 
            //std::cout<<"initializing branches"<<std::endl;

            std::vector<int> njetVecUp;
            std::vector<int> njetVecDown;
            std::vector<float> vbfVecUp;
            std::vector<float> vbfVecDown;
            std::vector<float> metVecUp;
            std::vector<float> metVecDown;
            std::vector<float> metphiVecUp;
            std::vector<float> metphiVecDown;
            std::vector<float> jetaVecUp;
            std::vector<float> jetaVecDown;
            std::vector<float> msvVecUp;
            std::vector<float> msvVecDown;
            std::vector<float> ptsvVecUp;
            std::vector<float> ptsvVecDown;
            std::vector<float> etasvVecUp;
            std::vector<float> etasvVecDown;
            std::vector<float> phisvVecUp;
            std::vector<float> phisvVecDown;

            for(unsigned int i=0; i<uncertNames.size(); i++){
                njetVecUp.push_back(-999);
                njetVecDown.push_back(-999);
                vbfVecUp.push_back(-999.);
                vbfVecDown.push_back(-999.);
                metVecUp.push_back(-999.);
                metVecDown.push_back(-999.);
                metphiVecUp.push_back(-999.);
                metphiVecDown.push_back(-999.);
                msvVecUp.push_back(-999.);
                msvVecDown.push_back(-999.);
                ptsvVecUp.push_back(-999.);
                ptsvVecDown.push_back(-999.);
                //jetaVecUp.push_back(-999.);
                //jetaVecDown.push_back(-999.);

            }
        
            for(unsigned int i=0; i<uncertNames.size(); i++){
                //up and down variables for calculations
                t->SetBranchAddress(("njet_"   +uncertNames[i]+"Up").c_str(),  &njetVecUp[i]  );
                t->SetBranchAddress(("njet_"   +uncertNames[i]+"Down").c_str(),&njetVecDown[i] );
                t->SetBranchAddress(("vbfMass_"+uncertNames[i]+"Up").c_str(),  &vbfVecUp[i]  );
                t->SetBranchAddress(("vbfMass_"+uncertNames[i]+"Down").c_str(),&vbfVecDown[i] );
                t->SetBranchAddress(("type1_pfMet_shiftedPt_Jet"    +uncertNames[i]+"Up").c_str(),  &metVecUp[i]);
                t->SetBranchAddress(("type1_pfMet_shiftedPt_Jet"    +uncertNames[i]+"Down").c_str(),&metVecDown[i]);    
                t->SetBranchAddress(("type1_pfMet_shiftedPhi_Jet" +uncertNames[i]+"Up").c_str(),  &metphiVecUp[i]);
                t->SetBranchAddress(("type1_pfMet_shiftedPhi_Jet" +uncertNames[i]+"Down").c_str(),&metphiVecDown[i]);
                //t->SetBranchAddress(("jeta_"   +uncertNames[i]+"Up").c_str(),  &jetaVecUp[i],    ("jeta_"+    uncertNames[i]+"Up/F"  ).c_str());
                //t->SetBranchAddress(("jeta_"   +uncertNames[i]+"Down").c_str(),&jetaVecDown[i],  ("jeta_"+    uncertNames[i]+"Down/F").c_str());

                //new svfit mass up and down variables
                t->Branch(("m_sv_"   +uncertNames[i]+"Up").c_str(),   &msvVecUp[i],    ("m_sv_"+    uncertNames[i]+"Up/F"  ).c_str());
                t->Branch(("m_sv_"   +uncertNames[i]+"Down").c_str(), &msvVecDown[i],  ("m_sv_"+    uncertNames[i]+"Down/F").c_str());
                t->Branch(("pt_sv_"  +uncertNames[i]+"Up").c_str(),   &ptsvVecUp[i],   ("pt_sv_"+   uncertNames[i]+"Up/F"  ).c_str());
                t->Branch(("pt_sv_"   +uncertNames[i]+"Down").c_str(),&ptsvVecDown[i], ("pt_sv_"+   uncertNames[i]+"Down/F").c_str());
                t->Branch(("eta_sv_"  +uncertNames[i]+"Up").c_str(),   &etasvVecUp[i],   ("eta_sv_"+   uncertNames[i]+"Up/F"  ).c_str());
                t->Branch(("eta_sv_"   +uncertNames[i]+"Down").c_str(),&etasvVecDown[i], ("eta_sv_"+   uncertNames[i]+"Down/F").c_str());
                t->Branch(("phi_sv_"  +uncertNames[i]+"Up").c_str(),   &phisvVecUp[i],   ("phi_sv_"+   uncertNames[i]+"Up/F"  ).c_str());
                t->Branch(("phi_sv_"   +uncertNames[i]+"Down").c_str(),&phisvVecDown[i], ("phi_sv_"+   uncertNames[i]+"Down/F").c_str());

            }

            t->Branch("m_sv",            &svFitMass,        "m_sv/F");
            t->Branch("pt_sv",           &svFitPt,          "pt_sv/F");
            t->Branch("phi_sv",          &svFitPhi,         "phi_sv/F");
            t->Branch("eta_sv",          &svFitEta,         "eta_sv/F");
            t->Branch("m_sv_muonESUp",   &svFitMass_MESUp,  "m_sv_muonESUp/F");
            t->Branch("m_sv_muonESDown", &svFitMass_MESDown,"m_sv_muonESDown/F");
            t->Branch("m_sv_DOWN",         &svFitMass_DOWN,     "m_sv_DOWN/F");
            t->Branch("pt_sv_DOWN",        &svFitPt_DOWN,       "pt_sv_DOWN/F");
            t->Branch("phi_sv_DOWN",     &svFitPhi_DOWN,    "phi_sv_DOWN/F");
            t->Branch("eta_sv_DOWN",     &svFitEta_DOWN,    "eta_sv_DOWN/F");

            t->Branch("m_sv_UP",         &svFitMass_UP,     "m_sv_UP/F");
            t->Branch("pt_sv_UP",        &svFitPt_UP,       "pt_sv_UP/F");
            t->Branch("phi_sv_UP",       &svFitPhi_UP,      "phi_sv_UP/F");
            t->Branch("eta_sv_UP",       &svFitEta_UP,      "eta_sv_UP/F");


            t->SetBranchAddress("pt_1",&pt1);
            t->SetBranchAddress("eta_1",&eta1);
            t->SetBranchAddress("phi_1",&phi1);
            t->SetBranchAddress("pt_2",&pt2);
            t->SetBranchAddress("eta_2",&eta2);
            t->SetBranchAddress("phi_2",&phi2);
            t->SetBranchAddress("tMass",&m2);
            //t->SetBranchAddress("l1_decayMode",&decayMode);
            //t->SetBranchAddress("l2_decayMode",&decayMode2);
            t->SetBranchAddress("decayMode_1",&decayMode);
            t->SetBranchAddress("decayMode_2",&decayMode2);
            //t->SetBranchAddress("mvacov00",&covMet00);
            //t->SetBranchAddress("mvacov01",&covMet01);
            //t->SetBranchAddress("mvacov10",&covMet10);
            //t->SetBranchAddress("mvacov11",&covMet11);
            t->SetBranchAddress("metcov00",&covMet00);
            t->SetBranchAddress("metcov01",&covMet01);
            t->SetBranchAddress("metcov10",&covMet10);
            t->SetBranchAddress("metcov11",&covMet11);
            //t->SetBranchAddress("mvamet",&mvamet);
            //t->SetBranchAddress("met",&mvamet);
            //t->SetBranchAddress("mvametphi",&mvametphi);
            //t->SetBranchAddress("metphi",&mvametphi);
            t->SetBranchAddress("njets", &njets);
            t->SetBranchAddress("type1_pfMetEt",&pfmet);
            t->SetBranchAddress("type1_pfMetPhi",&pfmetphi);

            t->SetBranchAddress("gen_match2",&gen_match_2);


            //extra variables for easy reading 
            float metup, metphiup;
            float metdown, metphidown;

            //loop over events

            //for(Int_t i=0;i<t->GetEntries();++i){
            for(Int_t i=0;i<100;++i){
                t->GetEntry(i);


                //edm::Handle<math::Error<2>::type> covHandle;
                TMatrixD covMET(2, 2);
                covMET[0][0] = covMet00;
                covMET[1][0] = covMet10;
                covMET[0][1] = covMet10; // (1,0) is the only one saved
                covMET[1][1] = covMet11;



                std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;

                measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay,pt1,eta1, phi1, 0.10566)); 
                measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay, pt2, eta2, phi2, m2, decayMode2));



            //JET Energy Scale Shifts ... JEC propagated to MET then we calculate SVFit

            for(unsigned int i=0; i<uncertNames.size(); i++){

 
                runSVFit(measuredTauLeptons, metVecUp[i]*TMath::Cos(metphiVecUp[i]), metVecUp[i]*TMath::Sin(metphiVecUp[i]), covMET, 0, msvVecUp[i], ptsvVecUp[i], etasvVecUp[i], phisvVecUp[i]);

                runSVFit(measuredTauLeptons, metVecDown[i]*TMath::Cos(metphiVecDown[i]), metVecDown[i]*TMath::Sin(metphiVecDown[i]), covMET, 0, msvVecDown[i], ptsvVecDown[i], etasvVecDown[i], phisvVecDown[i]);


        

            }



            //////////////now for the ES shifts

            float met = pfmet;
            float metphi = pfmetphi;

            float ES_UP_scale=1.0; // this value is for jet -> tau fakes
            if (eta1<-2.1) ES_UP_scale=1.027;
            else if (eta1>-2.1 && eta1<-1.2) ES_UP_scale=1.009;
            else if (eta1>-1.2 && eta1<1.2) ES_UP_scale=1.004;
            else if (eta1>1.2  && eta1<2.1) ES_UP_scale=1.009;
            else ES_UP_scale=1.017;

            double pt1_UP;
            pt1_UP = pt1 * ES_UP_scale;
            double metcorr_ex_UP, metcorr_ey_UP;
            double dx1_UP, dy1_UP;
            dx1_UP = pt1_UP * TMath::Cos( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
            dy1_UP = pt1_UP * TMath::Sin( phi1 ) * (( 1. / ES_UP_scale ) - 1.);
            metcorr_ex_UP = met*TMath::Cos(metphi) + dx1_UP;
            metcorr_ey_UP = met*TMath::Sin(metphi) + dy1_UP;

            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;

            measuredTauLeptonsUP.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt1_UP, eta1,  phi1, 0.10566));
            measuredTauLeptonsUP.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  pt2, eta2, phi2,  m2, decayMode2));

            runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMass_MESUp, svFitPt_MESUp, svFitEta_MESUp, svFitPhi_MESUp);

            //met*TMath::Cos(metphi), met*TMath::Sin(metphi)
            float ES_DOWN_scale=1.0; // jet
            if (eta1<-2.1) ES_DOWN_scale=1 - 0.027;
            else if (eta1<-1.2) ES_DOWN_scale=1-0.009;
            else if (eta1<1.2) ES_DOWN_scale=1-0.004;
            else if (eta1<2.1) ES_DOWN_scale=1-0.009;
            else ES_DOWN_scale=1-0.017;
            double pt1_DOWN;
            pt1_DOWN = pt1 * ES_DOWN_scale;
            double metcorr_ex_DOWN, metcorr_ey_DOWN;
            double dx1_DOWN, dy1_DOWN;
            dx1_DOWN = pt1_DOWN * TMath::Cos( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
            dy1_DOWN = pt1_DOWN * TMath::Sin( phi1 ) * (( 1. / ES_DOWN_scale ) - 1.);
            metcorr_ex_DOWN = met*TMath::Cos(metphi) + dx1_DOWN;
            metcorr_ey_DOWN =  met*TMath::Sin(metphi) + dy1_DOWN;

            std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;

            measuredTauLeptonsDOWN.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt1_DOWN, eta1,  phi1, 0.10566));
            measuredTauLeptonsDOWN.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  pt2, eta2, phi2,  m2, decayMode2));

            runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMass_MESDown, svFitPt_MESDown, svFitEta_MESDown, svFitPhi_MESDown);

            if (gen_match_2<=5){
                float ES_UP_scale=1.0; // this value is for jet -> tau fakes
                if (gen_match_2<5) ES_UP_scale=1.03; // for gen matched ele/muon
                /*
                   if (year==2016){
                   if (gen_match_2==5 && decayMode2==0) ES_UP_scale=1.010; // for real taus
                   if (gen_match_2==5 && decayMode2==0) ES_UP_scale=1.009; // for real taus
                   if (gen_match_2==5 && decayMode2==0) ES_UP_scale=1.011; // for real taus
                   }
                   else if (year==2017){

                   if (gen_match_2==5 && decayMode2==0) ES_UP_scale=1.008; // for real taus
                   if (gen_match_2==5 && decayMode2==0) ES_UP_scale=1.008; // for real taus
                   if (gen_match_2==5 && decayMode2==0) ES_UP_scale=1.009; // for real taus
                   }
                   else if (year==2018){
                   */
                if (gen_match_2==5 && decayMode2==0) ES_UP_scale=1.011; // for real taus
                if (gen_match_2==5 && decayMode2==0) ES_UP_scale=1.008; // for real taus
                if (gen_match_2==5 && decayMode2==0) ES_UP_scale=1.009; // for real taus
                //}

                double pt2_UP;
                double mass2_UP=m2;
                if (decayMode2!=0) mass2_UP = m2 * ES_UP_scale;
                pt2_UP = pt2 * ES_UP_scale;
                double metcorr_ex_UP, metcorr_ey_UP;
                double dx2_UP, dy2_UP;
                dx2_UP = pt2_UP * TMath::Cos( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
                dy2_UP = pt2_UP * TMath::Sin( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
                metcorr_ex_UP = met*TMath::Cos(metphi) + dx2_UP;
                metcorr_ey_UP =  met*TMath::Sin(metphi) + dy2_UP;

                std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;

                measuredTauLeptonsUP.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt1, eta1,  phi1, 0.10566));
                measuredTauLeptonsUP.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  pt2_UP, eta2, phi2,  mass2_UP, decayMode2));

                runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMass_UP, svFitPt_UP, svFitEta_UP, svFitPhi_UP);

                float ES_DOWN_scale=1.0; // jet
                if (gen_match_2<5) ES_DOWN_scale=0.97;  // elec/mu
                /*
                   if (year==2016){
                   if (gen_match_2==5 && decayMode2==0) ES_DOWN_scale=0.990; // for real taus
                   if (gen_match_2==5 && decayMode2==0) ES_DOWN_scale=0.991; // for real taus
                   if (gen_match_2==5 && decayMode2==0) ES_DOWN_scale=0.989; // for real taus
                   }
                   else if (year==2017){
                   if (gen_match_2==5 && decayMode2==0) ES_DOWN_scale=0.992; // for real taus
                   if (gen_match_2==5 && decayMode2==0) ES_DOWN_scale=0.992; // for real taus
                   if (gen_match_2==5 && decayMode2==0) ES_DOWN_scale=0.991; // for real taus
                   }
                   else if (year==2018){
                   */
                if (gen_match_2==5 && decayMode2==0) ES_DOWN_scale=0.989; // for real taus
                if (gen_match_2==5 && decayMode2==0) ES_DOWN_scale=0.992; // for real taus
                if (gen_match_2==5 && decayMode2==0) ES_DOWN_scale=0.991; // for real taus
                //}
                double pt2_DOWN;
                double mass2_DOWN = m2;
                pt2_DOWN = pt2 * ES_DOWN_scale;
                if (decayMode2!=0) mass2_DOWN = m2 * ES_DOWN_scale;
                double metcorr_ex_DOWN, metcorr_ey_DOWN;
                double dx2_DOWN, dy2_DOWN;
                dx2_DOWN = pt2_DOWN * TMath::Cos( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
                dy2_DOWN = pt2_DOWN * TMath::Sin( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
                metcorr_ex_DOWN = met*TMath::Cos(metphi) + dx2_DOWN;
                metcorr_ey_DOWN =  met*TMath::Sin(metphi) + dy2_DOWN;

                std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;

                measuredTauLeptonsDOWN.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt1, eta1,  phi1, 0.10566));

                measuredTauLeptonsDOWN.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  pt2_DOWN, eta2, phi2,  mass2_DOWN, decayMode2));

                runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMass_DOWN, svFitPt_DOWN, svFitEta_DOWN, svFitPhi_DOWN);

            }

            //nominal case
            //float met = pfmet;
            //float metphi = pfmetphi;

            //now for SVFit!
            //std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
            measuredTauLeptons.clear();
            measuredTauLeptonsUP.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToMuDecay, pt1, eta1,  phi1, 0.10566));
            measuredTauLeptonsUP.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  pt2, eta2, phi2,  m2, decayMode2));
            //measured leptons, measured met x, measured met y, covmet, 0, return values
            float t_svFitMass, t_svFitPt, t_svFitEta, t_svFitPhi;
            runSVFit(measuredTauLeptons, met*TMath::Cos(metphi), met*TMath::Sin(metphi), covMET, 0, t_svFitMass, t_svFitPt, t_svFitEta, t_svFitPhi);
            svFitMass = t_svFitMass; 
            svFitPt = t_svFitPt; 
            svFitPhi = t_svFitPhi; 
            svFitEta = t_svFitEta;
            //std::cout <<"nominal svFit Mass: "<<svFitMass<< " Pt: "<<svFitPt<< " Phi: "<<svFitPhi<< " Eta: "<<svFitEta<< " MET: "<<met<< " metphi: "<< metphi<<" NLeptons: "<< measuredTauLeptons.size()<<std::endl;




            t->Fill();
            }//end entries
            dir->cd();
            t->Write("",TObject::kOverwrite);
            strcpy(TreeToUse,stringA) ;

        }
    }
}


void runSVFit(std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons, double measuredMETx, double measuredMETy, TMatrixD &covMET, float num, float &svFitMass, float& svFitPt, float &svFitEta, float &svFitPhi){
  FastMTT aFastMTTAlgo;
  aFastMTTAlgo.run(measuredTauLeptons,measuredMETx,measuredMETy,covMET);
  LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
  svFitMass = ttP4.M(); // return value is in units of GeV
  svFitPt = ttP4.Pt();
  svFitEta = ttP4.Eta();
  svFitPhi = ttP4.Phi();
  //std::cout << "found mass = " << svFitMass << std::endl;

}



void CopyDir(TDirectory *source, optutl::CommandLineParser parser) {
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir; 
  if(source->GetName()!=parser.stringValue("inputFile")){
    adir = savdir->mkdir(source->GetName());
    std::cout<<"Source name is not outputfile name"<<std::endl;
    adir->cd();    
  }
  else{
    adir->cd();    
  }

  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    std::cout<<"My key is: "<<key->GetName()<<std::endl;
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir,parser);
      adir->cd();
    } else if (cl->InheritsFrom(TTree::Class())) {
      TTree *T = (TTree*)source->Get(key->GetName());
      adir->cd();
      TTree *newT = T->CloneTree(-1,"fast");
      newT->Write();
    } else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}
void CopyFile(const char *fname, optutl::CommandLineParser parser) {
  //Copy all objects and subdirs of file fname as a subdir of the current directory
  TDirectory *target = gDirectory;
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("Cannot copy file: %s\n",fname);
    target->cd();
    return;
  }
  target->cd();
  CopyDir(f,parser);
  delete f;
  target->cd();
}
void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) 
{
  //prepare files to be copied
  if(gSystem->AccessPathName(parser.stringValue("inputFile").c_str())) {
    gSystem->CopyFile("hsimple.root", parser.stringValue("inputFile").c_str());
  }

  fNew->cd();
  CopyFile(parser.stringValue("inputFile").c_str(),parser);
  fNew->ls();
  fNew->Close();

}
