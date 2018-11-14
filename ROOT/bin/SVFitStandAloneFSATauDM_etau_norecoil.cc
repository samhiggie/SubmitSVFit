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

//If doES       0 does not apply any ES shifts
//              1 applies ES shifts to TT channel, no effect on other channels
//

void copyFiles( optutl::CommandLineParser parser, TFile* fOld, TFile* fNew) ;
void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int doES) ;
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

   parser.parseArguments (argc, argv);

   std::cout << "EXTRA COMMANDS:"
    << "\n --- doES: " << parser.doubleValue("doES") << std::endl;

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
     readdir(fProduce,parser,TreeToUse,parser.doubleValue("doES"));

     fProduce->Close();
     f->Close();
   }
   else{
     TFile *f = new TFile(parser.stringValue("inputFile").c_str(),"UPDATE");
     readdir(f,parser,TreeToUse,parser.doubleValue("doES"));
     f->Close();
   }


} 


void readdir(TDirectory *dir, optutl::CommandLineParser parser, char TreeToUse[], int doES) 
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
      if (std::string(key->GetName()).find("_tree")){ readdir(subdir,parser,TreeToUse,parser.doubleValue("doES"));
      }

      dirsav->cd();
    }
    else if(k<1 && obj->IsA()->InheritsFrom(TTree::Class())) {
	k++;
	std::cout<<"ici!!!!"<<std::endl;

      TTree *t = (TTree*)obj;
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

      TBranch *newBranch1 = t->Branch("m_sv", &svFitMass, "m_sv/F");
      TBranch *newBranch2 = t->Branch("pt_sv", &svFitPt, "pt_sv/F");
      TBranch *newBranch3 = t->Branch("phi_sv", &svFitPhi, "phi_sv/F");
      TBranch *newBranch4 = t->Branch("eta_sv", &svFitEta, "eta_sv/F");
      TBranch *newBranch1D = t->Branch("m_sv_DOWN", &svFitMass_DOWN, "m_sv_DOWN/F");
      TBranch *newBranch2D = t->Branch("pt_sv_DOWN", &svFitPt_DOWN, "pt_sv_DOWN/F");
      TBranch *newBranch3D = t->Branch("phi_sv_DOWN", &svFitPhi_DOWN, "phi_sv_DOWN/F");
      TBranch *newBranch4D = t->Branch("eta_sv_DOWN", &svFitEta_DOWN, "eta_sv_DOWN/F");
      TBranch *newBranch1U = t->Branch("m_sv_UP", &svFitMass_UP, "m_sv_UP/F");
      TBranch *newBranch2U = t->Branch("pt_sv_UP", &svFitPt_UP, "pt_sv_UP/F");
      TBranch *newBranch3U = t->Branch("phi_sv_UP", &svFitPhi_UP, "phi_sv_UP/F");
      TBranch *newBranch4U = t->Branch("eta_sv_UP", &svFitEta_UP, "eta_sv_UP/F");

      int evt;
      int run, lumi;
      float pt1;
      float eta1;
      float phi1;
      float pt2;
      float eta2;
      float phi2;
      float m2;
      int gen_match_2;
      float decayMode2;
      float pfCovMatrix00;
      float pfCovMatrix10;
      float pfCovMatrix01;
      float pfCovMatrix11;

      // define MET
      double measuredMETx = 0.;
      double measuredMETy = 0.;
      float pfmet;
      float pfmetphi;
      TLorentzVector TMet(0,0,0,0);
      TMatrixD covMET(2, 2);
      TBranch *pt1branch;

      t->SetBranchAddress("evt",&evt);
      t->SetBranchAddress("run",&run);
      t->SetBranchAddress("lumi",&lumi);
      t->SetBranchAddress("pt_1",&pt1,&pt1branch);
      t->SetBranchAddress("eta_1",&eta1);
      t->SetBranchAddress("phi_1",&phi1);
      t->SetBranchAddress("pt_2",&pt2);
      t->SetBranchAddress("eta_2",&eta2);
      t->SetBranchAddress("phi_2",&phi2);
      t->SetBranchAddress("m_2",&m2);
      t->SetBranchAddress("gen_match_2",&gen_match_2);
      t->SetBranchAddress("l2_decayMode",&decayMode2);
      t->SetBranchAddress("met",&pfmet);
      t->SetBranchAddress("metphi",&pfmetphi);
      // FOR PF MET ANALYSIS
      t->SetBranchAddress("metcov00",&pfCovMatrix00);
      t->SetBranchAddress("metcov01",&pfCovMatrix01);
      t->SetBranchAddress("metcov10",&pfCovMatrix10);
      t->SetBranchAddress("metcov11",&pfCovMatrix11);

      printf("Found tree -> weighting\n");

      /*int verbosity = 1;
      ClassicSVfit svFitAlgo(verbosity);
      svFitAlgo.addLogM_fixed(true, 6.);
      svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");*/

      for(Int_t i=0;i<t->GetEntries();++i){
         t->GetEntry(i);

         TMet.SetPtEtaPhiM(pfmet,0,pfmetphi,0);
         measuredMETx = pfmet*TMath::Cos(pfmetphi);
         measuredMETy = pfmet*TMath::Sin(pfmetphi);

         covMET[0][0] =  pfCovMatrix00;
         covMET[1][0] =  pfCovMatrix10;
         covMET[0][1] =  pfCovMatrix01;
         covMET[1][1] =  pfCovMatrix11;

	 std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
	 measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, pt1, eta1,  phi1, 0.51100e-3)); 
         measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  pt2, eta2, phi2, m2, decayMode2));
         runSVFit(measuredTauLeptons, measuredMETx, measuredMETy, covMET, 0, svFitMass, svFitPt, svFitEta, svFitPhi);

	 if (doES){
             if (gen_match_2<=5){
                float ES_UP_scale=1.0; // this value is for jet -> tau fakes
                if (gen_match_2<5) ES_UP_scale=1.015; // for gen matched ele/muon
                if (gen_match_2==5) ES_UP_scale=1.008; // for real taus
                double pt2_UP;
                double mass2_UP=m2;
                if (decayMode2!=0) mass2_UP = m2 * ES_UP_scale;
                pt2_UP = pt2 * ES_UP_scale;
                double metcorr_ex_UP, metcorr_ey_UP;
                double dx2_UP, dy2_UP;
                dx2_UP = pt2_UP * TMath::Cos( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
                dy2_UP = pt2_UP * TMath::Sin( phi2 ) * (( 1. / ES_UP_scale ) - 1.);
                metcorr_ex_UP = measuredMETx + dx2_UP;
                metcorr_ey_UP = measuredMETy + dy2_UP;

                std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsUP;
                measuredTauLeptonsUP.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, pt1, eta1,  phi1, 0.51100e-3));
                measuredTauLeptonsUP.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  pt2_UP, eta2, phi2,  mass2_UP, decayMode2));
               runSVFit(measuredTauLeptonsUP, metcorr_ex_UP, metcorr_ey_UP, covMET, 0, svFitMass_UP, svFitPt_UP, svFitEta_UP, svFitPhi_UP);

                float ES_DOWN_scale=1.0; // jet
                if (gen_match_2==5) ES_DOWN_scale=0.992; // tau
                if (gen_match_2<5) ES_DOWN_scale=0.985;  // elec/mu
                double pt2_DOWN;
                double mass2_DOWN = m2;
                pt2_DOWN = pt2 * ES_DOWN_scale;
                if (decayMode2!=0) mass2_DOWN = m2 * ES_DOWN_scale;
                double metcorr_ex_DOWN, metcorr_ey_DOWN;
                double dx2_DOWN, dy2_DOWN;
                dx2_DOWN = pt2_DOWN * TMath::Cos( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
                dy2_DOWN = pt2_DOWN * TMath::Sin( phi2 ) * (( 1. / ES_DOWN_scale ) - 1.);
                metcorr_ex_DOWN = measuredMETx + dx2_DOWN;
                metcorr_ey_DOWN = measuredMETy + dy2_DOWN;

                std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptonsDOWN;
                measuredTauLeptonsDOWN.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToElecDecay, pt1, eta1,  phi1, 0.51100e-3));
                measuredTauLeptonsDOWN.push_back(classic_svFit::MeasuredTauLepton(classic_svFit::MeasuredTauLepton::kTauToHadDecay,  pt2_DOWN, eta2, phi2,  mass2_DOWN, decayMode2));
               runSVFit(measuredTauLeptonsDOWN, metcorr_ex_DOWN, metcorr_ey_DOWN, covMET, 0, svFitMass_DOWN, svFitPt_DOWN, svFitEta_DOWN, svFitPhi_DOWN);

	    }
	    else{
	       svFitMass_DOWN=svFitMass;
               svFitPt_DOWN=svFitPt;
               svFitEta_DOWN=svFitEta;
               svFitPhi_DOWN=svFitPhi;
               svFitMass_UP=svFitMass;
               svFitPt_UP=svFitPt;
               svFitEta_UP=svFitEta;
               svFitPhi_UP=svFitPhi;
	    }
	 }

         newBranch1->Fill();
         newBranch2->Fill();
         newBranch3->Fill();
         newBranch4->Fill();
         newBranch1U->Fill();
         newBranch2U->Fill();
         newBranch3U->Fill();
         newBranch4U->Fill();
         newBranch1D->Fill();
         newBranch2D->Fill();
         newBranch3D->Fill();
         newBranch4D->Fill();
    }
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

