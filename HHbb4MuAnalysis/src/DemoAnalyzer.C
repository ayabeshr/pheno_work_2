#define DemoAnalyzer_cxx

R__LOAD_LIBRARY(libPhysics)
//R__LOAD_LIBRARY(/home/aya/programs/Delphes-3.4.2/libDelphes.so)
R__ADD_LIBRARY_PATH(/home/aya/programs/Delphes-3.4.2)
R__LOAD_LIBRARY(libDelphes)

//#include "/HEP_DATA/aya/DemoAnalyzer.h"
#include "/home/aya/DemoAnalyzer.h"
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>  
#include <TMath.h>
#include "THStack.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>


void DemoAnalyzer::Loop()
{
//   In a ROOT session, you can do:
//      root> .L DemoAnalyzer.C
//      root> DemoAnalyzer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   
   
   // TH1D
   // various muon combinations bet 1,2,3,4 muons gives various Z masses, choose combination wich is nearer to Z mass 91 GeV to be Za, and the other with mass < 91 GeV take it to be Zb of signal and call them Za, Zb
   
   // Gen_Particles
   TH1D *h_Gen_Muons_pt;
   TH1D *h_Gen_Muons_eta;
   TH1D *h_Gen_Muons_phi;
   TH1D *h_Gen_Bjet_pt;
   TH1D *h_Gen_Bjet_eta;
   TH1D *h_Gen_Bjet_phi;	

   // size of available objects
   TH1D *h_muons_size_Loose;
   TH1D *h_muons_size_Tight;
   TH1D *h_jets_size;
   TH1D *h_MET_size;
   
   // MET    
   TH1D *h_MET;
   TH1D *h_eta_MET;
   TH1D *h_phi_MET;
   
   // All Muons   
   TH1D *h_pt_allMuons_Loose;
   TH1D *h_eta_allMuons_Loose;
   TH1D *h_phi_allMuons_Loose;
   TH1D *h_pt_allMuons_Tight;
   TH1D *h_eta_allMuons_Tight;
   TH1D *h_phi_allMuons_Tight;
   
   // All Jets
   TH1D *h_pt_allJets;
   TH1D *h_eta_allJets;
   TH1D *h_phi_allJets;
	
   // All B-Jets	
   TH1D *h_pt_all_B_Jets;
   TH1D *h_eta_all_B_Jets;
   TH1D *h_phi_all_B_Jets;
   TH1D *h_bjet_size; 
   
   TH1D *h_mZ12;
   TH1D *h_mZ34;
   TH1D *h_mZ13;
   TH1D *h_mZ24;
   TH1D *h_mZ14;
   TH1D *h_mZ23;
      
   TH1D *h_DR_mu1mu2;     // define those histo
   TH1D *h_DR_mu3mu4;
   TH1D *h_DR_mu1mu3;
   TH1D *h_DR_mu2mu4;
   TH1D *h_DR_mu1mu4;
   TH1D *h_DR_mu2mu3;
      
   TH1D *h_mZa_4mu;             // Za and Zb are Z from h1->Za Zb
   TH1D *h_mZb_4mu;
   TH1D *h_pt_Za;
   TH1D *h_pt_Zb;
   TH1D *h_eta_Za;
   TH1D *h_eta_Zb;
   TH1D *h_phi_Za;
   TH1D *h_phi_Zb;
      
   // h1 -> Za Zb
   TH1D *h_mh1_ZaZb;
   TH1D *h_pt_h1_ZaZb;
   TH1D *h_eta_h1_ZaZb;
   TH1D *h_phi_h1_ZaZb;
      
   TH1D *h_mb_jet_1;
   TH1D *h_mb_jet_2;
   TH1D *h_pt_b_jet_1;
   TH1D *h_pt_b_jet_2;
   TH1D *h_eta_b_jet_1;
   TH1D *h_eta_b_jet_2;
   TH1D *h_phi_b_jet_1;
   TH1D *h_phi_b_jet_2;
   TH1D *h_DR_b1b2 ;
   
   // h1 -> b b~
   TH1D *h_mh1_b1b2;
   TH1D *h_pt_h1_b1b2;
   TH1D *h_eta_h1_b1b2;
   TH1D *h_phi_h1_b1b2;
   
   // h2 -> h1 h1 
   TH1D *h_mh2_h1h1;
   TH1D *h_pt_h2_h1h1;
   TH1D *h_eta_h2_h1h1;
   TH1D *h_phi_h2_h1h1;

   // Gen Particles	
   h_Gen_Muons_pt = new TH1D("h_Gen_Muons_pt", "h_Gen_Muons_pt", 500, 0., 500.);
   h_Gen_Muons_pt->GetXaxis()->SetTitle("p_{T} [GeV/C]");
   h_Gen_Muons_pt->GetYaxis()->SetTitle("Number of Events");
   
   h_Gen_Muons_eta = new TH1D("h_Gen_Muons_eta", "h_Gen_Muons_eta", 20, -10., 10.);
   h_Gen_Muons_eta->GetXaxis()->SetTitle("#eta");
   h_Gen_Muons_eta->GetYaxis()->SetTitle("Number of Events");
   
   h_Gen_Muons_phi = new TH1D("h_Gen_Muons_phi", "h_Gen_Muons_phi", 20, -10., 10.);
   h_Gen_Muons_phi->GetXaxis()->SetTitle("#phi");
   h_Gen_Muons_phi->GetYaxis()->SetTitle("Number of Events");
   
   h_Gen_Bjet_pt = new TH1D("h_Gen_Bjet_pt", "h_Gen_Bjet_pt", 500, 0., 500.);
   h_Gen_Bjet_pt->GetXaxis()->SetTitle("p_{T} [GeV/C]");
   h_Gen_Bjet_pt->GetYaxis()->SetTitle("Number of Events");
   
   h_Gen_Bjet_eta = new TH1D("h_Gen_Bjet_eta", "h_Gen_Bjet_eta", 20, -10., 10.);
   h_Gen_Bjet_eta->GetXaxis()->SetTitle("#eta");
   h_Gen_Bjet_eta->GetYaxis()->SetTitle("Number of Events");
   
   h_Gen_Bjet_phi = new TH1D("h_Gen_Bjet_phi", "h_Gen_Bjet_phi", 20, -10., 10.);
   h_Gen_Bjet_phi->GetXaxis()->SetTitle("#phi");
   h_Gen_Bjet_phi->GetYaxis()->SetTitle("Number of Events");	
      
   // Size of all Muons
   h_muons_size_Loose = new TH1D("h_muonsloose_size", "h_muonsloose_size", 10, 0., 10.);
   h_muons_size_Loose->GetXaxis()->SetTitle("Number of Muons");
   h_muons_size_Loose->GetYaxis()->SetTitle("Number of Events"); 
   
   h_muons_size_Tight = new TH1D("h_muonstight_size", "h_muonstight_size", 10, 0., 10.);
   h_muons_size_Tight->GetXaxis()->SetTitle("Number of Muons");
   h_muons_size_Tight->GetYaxis()->SetTitle("Number of Events"); 
   
   // Size of all Jets
   h_jets_size = new TH1D("h_jets_size", "h_jets_size", 10, 0., 10.);
   h_jets_size->GetXaxis()->SetTitle("Number of Jets");
   h_jets_size->GetYaxis()->SetTitle("Number of Events"); 
   
   // Size of MET 
   h_MET_size = new TH1D("h_MET_size", "h_MET_size", 10, 0., 10.);
   h_MET_size->GetXaxis()->SetTitle("Number of MET");
   h_MET_size->GetYaxis()->SetTitle("Number of Events"); 
   
   // MET
   h_MET = new TH1D("h_MET", "h_MET", 100, 0., 200.);
   h_MET->GetXaxis()->SetTitle("MET[GeV]");
   h_MET->GetYaxis()->SetTitle("Number of Events");
   
   // Eta MET
   h_eta_MET = new TH1D("h_eta_MET", "h_eta_MET", 16, -8., 8.);
   h_eta_MET->GetXaxis()->SetTitle("#eta");
   h_eta_MET->GetYaxis()->SetTitle("Number of Events");
   
   // Phi MET 
   h_phi_MET = new TH1D("h_phi_MET", "h_phi_MET", 20, -10., 10.);
   h_phi_MET->GetXaxis()->SetTitle("#phi");
   h_phi_MET->GetYaxis()->SetTitle("Number of Events");
   
   // pT of all Muons
   h_pt_allMuons_Loose = new TH1D("h_pt_allMuonsLoose", "h_pt_allMuonsLoose", 500, 0., 500.);
   h_pt_allMuons_Loose->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_allMuons_Loose->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_allMuons_Tight = new TH1D("h_pt_allMuonsTight", "h_pt_allMuonsTight", 500, 0., 500.);
   h_pt_allMuons_Tight->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_allMuons_Tight->GetYaxis()->SetTitle("Number of Events");
   
   // Eta of all Muons
   h_eta_allMuons_Loose = new TH1D("h_eta_allMuonsLoose", "h_eta_allMuonsLoose", 16, -8., 8.);
   h_eta_allMuons_Loose->GetXaxis()->SetTitle("#eta");
   h_eta_allMuons_Loose->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_allMuons_Tight = new TH1D("h_eta_allMuonsTight", "h_eta_allMuonsTight", 16, -8., 8.);
   h_eta_allMuons_Tight->GetXaxis()->SetTitle("#eta");
   h_eta_allMuons_Tight->GetYaxis()->SetTitle("Number of Events");
   
   // Phi of all Muons
   h_phi_allMuons_Loose= new TH1D("h_phi_allMuonsLoose", "h_phi_allMuonsLoose", 20, -10., 10.);
   h_phi_allMuons_Loose->GetXaxis()->SetTitle("#phi");
   h_phi_allMuons_Loose->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_allMuons_Tight= new TH1D("h_phi_allMuonsTight", "h_phi_allMuonsTight", 20, -10., 10.);
   h_phi_allMuons_Tight->GetXaxis()->SetTitle("#phi");
   h_phi_allMuons_Tight->GetYaxis()->SetTitle("Number of Events");

   // pT of all Jets
   h_pt_allJets = new TH1D("h_pt_allJets", "h_pt_allJets", 500, 0., 500.);
   h_pt_allJets->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_allJets->GetYaxis()->SetTitle("Number of Events");
   
   // Eta of all Jets
   h_eta_allJets = new TH1D("h_eta_allJets", "h_eta_allJets", 16, -8., 8.);
   h_eta_allJets->GetXaxis()->SetTitle("#eta");
   h_eta_allJets->GetYaxis()->SetTitle("Number of Events");
   
   // Phi of all Jets
   h_phi_allJets= new TH1D("h_phi_allJets", "h_phi_allJets", 20, -10., 10.);
   h_phi_allJets->GetXaxis()->SetTitle("#phi");
   h_phi_allJets->GetYaxis()->SetTitle("Number of Events");
	
   // pT of all B-Jets
   h_pt_all_B_Jets = new TH1D("h_pt_all_B_Jets", "h_pt_all_B_Jets", 500, 0., 500.);
   h_pt_all_B_Jets->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_all_B_Jets->GetYaxis()->SetTitle("Number of Events");
   
   // Eta of all B-Jets
   h_eta_all_B_Jets = new TH1D("h_eta_all_B_Jets", "h_eta_all_B_Jets", 16, -8., 8.);
   h_eta_all_B_Jets->GetXaxis()->SetTitle("#eta");
   h_eta_all_B_Jets->GetYaxis()->SetTitle("Number of Events");
   
   // Phi of all B-Jets
   h_phi_all_B_Jets = new TH1D("h_phi_all_B_Jets", "h_phi_all_B_Jets", 20, -10., 10.);
   h_phi_all_B_Jets->GetXaxis()->SetTitle("#phi");
   h_phi_all_B_Jets->GetYaxis()->SetTitle("Number of Events");
   
   h_bjet_size = new TH1D("h_bjet_size", "h_bjet_size", 20, 0., 20.);
   h_bjet_size->GetXaxis()->SetTitle("Number of B_Jets");
   h_bjet_size->GetYaxis()->SetTitle("Number of Events");	
	

   // Combinations 1234
   
   // Z12
   h_mZ12 = new TH1D("h_mZ12", "h_mZ12", 75, 0., 150.);
   h_mZ12->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ12->GetYaxis()->SetTitle("Number of Events");
   
   // Z34
   h_mZ34 = new TH1D("h_mZ34", "h_mZ34", 75, 0., 150.);
   h_mZ34->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ34->GetYaxis()->SetTitle("Number of Events");
   
   // Combinations 1324
   
   // Z13
   h_mZ13 = new TH1D("h_mZ13", "h_mZ13", 75, 0., 150.);
   h_mZ13->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ13->GetYaxis()->SetTitle("Number of Events");
   
   // Z24
   h_mZ24 = new TH1D("h_mZ24", "h_mZ24", 75, 0., 150.);
   h_mZ24->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ24->GetYaxis()->SetTitle("Number of Events");
   
   // Combinations 1423
   
   // Z14
   h_mZ14 = new TH1D("h_mZ14", "h_mZ14", 75, 0., 150.);
   h_mZ14->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ14->GetYaxis()->SetTitle("Number of Events");
   
   // Z23
   h_mZ23 = new TH1D("h_mZ23", "h_mZ23", 75, 0., 150.);
   h_mZ23->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZ23->GetYaxis()->SetTitle("Number of Events");
   
   // Za 
   
   // Za Invariant Mass of dimuons closest to Z mass ~ 91. GeV
   h_mZa_4mu = new TH1D("h_mZa_4mu", "h_mZa_4mu_closest to Z Mass", 120, 0., 120.);
   h_mZa_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZa_4mu->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_Za = new TH1D("h_pt_Za", "h_pt_Za", 500, 0., 500.); 
   h_pt_Za->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_Za->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_Za = new TH1D("h_eta_Za", "h_eta_Za", 16, -8., 8.); 
   h_eta_Za->GetXaxis()->SetTitle("#eta");
   h_eta_Za->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_Za = new TH1D("h_phi_Za", "h_phi_Za", 20, -10., 10.); 
   h_phi_Za->GetXaxis()->SetTitle("#phi");
   h_phi_Za->GetYaxis()->SetTitle("Number of Events");
     
   // Zb
   
   // Za Invariant Mass of dimuons not close to Z mass ~ 91. GeV 
   h_mZb_4mu = new TH1D("h_mZb_4mu", "h_mZb_4mu_Not Close to Z Mass", 120, 0., 120.);
   h_mZb_4mu->GetXaxis()->SetTitle("Invariant mass of dimuon [GeV/c^{2}]"); 
   h_mZb_4mu->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_Zb = new TH1D("h_pt_Zb", "h_pt_Zb", 500, 0., 500.); 
   h_pt_Zb->GetXaxis()->SetTitle("p_{T} [GeV/c]");
   h_pt_Zb->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_Zb = new TH1D("h_eta_Zb", "h_eta_Zb", 16, -8., 8.); 
   h_eta_Zb->GetXaxis()->SetTitle("#eta");
   h_eta_Zb->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_Zb = new TH1D("h_phi_Zb", "h_phi_Zb", 20, -10., 10.); 
   h_phi_Zb->GetXaxis()->SetTitle("#phi");
   h_phi_Zb->GetYaxis()->SetTitle("Number of Events");
     
   h_mb_jet_1 = new TH1D("h_mb_jet_1", "h_mb_jet_1", 100, 0., 100.);
   h_mb_jet_1->GetXaxis()->SetTitle(" mass of b1 jet [GeV/C^{2}]");
   h_mb_jet_1->GetYaxis()->SetTitle("Number of Events");
      
   h_mb_jet_2 = new TH1D("h_mb_jet_2", "h_mb_jet_2", 100, 0., 100.);
   h_mb_jet_2->GetXaxis()->SetTitle(" mass of b2 jet [GeV/C^{2}]");
   h_mb_jet_2->GetYaxis()->SetTitle("Number of Events");
      
   h_pt_b_jet_1 = new TH1D("h_pt_b_jet_1", "h_pt_b_jet_1", 500, 0., 500.);
   h_pt_b_jet_1->GetXaxis()->SetTitle("p_{T} (GeV/C)");
   h_pt_b_jet_1->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_b_jet_2 = new TH1D("h_pt_b_jet_2", "h_pt_b_jet_2", 500, 0., 500.);
   h_pt_b_jet_2->GetXaxis()->SetTitle("p_{T} (GeV/C)");
   h_pt_b_jet_2->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_b_jet_1 = new TH1D("h_eta_b_jet_1", "h_eta_b_jet_1", 20, -8., 8.);
   h_eta_b_jet_1->GetXaxis()->SetTitle("#eta");
   h_eta_b_jet_1->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_b_jet_2 = new TH1D("h_eta_b_jet_2", "h_eta_b_jet_2", 20, -8., 8.);
   h_eta_b_jet_2->GetXaxis()->SetTitle("#eta");
   h_eta_b_jet_2->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_b_jet_1 = new TH1D("h_phi_b_jet_1", "h_phi_b_jet_1", 20, -10., 10.);
   h_phi_b_jet_1->GetXaxis()->SetTitle("#phi");
   h_phi_b_jet_1->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_b_jet_2 = new TH1D("h_phi_b_jet_2", "h_phi_b_jet_2", 20, -10., 10.);
   h_phi_b_jet_2->GetXaxis()->SetTitle("#phi");
   h_phi_b_jet_2->GetYaxis()->SetTitle("Number of Events");
     
   h_DR_b1b2 = new TH1D("h_DR_b1b2", "DR between b1 b2 jets", 20, 0., 20.);
   h_DR_b1b2->GetXaxis()->SetTitle("#Delta R");
   h_DR_b1b2->GetYaxis()->SetTitle("Number of Events");
   
   // DR seperation bet 4Muons  
   h_DR_mu1mu2 = new TH1D("h_DR_mu1mu2", "DR between muon1 muon2", 20, 0., 20.);
   h_DR_mu1mu2->GetXaxis()->SetTitle("#Delta R_{12}");
   h_DR_mu1mu2->GetYaxis()->SetTitle("Number of Events");
      
   h_DR_mu3mu4 = new TH1D("h_DR_mu3mu4", "DR between muon3 muon4", 20, 0., 10.);
   h_DR_mu3mu4->GetXaxis()->SetTitle("#Delta R_{34}");
   h_DR_mu3mu4->GetYaxis()->SetTitle("Number of Events");
      
   h_DR_mu1mu3 = new TH1D("h_DR_mu1mu3", "DR between muon1 muon3", 20, 0., 10.);
   h_DR_mu1mu3->GetXaxis()->SetTitle("#Delta R_{13}");
   h_DR_mu1mu3->GetYaxis()->SetTitle("Number of Events");
      
   h_DR_mu2mu4 = new TH1D("h_DR_mu2mu4", "DR between muon2 muon4", 20, 0., 10.);
   h_DR_mu2mu4->GetXaxis()->SetTitle("#Delta R_{24}");
   h_DR_mu2mu4->GetYaxis()->SetTitle("Number of Events");
   
   h_DR_mu1mu4 = new TH1D("h_DR_mu1mu4", "DR between muon1 muon4", 20, 0., 10.);
   h_DR_mu1mu4->GetXaxis()->SetTitle("#Delta R_{14}");
   h_DR_mu1mu4->GetYaxis()->SetTitle("Number of Events");
   
   h_DR_mu2mu3 = new TH1D("h_DR_mu2mu3", "DR between muon2 muon3", 20, 0., 10.);
   h_DR_mu2mu3->GetXaxis()->SetTitle("#Delta R_{23}");
   h_DR_mu2mu3->GetYaxis()->SetTitle("Number of Events");
   
   
   // DR seperation bet 2 b jets 
   h_DR_b1b2 = new TH1D("h_DR_b1b2", "h_DR_b1b2", 20, 0., 10.); 
   h_DR_b1b2->GetXaxis()->SetTitle("#Delta R");
   h_DR_b1b2->GetYaxis()->SetTitle("Number of Events");
   
   // h1 -> Za Zb
   h_mh1_ZaZb = new TH1D("h_mh1_ZaZb", "h_mh1_ZaZb", 100, 0., 200.);
   h_mh1_ZaZb->GetXaxis()->SetTitle("Invariant mass of 4 muons (GeV/c^{2})");
   h_mh1_ZaZb->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_h1_ZaZb = new TH1D("h_pt_h1_ZaZb", "h_pt_h1_ZaZb", 500, 0., 500.);
   h_pt_h1_ZaZb->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_pt_h1_ZaZb->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_h1_ZaZb = new TH1D("h_eta_h1_ZaZb", "h_eta_h1_ZaZb", 16, -8., 8.);
   h_eta_h1_ZaZb->GetXaxis()->SetTitle("#eta");
   h_eta_h1_ZaZb->GetYaxis()->SetTitle("Number of Events");
      
   h_phi_h1_ZaZb = new TH1D("h_phi_h1_ZaZb", "h_phi_h1_ZaZb", 20, -10., 10.);
   h_phi_h1_ZaZb->GetXaxis()->SetTitle("#phi");
   h_phi_h1_ZaZb->GetYaxis()->SetTitle("Number of Events");
   
   // h1 -> b b~
   h_mh1_b1b2 = new TH1D("h_mh1_b1b2", "h_mh1_b1b2", 500, 0., 1000.);
   h_mh1_b1b2->GetXaxis()->SetTitle("Invariant mass of 2 b-jets (GeV/c^{2})");
   h_mh1_b1b2->GetYaxis()->SetTitle("Number of Events");
   
   h_pt_h1_b1b2 = new TH1D("h_pt_h1_b1b2", "h_pt_h1_b1b2", 500, 0., 500.);
   h_pt_h1_b1b2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_pt_h1_b1b2->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_h1_b1b2 = new TH1D("h_eta_h1_b1b2", "h_eta_h1_b1b2", 16, -8., 8.);
   h_eta_h1_b1b2->GetXaxis()->SetTitle("#eta");
   h_eta_h1_b1b2->GetYaxis()->SetTitle("Number of Events"); 
      
   h_phi_h1_b1b2 = new TH1D("h_phi_h1_b1b2", "h_phi_h1_b1b2", 20, -10., 10.);
   h_phi_h1_b1b2->GetXaxis()->SetTitle("#phi");
   h_phi_h1_b1b2->GetYaxis()->SetTitle("Number of Events"); 
   
   // h2 -> h1 h1 
   h_mh2_h1h1 = new TH1D("h_mh2_h1h1", "h_mh2_h1h1", 450, 100., 1000.);
   h_mh2_h1h1->GetXaxis()->SetTitle("Invariant mass of 2 b-jets + 4 muons (GeV/C^{2})");
   h_mh2_h1h1->GetYaxis()->SetTitle("Number of Events"); 
   
   h_pt_h2_h1h1 = new TH1D("h_pt_h2_h1h1", "h_pt_h2_h1h1", 700, 0., 700.);
   h_pt_h2_h1h1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   h_pt_h2_h1h1->GetYaxis()->SetTitle("Number of Events");
   
   h_eta_h2_h1h1 = new TH1D("h_eta_h2_h1h1", "h_eta_h2_h1h1", 16, -8., 8.);
   h_eta_h2_h1h1->GetXaxis()->SetTitle("#eta");
   h_eta_h2_h1h1->GetYaxis()->SetTitle("Number of Events");
   
   h_phi_h2_h1h1 = new TH1D("h_phi_h2_h1h1", "h_phi_h2_h1h1", 16, -8., 8.);
   h_phi_h2_h1h1->GetXaxis()->SetTitle("#phi");
   h_phi_h2_h1h1->GetYaxis()->SetTitle("Number of Events");
      
      
   // Declaring trees 
   TTree *muon1 = new TTree("muon1", "1st muon");
   TTree *muon2 = new TTree("muon2", "2nd muon");
   TTree *muon3 = new TTree("muon3", "3rd muon");
   TTree *muon4 = new TTree("muon4", "4th muon");
   
   
   //================================================================================================//
   //                                   DATASET SIGNAL H->hh->bb4Mu  NoPU                            //
   //================================================================================================//  
   //TFile *indata = new TFile("/media/aya/PACKUP/Aya/signal_20210831/BMP1_hh_bb_4Mu.root", "READ");  
   //TFile *indata = new TFile("/home/aya/programs/Delphes-3.4.2/sm_h_4Mu.root", "READ");
   //TFile *indata = new TFile("/home/aya/programs/Delphes-3.4.2/sm_h_bb.root", "READ");
   //TFile *indata = new TFile("/home/aya/Desktop/Pheno_Work/analysis/gg.root", "READ");
   //TFile *indata = new TFile("/HEP_DATA/aya/gghhbb4M_offshellsyntax_20210911.root", "READ");  
   //TFile *indata = new TFile("/media/aya/LinuxSpace/Pheno_Work/gghhbb4M_offshellsyntax_20210911.root", "READ");
   TFile *indata = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/analysis/gghhbb4M_AllDiagrams_20210915.root", "READ");
   
   
   //================================================================================================//
   //                                   DATASET Background 14 TeV                                    //
   //================================================================================================// 
   
   // DY 
   //TFile *input_file = new TFile("/media/aya/LinuxSpace/MyWork_Final_Samples/DY_BG_14TeV_SMfull_GEN/DY_BG_14TeV_SMfull_pythia8_CMSPhaseII-0PU_GEN-SIM.root", "READ");
    
   // ttbar
   //TFile *input_file = new TFile("/run/media/Aya/LinuxSpace/MyWork_Final_Samples/TTbar_BG_14TeV_SMfull_CMS_PhaseII_0PU_GEN_pythia8/TTbar_BG_14TeV_SMfull_pythia8_CMS_PhaseII_0PU_GEN-SIM.root", "READ");
   
   // ZZ_4Mu
   //TFile *input_file = new TFile("/media/aya/LinuxSpace/MyWork_Final_Samples/ZZto4Mu_BG_14TeV_SMfull_Pythia8_GEN/ZZto4MU_BG_14TeV_SMfull_pythia8_CMSPhaseII-0PU_GEN_SIM.root", "READ");
   
   
   //================================================================================================//
   //                                         Output Root files                                      //
   //================================================================================================//
   
   // DY 
   //TFile *op_file = new TFile("/home/aya/Desktop/Pheno_Work/analysis/Results/output_demo_DY.root", "RECREATE"); 
   //TFile *op_file = new TFile("/media/aya/PACKUP/Aya/output_demo_DY_2.root", "RECREATE"); 
   //TFile *op_file = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/analysis/output_demo_DY_20210916.root", "RECREATE"); 
   
   // ttbar
   //TFile *op_file = new TFile("/HEP_DATA/aya/Results/output_demo_ttbar.root", "RECREATE");
	
   //TFile *op_file = new TFile("/home/aya/Desktop/Pheno_Work/analysis/Results/output_demo_BMP1_hh_bb_4Mu.root", "RECREATE"); 
   //TFile *op_file = new TFile("/home/aya/Desktop/Pheno_Work/analysis/Results/output_sm_h_4Mu.root", "RECREATE"); 
   //TFile *op_file = new TFile("/home/aya/Desktop/Pheno_Work/analysis/Results/output_sm_h_bb.root", "RECREATE"); 
   //TFile *op_file = new TFile("/media/aya/PACKUP/Aya/output_gg_2.root", "RECREATE");
   //TFile *op_file = new TFile("/HEP_DATA/aya/Results/output_demo_gghhbb4M_20210911.root", "RECREATE"); 
   //TFile *op_file = new TFile("/home/aya/Desktop/Pheno_Work/analysis/Results/output_demo_gghhbb4M_20210914.root", "RECREATE"); 
   //TFile *op_file = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/output_demo_ZZ4M_20210914.root", "RECREATE"); 
   //TFile *op_file = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/analysis/output_demo_gghhbb4M_AllDiagrams_20210915.root", "RECREATE"); 
   TFile *op_file = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/analysis/output_demo_2_gghhbb4M_AllDiagrams_20210915.root", "RECREATE"); // adding control region cuts
   //TFile *op_file = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/analysis/output_demo_weighted_gghhbb4M_AllDiagrams_20210915.root", "RECREATE"); 
   
   
   //------------------------WEIGHT Calculation---------------------------
  
   float Lumi_data = 3.e+03;    // in 1/fb
	
   //-----------------------------------//
   //  Lumi_mc = nEvents/xsection(fb);  //
   //-----------------------------------//	
	
   //float Lumi_mc = 100000./42.32e-11;      // BMP1
   //float Lumi_mc = 10000./2.186e-04;       // sm_h_4Mu 
   //float Lumi_mc = 10000./ 220.2;          // sm_h_bb
   //float Lumi_mc = 10000./ 3.044e-05 ;     // gg_h_zz
   //float Lumi_mc = 5000./898225.;          // DY
   //float Lumi_mc = 1000000./7559.171362345;  // ttbar
   
   //float Lumi_mc = 1.e+05/ 3.62e-05; // gghhbb4M_offshellsyntax_alldiagrams
   //float wt = Lumi_data/Lumi_mc;
   float wt = 1.;    // examine plots with unweighted events
   
   /*===================================================================================*/  
  /*------------------------------Looping over ALL Events-----------------------------*/
  /*==================================================================================*/ 
   
   
   Long64_t nentries = fChain->GetEntriesFast();
   //Long64_t nentries = 10;
   //Long64_t nentries = 5000;
   
   //Long64_t nSelectedEvents = 0;
   
   //double Efficiency; 
    
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry = 0; jentry < nentries; jentry++) 
   {
      cout << "******START EVENT LOOP!******    ,    Event nb = " << jentry << endl; 
	  Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       
       // Loop overall Jets
       cout << "start loop overall jets" << endl;
       cout << "Jet size = " << Jet_size << endl;
       
       vector<Float_t> v_jet_pt;
       vector<Float_t> v_jet_eta;
       vector<UInt_t>  v_jet_bTag;
       v_jet_pt.clear();
       v_jet_eta.clear();
       v_jet_bTag.clear();
       
       for (Int_t i = 0; i < Jet_size; i++){
	 
	        Float_t jet_pt = Jet_PT[i];
			Float_t jet_eta = Jet_Eta[i];
			UInt_t jet_bTag = Jet_BTag[i];
			
			v_jet_pt.push_back(jet_pt);
            v_jet_eta.push_back(jet_eta);   
            v_jet_bTag.push_back(jet_bTag);      
	        
	        h_jets_size->Fill(i);
            h_pt_allJets->Fill(Jet_PT[i], wt);
	        h_eta_allJets->Fill(Jet_Eta[i], wt);
	        h_phi_allJets->Fill(Jet_Phi[i], wt);
	       
       }
       cout << "end loop overall jets" << endl;
       
       for ( int i = 0; i < v_jet_pt.size(); i++ ) { cout << "PT of Jets in GeV = " << v_jet_pt[i] << endl; }
       for ( int i = 0; i < v_jet_bTag.size(); i++){ cout << "BTag Discriminator = " << v_jet_bTag[i] << endl; }
       for ( int i = 0; i < v_jet_eta.size(); i++) { cout << "Eta of Jets = " << v_jet_eta[i] << endl; }
      
       // Loop over MET
       cout << "start loop overall MET" << endl;
       for (Int_t i = 0; i < MissingET_size; i++){
		  
	    h_MET_size->Fill(i);
            h_MET->Fill(MissingET_MET[i], wt);
            h_eta_MET->Fill(MissingET_Eta[i], wt);
            h_phi_MET->Fill(MissingET_Phi[i], wt);
	       
       }
       cout << "end loop overall MET" << endl;
      
       // Loop over all Muons Loose
       cout << "start loop overall loose muons" << endl;
       //cout << "Muon_Loose size = " << MuonLoose_size << endl;
       
       vector<Float_t> v_muon_pt;
       vector<Float_t> v_muon_eta;
       v_muon_pt.clear();
       v_muon_eta.clear();
       
       for (Int_t i = 0; i < MuonLoose_size; i++){
		   
		    Float_t muon_pt = MuonLoose_PT[i];
			Float_t muon_eta = MuonLoose_Eta[i];
			v_muon_pt.push_back(muon_pt);
			v_muon_eta.push_back(muon_eta);
		  
	        h_muons_size_Loose->Fill(i);
            h_pt_allMuons_Loose->Fill(MuonLoose_PT[i], wt);
	        h_eta_allMuons_Loose->Fill(MuonLoose_Eta[i], wt);
	        h_phi_allMuons_Loose->Fill(MuonLoose_Phi[i], wt);
		
       }  
       cout << "end loop overall Muon_loose" << endl;
       
       for ( int i = 0; i < v_muon_pt.size(); i++) { cout << "PT of Muons in GeV = " << v_muon_pt[i] << endl; }
       for ( int i = 0; i < v_muon_eta.size(); i++) { cout << "Eta of Muons = " << v_muon_eta[i] << endl; }
       
       
       cout << "start loop overall tight muons" << endl;
       for (Int_t i = 0; i < MuonTight_size; i++){
		  
	       h_muons_size_Tight->Fill(i);
           h_pt_allMuons_Tight->Fill(MuonTight_PT[i], wt);
	       h_eta_allMuons_Tight->Fill(MuonTight_Eta[i], wt);
	       h_phi_allMuons_Tight->Fill(MuonTight_Phi[i], wt);
		
       }  
       cout << "end loop overall Muon_tight" << endl;
	   
	   
       //------------------------GEN PARTICLES------------------------------ 	 
		 
	   // Get Thresholds for pT, eta for generated objects 
	 
	  // Looping overall Gen_Particles 	 
	  for (Int_t i = 0; i < Particle_size; i++){
		 
	        // Get particle pdg id 
	        int p_id = Particle_PID[i];
		 
	        // Check for Muons
	        if ( p_id == 13 ) { // pdg_id = 13 for Muon
			 
		         float gen_muon_pt = Particle_PT[i];
	             float gen_muon_eta = Particle_Eta[i];
		         float gen_muon_phi = Particle_Phi[i];

                 h_Gen_Muons_pt->Fill(gen_muon_pt, wt); 
                 h_Gen_Muons_eta->Fill(gen_muon_eta, wt); 
                 h_Gen_Muons_phi->Fill(gen_muon_phi, wt); 
              
                /*muon_pt->Fill(); 
		         muon_eta->Fill();
		         muon_phi->Fill(); */
			 
	        } // end if on id 13 	
	     
	        // Check for b quarks
	        if ( p_id == 5 ) { // pdg_id = 5 for b quark 
			 
		     float gen_b_pt = Particle_PT[i];
		     float gen_b_eta = Particle_Eta[i];
		     float gen_b_phi = Particle_Phi[i]; 
			 
		     h_Gen_Bjet_pt->Fill(gen_b_pt, wt); 
                     h_Gen_Bjet_eta->Fill(gen_b_eta, wt); 
                     h_Gen_Bjet_phi->Fill(gen_b_phi, wt);
			  
	         } // end if on id 5
		 
           } // end loop over gen particles    
           
           
	   
	   //=====================================================================//
       //                        Start 4Muons , ZZ Selections                 //
       //=====================================================================//
	      
	    vector<Int_t> v_muon_idx; // saves muon index if fullfill Object Selection 
	    v_muon_idx.clear();
	    cout << "Original Number of Muons: " << MuonLoose_size << endl;
	    
	   for ( Int_t i = 0; i < MuonLoose_size; i++ ){
			  
		    float muonL_pt = MuonLoose_PT[i];
		    float muonL_eta = fabs(MuonLoose_Eta[i]);
		     
			// 1st Selection: on pT, eta of reconstructed Muons (Object_Selection) 
			
			  if ( muonL_pt > 5. ){ 
				
				if ( muonL_eta < 2.4 ) {
					
					int mu_idx = i;
					v_muon_idx.push_back(mu_idx);
				  
			       
	            } //  end if muonL_eta < 2.4
            } // end if muonL_pt > 5
        } // end loop overall loose muons		       
         
      //  cout << "Number of Muons after Object Selection: " << v_muon_idx.size() << endl;
        
        // 2nd Selection: on number of Muons per event
        if ( v_muon_idx.size() > 3 ){ // having at least 4 Muons per event 
			  
	        // TLorentzVector declarations 
            TLorentzVector mu1, mu2, mu3, mu4, Z12, Z34, Z13, Z24, Z14, Z23, Za, Zb, b1, b2, h1, h2, H;
      
            // Set Pt, eta, phi mass for 4 muons TLV
            double muon_mass = 0.105658375;  // in GeV
            //double m_mu1, m_mu2, m_mu3, m_mu4;
            Float_t pt_mu1, pt_mu2, pt_mu3, pt_mu4;
            Float_t eta_mu1, eta_mu2, eta_mu3, eta_mu4;
            Float_t phi_mu1, phi_mu2, phi_mu3, phi_mu4;
      
            // Initialize variables
            pt_mu1 = -9999.; pt_mu2 = -9999.; pt_mu3 = -9999.; pt_mu4 = -9999.;
            eta_mu1 = -9999.; eta_mu2 = -9999.; eta_mu3 = -9999.; eta_mu4 = -9999.;
            phi_mu1 = -9999.; phi_mu2 = -9999.; phi_mu3 = -9999.; phi_mu4 = -9999.;
                
            int mu1_idx = v_muon_idx[0];
            int mu2_idx = v_muon_idx[1];
            int mu3_idx = v_muon_idx[2];
            int mu4_idx = v_muon_idx[3];
      
            // In Delphes tree Muon_PT branch pT are sorted from highest to least one
            mu1.SetPtEtaPhiM(MuonLoose_PT[mu1_idx], MuonLoose_Eta[mu1_idx], MuonLoose_Phi[mu1_idx], muon_mass); 
            mu2.SetPtEtaPhiM(MuonLoose_PT[mu2_idx], MuonLoose_Eta[mu2_idx], MuonLoose_Phi[mu2_idx], muon_mass);
            mu3.SetPtEtaPhiM(MuonLoose_PT[mu3_idx], MuonLoose_Eta[mu3_idx], MuonLoose_Phi[mu3_idx], muon_mass);
            mu4.SetPtEtaPhiM(MuonLoose_PT[mu4_idx], MuonLoose_Eta[mu4_idx], MuonLoose_Phi[mu4_idx], muon_mass);  
      
            // Leading Muon > 20 GeV (Muon with highest pT)
            pt_mu1 = mu1.Pt();
            eta_mu1 = mu1.Eta();
            phi_mu1 = mu1.Phi();
                
            // subleading Muon > 10 GeV (Muon with second-highest pT)    
            pt_mu2 = mu2.Pt();
            eta_mu2 = mu2.Eta();
            phi_mu2 = mu2.Phi();
                
            pt_mu3 = mu3.Pt();
            eta_mu3 = mu3.Eta();
            phi_mu3 = mu3.Phi();
                
            pt_mu4 = mu2.Pt();
            eta_mu4 = mu2.Eta();
            phi_mu4 = mu2.Phi();
                
                //TTree *muon1 = new TTree("muon1", "1st muon");
                /* muon1->Branch("pt", &pt_mu1, "pt_mu1/F");
                muon1->Branch("eta", &eta_mu1, "eta_mu1/F");
                muon1->Branch("phi", &phi_mu1, "phi_mu1/F");
              
                // TTree *muon2 = new TTree("muon2", "2nd muon");
                muon2->Branch("pt", &pt_mu2, "pt_mu2/F");
                muon2->Branch("eta", &eta_mu2, "eta_mu2/F");
                muon2->Branch("phi", &phi_mu2, "phi_mu2/F");
              
                //TTree *muon3 = new TTree("muon3", "3rd muon");
                muon3->Branch("pt", &pt_mu3, "pt_mu3/F");
                muon3->Branch("eta", &eta_mu3, "eta_mu3/F");
                muon3->Branch("phi", &phi_mu3, "phi_mu3/F");
              
                //TTree *muon4 = new TTree("muon4", "4th muon");
                muon4->Branch("pt", &pt_mu4, "pt_mu4/F");
                muon4->Branch("eta", &eta_mu4, "eta_mu4/F");
                muon4->Branch("phi", &phi_mu4, "phi_mu4/F");
              
                muon1->Fill();
                muon2->Fill();
                muon3->Fill();
                muon4->Fill(); */
              
              
            // Calculate DR bet each 2 muons for all possible combinations 1234, 1324, 1423
            double DR_mu12, DR_mu34, DR_mu13, DR_mu24, DR_mu14, DR_mu23; 
      
            // Initialize variables 
            DR_mu12 = -9999.; DR_mu34 = -9999.; DR_mu13 = -9999.; DR_mu24 = -9999.; DR_mu14 = -9999.; DR_mu23 = -9999.; 
      
            DR_mu12 = TMath::Sqrt(TMath::Power((eta_mu1 - eta_mu2), 2) + TMath::Power((phi_mu1 - phi_mu2), 2));
            DR_mu34 = TMath::Sqrt(TMath::Power((eta_mu3 - eta_mu4), 2) + TMath::Power((phi_mu3 - phi_mu4), 2));
            DR_mu13 = TMath::Sqrt(TMath::Power((eta_mu1 - eta_mu3), 2) + TMath::Power((phi_mu1 - phi_mu3), 2));
            DR_mu24 = TMath::Sqrt(TMath::Power((eta_mu2 - eta_mu4), 2) + TMath::Power((phi_mu2 - phi_mu4), 2));
            DR_mu14 = TMath::Sqrt(TMath::Power((eta_mu1 - eta_mu4), 2) + TMath::Power((phi_mu1 - phi_mu4), 2));
            DR_mu23 = TMath::Sqrt(TMath::Power((eta_mu2 - eta_mu3), 2) + TMath::Power((phi_mu2 - phi_mu3), 2));
      
            h_DR_mu1mu2->Fill(DR_mu12, wt);
            h_DR_mu3mu4->Fill(DR_mu34, wt);
            h_DR_mu1mu3->Fill(DR_mu13, wt);
            h_DR_mu2mu4->Fill(DR_mu24, wt);
            h_DR_mu1mu4->Fill(DR_mu14, wt);
            h_DR_mu2mu3->Fill(DR_mu23, wt); 
       
       
            /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
             ^                                         ^ 
             ^            Determine Za, Zb             ^
             ^                                         ^
             ^              for h -> Z Z               ^
             ^                                         ^ 
             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
      
      
            // 4 Muons various combination pairs 1234, 1324, 1423
            double mZ = 91.1876; // in GeV 
            double mZ12, mZ34, mZ13, mZ24, mZ14, mZ23, dZc1, dZc2, dZc3;
            double pt_Z12, pt_Z34, pt_Z13, pt_Z24, pt_Z14, pt_Z23;
            double eta_Z12, eta_Z34, eta_Z13, eta_Z24, eta_Z14, eta_Z23;
            double phi_Z12, phi_Z34, phi_Z13, phi_Z24, phi_Z14, phi_Z23;
            double dZ12, dZ23, dZ34, dZ13, dZ14, dZ24;
            double mZa, mZb, ptZa, ptZb, etaZa, etaZb, phiZa, phiZb;
      
            // Initialize variables
            mZ12 = -9999.; mZ34 = -9999.; mZ13 = -9999.; mZ24 = -9999.; mZ14 = -9999.; mZ23 = -9999.; dZc1 = -9999.; dZc2 = -9999.; dZc3 = -9999.;
            pt_Z12 = -9999.; pt_Z34 = -9999.; pt_Z13 = -9999.; pt_Z24 = -9999.; pt_Z14 = -9999.; pt_Z23 = -9999.;
            eta_Z12 = -9999.; eta_Z34 = -9999.; eta_Z13 = -9999.; eta_Z24 = -9999.; eta_Z14 = -9999.; eta_Z23 = -9999.;
            phi_Z12 = -9999.; phi_Z34 = -9999.; phi_Z13 = -9999.; phi_Z24 = -9999.; phi_Z14 = -9999.; phi_Z23 = -9999.;
            dZ12 = -9999.; dZ23 = -9999.; dZ34 = -9999.; dZ13 = -9999.; dZ14 = -9999.; dZ24 = -9999.;
            mZa = -9999.; mZb = -9999.; ptZa = -9999.; ptZb = -9999.; etaZa = -9999.; etaZb = -9999.; phiZa = -9999.; phiZb = -9999.;
      
      
            // 3rd Selection: on 4Muons charge 
            if ( MuonLoose_Charge[0] + MuonLoose_Charge[1] + MuonLoose_Charge[2] + MuonLoose_Charge[3] == 0){ //Sure that muon pairs are of opposite signs 

		       // Start 4 muon combination 1234 
		       if ( MuonLoose_Charge[0] + MuonLoose_Charge[1] == 0){  // mu1, mu2 
			  
                   if ( MuonLoose_Charge[2] + MuonLoose_Charge[3] == 0){  // mu3, mu4 
				      
				      // 4th Selection: on highest and second highest pT Muons
			          if ( ( pt_mu1 > 20. ) && ( pt_mu2 > 10. ) ){ 
						  
						  double m12 = (mu1 + mu2).M();
						  double m34 = (mu3 + mu4).M(); 
						  
						  // 5th Selection: on mass of each oppositely charged muon pairs should be mll > 4 GeV  
						  if ( m12 > 4. ) { 
						 
			                  Z12 = mu1 + mu2;
                              mZ12 = Z12.M();
                              pt_Z12 = Z12.Pt();
                              eta_Z12 = Z12.Eta();
                              phi_Z12 = Z12.Phi();
                              
				          } // end if m12
				          
				          if ( m34 > 4. ) {
				             
				              Z34 = mu3 + mu4;
                              mZ34 = Z34.M();
                              pt_Z34 = Z34.Pt();
                              eta_Z34 = Z34.Eta();
                              phi_Z34 = Z34.Phi();
			          
					      } // end if m34
			           
			              if (mZ12 > 0.) h_mZ12->Fill(mZ12, wt);
			              if (mZ34 > 0.) h_mZ34->Fill(mZ34, wt);
			         
			         } // end if on 4th selection 
		           } // end if on mu3,4 charge
		        } // end if on mu1,2 charge
		  
		        dZ12 = fabs(mZ12 - mZ);
		        dZ34 = fabs(mZ34 - mZ);
		    
		        // condition ? result_if_true : result_if_false  -> syntax for using ? conditional operator 
		        //dZc1 = ( dZ12 < dZ34 ) ? dZ12 : dZ34;
		          
		        if ( dZ12 < dZ34 ) {
				     
			        dZc1 = dZ12;
		            cout << "dZ mass for combination 1234 = dZ12 = " << dZc1 << endl;
	            }	
		    
		        else 
		        {
			        dZc1 = dZ34;
		            cout << "dZ mass for combination 1234 = dZ34 = " << dZc1 << endl;
		        }
		    
		        // Start 4 muon combination 1324 
		        if ( MuonLoose_Charge[0] + MuonLoose_Charge[2] == 0){  // mu1, mu3
				
		           if ( MuonLoose_Charge[1] + MuonLoose_Charge[3] == 0){ // mu2, mu4
					   
					   if ( ( pt_mu1 > 20. ) && ( pt_mu2 > 10. ) ) { 
						   
						   double m13 = (mu1 + mu3).M();
						   double m24 = (mu2 + mu4).M();
						   
						   if ( m13 > 4. ) { 
					
			                   Z13 = mu1 + mu3;
			                   mZ13 = Z13.M();
			                   pt_Z13 = Z13.Pt();
                               eta_Z13 = Z13.Eta();
                               phi_Z13 = Z13.Phi();
                               
						   } // end if m13
						   
						   if ( m24 > 4. ) { 
					
			                   Z24 = mu2 + mu4;
			                   mZ24 = Z24.M();
			                   pt_Z24 = Z24.Pt();
                               eta_Z24 = Z24.Eta();
                               phi_Z24 = Z24.Phi();
                               
						   } // end if m24
					
			               if (mZ13 > 0.) h_mZ13->Fill(mZ13, wt);
			               if (mZ24 > 0.) h_mZ24->Fill(mZ24, wt);
			              
					  } // end if on selection on two highest pT muons
			       } // end if on mu1,3 charge
	            } // end if on mu2,4 charge
		    
		        dZ13 = fabs(mZ13 - mZ);
		        dZ24 = fabs(mZ24 - mZ);
		
		        //dZc2 = ( dZ13 < dZ24 ) ? dZ13 : dZ24; 
		      
		        if ( dZ13 < dZ24 ) {
			
			        dZc2 = dZ13;
			        cout << "dZ mass for combination 1324 = dZ13 = " << dZc2 << endl;
		        }	
		          
		        else 
		        {
			        dZc2 = dZ24;
		            cout << "dZ mass for combination 1324 = dZ24 = " << dZc2 << endl;
		        }
		    
		        // Start 4 muon combination 1423 
		        if ( MuonLoose_Charge[0] + MuonLoose_Charge[3] == 0){  // mu1, mu4
				
		           if ( MuonLoose_Charge[1] + MuonLoose_Charge[2] == 0){ // mu2, mu3
					   
					   if ( ( pt_mu1 > 20. ) && ( pt_mu2 > 10. ) ) { 
						   
						   double m14 = (mu1 + mu4).M();
						   double m23 = (mu2 + mu3).M();
						   
						   if ( m14 > 4. ) {
							   
					           Z14 = mu1 + mu4;
			                   mZ14 = Z14.M();
			                   pt_Z14 = Z14.Pt();
			                   eta_Z14 = Z14.Eta();
			                   phi_Z14 = Z14.Phi();
			                   
						   } // end if m14
						   
						   if ( m23 > 4. ) {
					
			                   Z23 = mu2 + mu3;
			                   mZ23 = Z23.M();
			                   pt_Z23 = Z23.Pt();
		                       eta_Z23 = Z23.Eta();
			                   phi_Z23 = Z23.Phi();
			                   
						   } // end if m23
					
			               if (mZ14 > 0.) h_mZ14->Fill(mZ14, wt);
			               if (mZ23 > 0.) h_mZ23->Fill(mZ23, wt);
			              
					   } // end if on selection on two highest pT muons
			        } // end if on mu1,4 charge
		         } // end if on mu2,3 charge
		    
		         dZ14 = fabs(mZ14 - mZ);
		         dZ23 = fabs(mZ23 - mZ);
		
		         //dZc3 = ( dZ14 < dZ23 ) ? dZ14 : dZ23; 
		      
		         if ( dZ14 < dZ23 ) {
				      
			         dZc3 = dZ14;
		             cout << "dZ mass for combination 1423 = dZ14 = " << dZc3 << endl;
		         }	
		        
		         else 
		         {
			         dZc3 = dZ23;
		             cout << "dZ mass for combination 1423 = dZ23 = " << dZc3 << endl;
		         } 
       
       
                 // Choose dZc of the least value bet dZc1, dZc2, dZc3 to be Za, Zb combination
		      
		         if ( dZc1 < dZc2 && dZc1 < dZc3 ){  // dZc1 < dZc2, dZc3
			    
			        if ( dZ12 < dZ34 ){
					
			            Za = Z12;
			            mZa = mZ12;
                        ptZa = pt_Z12;
                        etaZa = eta_Z12;
                        phiZa = phi_Z12;
					
				        Zb = Z34;
			            mZb = mZ34;
                        ptZb = pt_Z34;
                        etaZb = eta_Z34;
			            phiZb = phi_Z34;
					
			        } 
				
			        else
			        { 
			            Za = Z34;
			            mZa = mZ34;
                        ptZa = pt_Z34;
                        etaZa = eta_Z34;
                        phiZa = phi_Z34;
                    
			            Zb = Z12;
			            mZb = mZ12;
                        ptZb = pt_Z12;
                        etaZb = eta_Z12;
                        phiZb = phi_Z12;
					
			         } 
			
		         } // end if dZc1   
			
		         else if ( dZc2 < dZc1 && dZc2 < dZc3 ){  // dZc2 < dZc1, dZc3
					 
				     if ( dZ13 < dZ24 ){
					
				         Za = Z13;
				         mZa = mZ13;
                         ptZa = pt_Z13;
                         etaZa = eta_Z13;
                         phiZa = phi_Z13;
					
				         Zb = Z24;
			             mZb = mZ24;
                         ptZb = pt_Z24;
                         etaZb = eta_Z24;
                         phiZb = phi_Z24;
                    
				     }
				 
				     else
				     {      
				         Za = Z24;
				         mZa = mZ24;
                         ptZa = pt_Z24;
                         etaZa = eta_Z24;
                         phiZa = phi_Z24;
					
				         Zb = Z13;
				         mZb = mZ13;
                         ptZb = pt_Z13;
                         etaZb = eta_Z13;
                         phiZb = phi_Z13;
					
				     }
				
			     } // end else if dZc2
		    
		         else 
			     {  // dZc3 < dZc1, dZc2
					
			        if ( dZ14 < dZ23 ){
					
			             Za = Z14;
			             mZa = mZ14;
                         ptZa = pt_Z14;
                         etaZa = eta_Z14;
                         phiZa = phi_Z14;
					
			             Zb = Z23;
			             mZb = mZ23;
                         ptZb = pt_Z23;
                         etaZb = eta_Z23;
                         phiZb = phi_Z23;
					
		            }
				
                    else
		            {
		                 Za = Z23;
			             mZa = mZ23;
                         ptZa = pt_Z23;
                         etaZa = eta_Z23;
                         phiZa = phi_Z23;
					
			             Zb = Z14;
			             mZb = mZ14;
                         ptZb = pt_Z14;
                         etaZb = eta_Z14;
                         phiZb = phi_Z14;
					
		            }
				
			      } // end else dZc3
			
                  
                  // 6th Selection: on Za and Zb masses, Za > 40 GeV (closest to nomial Z mass) & Zb > 12 GeV
                  
                  if ( mZa > 40. && mZa < 120.){
					  
	                 if ( mZb > 12. && mZb < 120. ){
						 
					     h_mZa_4mu->Fill(mZa, wt);            
                         h_mZb_4mu->Fill(mZb, wt);
                         h_pt_Za->Fill(ptZa, wt);
                         h_pt_Zb->Fill(ptZb, wt);
                         h_eta_Za->Fill(etaZa, wt);
                         h_eta_Zb->Fill(etaZb, wt);
		                 h_phi_Za->Fill(phiZa, wt);
		                 h_phi_Zb->Fill(phiZb, wt); 
		            
		            
		                 //============================
		                 // Reconstruct h1 from Za, Zb
		                 //============================
            
                         double mh1_ZaZb, pt_h1_ZaZb, eta_h1_ZaZb, phi_h1_ZaZb;
                   
                         // Initialize Variables
                         mh1_ZaZb = -9999.; pt_h1_ZaZb = -9999.; eta_h1_ZaZb = -9999.; phi_h1_ZaZb = -9999.;
                   
                         h1 = Za + Zb;
                         mh1_ZaZb = h1.M(); // get invariant mass of 4Muons 
                         
                         // 7th Selection: on invariant mass of 4Muons, in Signal Region 115 ≤ m4l ≤ 135 GeV  
                         if ( ( mh1_ZaZb >= 115. ) && ( mh1_ZaZb <= 135. ) ) {
                        // if ( ( mh1_ZaZb >= 115. ) && ( mh1_ZaZb <= 135. ) ) continue;
                         
                         // 8th Selection: on control region or side bands  m4l < 115 or m4l > 135 GeV
						// if ( ( mh1_ZaZb < 115. ) || ( mh1_ZaZb > 135. ) ) continue;	
						  
                             pt_h1_ZaZb = h1.Pt();
                             eta_h1_ZaZb = h1.Eta();
                             phi_h1_ZaZb = h1.Phi();
            
                             h_mh1_ZaZb->Fill(mh1_ZaZb, wt);
                             h_pt_h1_ZaZb->Fill(pt_h1_ZaZb, wt);
                             h_eta_h1_ZaZb->Fill(eta_h1_ZaZb, wt);
                             h_phi_h1_ZaZb->Fill(phi_h1_ZaZb, wt);
		            
					         
					         /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                              ^                                         ^ 
                              ^            Determine b1, b2             ^
                              ^                                         ^
                              ^              for h -> b1 b2             ^
                              ^                                         ^ 
                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
       
                            bool found_bjet = false;
                            int nbjets = 0;            // total nb of b-jets found per event
                            vector<Int_t> bjet_indx;
                            bjet_indx.clear();
            
                            // Loop overall jets
		                    cout << "Jet_Size =  " << Jet_size << " for Event nb  " <<  jentry << endl;
                            
                            for ( Int_t i = 0; i < Jet_size; i++ ){
					   
                                 UInt_t jet_bTag = Jet_BTag[i];  // save BTag value for each jet
                                  
				                 
                                 // 1st Selection: select loose B-jets WP with BTag values >= 4 (4, 5, 6, 7) stored at bit 0 in DelphesCard
                                 if ( ( jet_bTag == 4 ) || ( jet_bTag > 4 ) ) {
									 
								 //	1st Selection: select Tight B-jets WP with BTag values 1, 3, 5, 7  stored at bit 2 in DelphesCard
								 // if ( ( jet_bTag == 1 ) || ( jet_bTag == 3 ) || ( jet_bTag == 5 ) || ( jet_bTag == 7 ) ) {
									 
									 found_bjet = true;
									 
									 if ( found_bjet ){
										 
										 nbjets++;
										 
										 int bj_indx = i;
					
					                     cout << "Jet [" << bj_indx
                                              << "]  has BTag discriminator = " << Jet_BTag[i]
                                              <<  " and Jet_PT = " << Jet_PT[i] <<  " GeV/C" << endl; 
									 
									     bjet_indx.push_back(bj_indx);  // filling is sorted with order of highest pT at first element then decreases 
									 
									 } // end if found_bjet
								 } // end if jet_bTagS
							 } // end loop overall jets
							 
							 h_bjet_size->Fill(nbjets, wt);
							 
							 	
							 // 2nd Selection: on having at least 2 b-jets/event 
							 if ( nbjets > 1 ){ 
								 
								 cout << "========================================" << endl;
		 	                     cout << " number of b-jets/event =  " << nbjets    << endl;
		                         cout << "========================================" << endl;
							 		 
								 // save index of b-jets after applying series of selections
								 vector<Int_t> bjet_indx_AfterSel;
								 bjet_indx_AfterSel.clear();
								 
								 for ( Int_t i = 0; i < bjet_indx.size(); i++ ) { cout << "bjet index for v_element [" << i << "] is : "<< bjet_indx[i] << endl;}
							    
							     // Loop overall b-jets/event 
								 for ( Int_t i = 0; i < bjet_indx.size(); i++ ) { // loop over bjet_indx vector elements	
									 
									 int bj_indx_AfterSelec = bjet_indx[i];
									 double bjet_pt = Jet_PT[bj_indx_AfterSelec];
									 double bjet_eta = Jet_Eta[bj_indx_AfterSelec];
									 double abs_bjet_eta = fabs(Jet_Eta[bj_indx_AfterSelec]);
									 double bjet_phi = Jet_Phi[bj_indx_AfterSelec]; 
									 
					
									 double DEta_b_mu1_sqr = TMath::Power((bjet_eta - eta_mu1), 2);
									 double DEta_b_mu2_sqr = TMath::Power((bjet_eta - eta_mu2), 2);
									 double DEta_b_mu3_sqr = TMath::Power((bjet_eta - eta_mu3), 2);
									 double DEta_b_mu4_sqr = TMath::Power((bjet_eta - eta_mu4), 2);
									 double DPhi_b_mu1_sqr = TMath::Power((bjet_phi - phi_mu1), 2);
									 double DPhi_b_mu2_sqr = TMath::Power((bjet_phi - phi_mu2), 2);
									 double DPhi_b_mu3_sqr = TMath::Power((bjet_phi - phi_mu3), 2);
									 double DPhi_b_mu4_sqr = TMath::Power((bjet_phi - phi_mu4), 2);
									 
									 
									 double DR_b_mu1 = TMath::Sqrt( DEta_b_mu1_sqr + DPhi_b_mu1_sqr );
									 double DR_b_mu2 = TMath::Sqrt( DEta_b_mu2_sqr + DPhi_b_mu2_sqr );
									 double DR_b_mu3 = TMath::Sqrt( DEta_b_mu3_sqr + DPhi_b_mu3_sqr );
									 double DR_b_mu4 = TMath::Sqrt( DEta_b_mu4_sqr + DPhi_b_mu4_sqr );
                
									 
									 // 3rd Selection: on DeltaR(b-jet,lepton) of ZZ candidates > 0.3
									 if ( ( DR_b_mu1 > 0.3) && ( DR_b_mu2 > 0.3) && ( DR_b_mu3 > 0.3) && ( DR_b_mu4 > 0.3) ) {
										 
										 cout << "......Selection on DR between (bjet, lepton) DONE......"  << endl;
										 // 4th Selection: on b-jet pT > 20 GeV
										 if ( bjet_pt > 20.) {   
											 
											 cout << "......Selection on bjet PT DONE......"  << endl;
											 // 5th Selection on b-jet abs_eta < 2.4
											 if ( abs_bjet_eta < 2.4) {
												 
											     // save b-jet index for further selection 
											     bjet_indx_AfterSel.push_back(bj_indx_AfterSelec);
											     cout << "......Selection on bjet abs Eta DONE......"  << endl;
											} 
											
											 
										 } // end if bjet_pt 	 
								      } // end if DR 
				                   } // end loop over vector elements 
                                   
                                
                                
                                 // store BTag scores for b-jets survived the above selections 
					             vector<Int_t> v_BTag_scores; 
								 v_BTag_scores. clear(); 
								   
								 cout << "......Start Final Selection on BTagScore......" << endl;
								   
					             // make sure that we still have at least 2 bjets after above selections
					             if ( bjet_indx_AfterSel.size() > 1 ) {  
									    
									    cout << "bjet_indx_AfterSel size is : " <<  bjet_indx_AfterSel.size() << endl; 
                                        cout << "bjet indx = " << bjet_indx_AfterSel[0] << endl;
									    
									    // save BTag scores for selected bjets
									    for ( Int_t i = 0; i < bjet_indx_AfterSel.size(); i++ ) {  
											
											int bjet_indx_final = bjet_indx_AfterSel[i];
										    Int_t BTag_score = Jet_BTag[bjet_indx_final];
										   
										    cout << "B-jet [" << bjet_indx_final 
										         << "] has BTag_Score = " << BTag_score << endl;
										         
										    v_BTag_scores.push_back(BTag_score);
					              
								        } //end loop on bjet_indx_AfterSel.size()
								        
								        for ( Int_t i = 0; i < v_BTag_scores.size(); i++ ) { cout << " BTag_Scores = " << v_BTag_scores[i] << " "; }
		                         
		                
		                                // sort elements of BTagScores vector in descending order (starting with highest score)
		                                sort(v_BTag_scores.begin(), v_BTag_scores.end(), greater<int>());
		                         
		                                for ( Int_t i = 0; i < v_BTag_scores.size(); i++ ) { cout << "sorted BTag_Scores = " << v_BTag_scores[i] << endl; }  
		                         
		                                
		                                // 5th Selection: on BTag Scores of b-jets (select 2 b-jets with highest BTag score)
		                                int max_BTag_score_b1 = v_BTag_scores[0];  // highest bTag score for 1st b-jet
		                                int max2_BTag_score_b2 = v_BTag_scores[1];  // second highest bTag score for 2nd b-jet
				    
				    
				                        // Get index of 2 b-jets with highest scores
				                        Int_t btag, b_id, signal_bjet_1_indx, signal_bjet_2_indx;
				                 
									 
				                        // save signal 2 bjets index with highest BTag score   
				                        for ( Int_t i = 0; i < bjet_indx_AfterSel.size(); i++ ) { 
											
											b_id =  bjet_indx_AfterSel[i];
									        btag =  Jet_BTag[b_id];
									        
									        if ( btag == max_BTag_score_b1 ) signal_bjet_1_indx = b_id; 
				                          
				                            else if ( btag == max2_BTag_score_b2 ) signal_bjet_2_indx = b_id;
				                          
				                            else { cout << "btag for bjet [" << b_id << "] doesn't match any of highest scores" << endl;}
				                     
							            }
							     
							            cout << "2 B-jets of signal are ( " << signal_bjet_1_indx << ", " << signal_bjet_2_indx 
							                 << " ) with Highest BTag Scores " << max_BTag_score_b1 << ", " << max2_BTag_score_b2
							                 << " respectively!" << endl;
							          
							          
							            // Set TLorentzVectors for selected 2 b-jets of signal
							            b1.SetPtEtaPhiM( Jet_PT[signal_bjet_1_indx], 
							                             Jet_Eta[signal_bjet_1_indx], 
							                             Jet_Phi[signal_bjet_1_indx], 
							                             Jet_Mass[signal_bjet_1_indx] );
							                     
							            b2.SetPtEtaPhiM( Jet_PT[signal_bjet_2_indx], 
							                             Jet_Eta[signal_bjet_2_indx], 
							                             Jet_Phi[signal_bjet_2_indx], 
							                             Jet_Mass[signal_bjet_2_indx] );
							     
							     
							            double b1_mass =  b1.M();
							            double b1_pt   =  b1.Pt();
							            double b1_eta  =  b1.Eta();
							            double b1_phi  =  b1.Phi();
							     
							            double b2_mass =  b1.M();
							            double b2_pt   =  b1.Pt();
							            double b2_eta  =  b1.Eta();
							            double b2_phi  =  b1.Phi();
							     
							            double DeltaEta_b1b2_sqr = TMath::Power( ( b1_eta - b2_eta ), 2); 
							            double DeltaPhi_b1b2_sqr = TMath::Power( ( b1_phi - b2_phi ), 2); 
							            double DR_b1b2 = TMath::Sqrt( DeltaEta_b1b2_sqr + DeltaPhi_b1b2_sqr );
							     
							            h_mb_jet_1->Fill(b1_mass, wt);
                                        h_mb_jet_2->Fill(b2_mass, wt);
							            h_pt_b_jet_1->Fill(b1_pt, wt);
                                        h_pt_b_jet_2->Fill(b2_pt, wt);
                                        h_eta_b_jet_1->Fill(b1_eta, wt);
                                        h_eta_b_jet_2->Fill(b2_eta, wt);
                                        h_phi_b_jet_1->Fill(b1_phi, wt);
                                        h_phi_b_jet_2->Fill(b2_phi, wt);
                                        //h_DR_b1b2->Fill(DR_b1b2, wt); 
                                 
                                 
							            //========================================//
		                                //       Reconstruct h2 from b1, b2       //
		                                //========================================//
							     
							            // TLorentzVector for 2nd SM higgs
							            h2 = b1 + b2;
							     
							            double h2_b1b2_mass =  h2.M();
							            double h2_b1b2_pt   =  h2.Pt();
							            double h2_b1b2_eta  =  h2.Eta();
							            double h2_b1b2_phi  =  h2.Phi();
							     
							            h_mh1_b1b2->Fill(h2_b1b2_mass, wt);
                                        h_pt_h1_b1b2->Fill(h2_b1b2_pt, wt);
                                        h_eta_h1_b1b2->Fill(h2_b1b2_eta, wt);
                                        h_phi_h1_b1b2->Fill(h2_b1b2_phi, wt);
							     
							     
                                        //======================================//
		                                //                                       //
		                                //      Reconstruct BSM H from SM h      // 
		                                //               H -> h h                //
		                                //                                       //
		                                //=======================================//
                
                                        // TLorentzVector for BSM Heavy Higgs 
                                        H = h1 + h2;
                                 
                                        double H_bb4Mu_mass =  H.M();
                                        double H_bb4Mu_pt   =  H.Pt();
                                        double H_bb4Mu_eta  =  H.Eta();
                                        double H_bb4Mu_phi  =  H.Phi(); 
                                 
                                        h_mh2_h1h1->Fill(H_bb4Mu_mass, wt);
                                        h_pt_h2_h1h1->Fill(H_bb4Mu_pt, wt);
                                        h_eta_h2_h1h1->Fill(H_bb4Mu_eta, wt);
                                        h_phi_h2_h1h1->Fill(H_bb4Mu_phi, wt); 
								       
		                         } // end if bjet_indx_AfterSel.size() > 1 
				             } // end if nbjets > 1      
		                 } // end if on m4l invariant mass selection
		            } // end if mZb  
                  } // end if mZa
                } // end if on 4 muons charge 
              } // end if MuonLoose_size > 3
		     		  
 
   } // end loop overall events
   
  // Efficiency = nSelectedEvents/nentries;
   
   cout << "***Analysis Loop Ends!***" << endl; 
   
   
  /* cout << " Total Number of Events = " << nentries 
        << ", Number of Selected Events = " << nSelectedEvents
        << ", Efficiency = Nsel/Ngen = " << Efficiency << endl;
  */    
   
   cout << "Writing Histograms!" << endl;  
   
   h_Gen_Muons_pt->Write();
   h_Gen_Muons_eta->Write();
   h_Gen_Muons_phi->Write();
   h_Gen_Bjet_pt->Write();
   h_Gen_Bjet_eta->Write();
   h_Gen_Bjet_phi->Write();
   h_muons_size_Loose->Write();
   h_pt_allMuons_Loose->Write();
   h_eta_allMuons_Loose->Write();
   h_phi_allMuons_Loose->Write();
   h_muons_size_Tight->Write();
   h_pt_allMuons_Tight->Write();
   h_eta_allMuons_Tight->Write();
   h_phi_allMuons_Tight->Write();
   h_jets_size->Write();
   h_pt_allJets->Write();
   h_eta_allJets->Write();
   h_phi_allJets->Write();
   h_pt_all_B_Jets->Write();
   h_eta_all_B_Jets->Write();
   h_phi_all_B_Jets->Write();
   h_bjet_size->Write();	
   h_MET_size->Write();
   h_MET->Write();
   h_eta_MET->Write();
   h_eta_MET->Write();  
   h_DR_mu1mu2->Write();
   h_DR_mu3mu4->Write();
   h_DR_mu1mu3->Write();
   h_DR_mu2mu4->Write();
   h_DR_mu1mu4->Write();
   h_DR_mu2mu3->Write();
   h_mZ12->Write();
   h_mZ34->Write();
   h_mZ13->Write();
   h_mZ24->Write();
   h_mZ14->Write();
   h_mZ23->Write();
   h_mZa_4mu->Write();            
   h_mZb_4mu->Write();
   h_pt_Za->Write();
   h_pt_Zb->Write();
   h_eta_Za->Write();
   h_eta_Zb->Write();
   h_phi_Za->Write();
   h_phi_Zb->Write();
   h_mh1_ZaZb->Write();
   h_pt_h1_ZaZb->Write();
   h_eta_h1_ZaZb->Write();
   h_phi_h1_ZaZb->Write();
   h_mb_jet_1->Write();
   h_mb_jet_2->Write();
   h_pt_b_jet_1->Write();
   h_pt_b_jet_2->Write();
   h_eta_b_jet_1->Write();
   h_eta_b_jet_2->Write();
   h_phi_b_jet_1->Write();
   h_phi_b_jet_2->Write();
  //h_DR_b1b2->Write(); 
   h_mh1_b1b2->Write();
   h_pt_h1_b1b2->Write();
   h_eta_h1_b1b2->Write();
   h_phi_h1_b1b2->Write();
   h_mh2_h1h1->Write();
   h_pt_h2_h1h1->Write();
   h_eta_h2_h1h1->Write();
   h_phi_h2_h1h1->Write();
   
   
   /*muon1->Write();
   muon2->Write();
   muon3->Write();
   muon4->Write(); */
   
   //std::cout << "Writing Histograms ends!" << endl;
   
   op_file->Write();
   std::cout << "write output file! " << endl;
   
   op_file->Close();   
   
   std::cout << "saving..." << endl;    
   std::cout << "ROOT file: " << op_file << " has been created sucessfully!" << endl;
   
   
}
