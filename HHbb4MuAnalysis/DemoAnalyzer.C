////////////////////////////////////////////////////////
//      To run macro in a ROOT session do:            //
//     gSystem->Load("/path/to/libDelphes.so");       //
//     .L path/to/DemoAnalyzer.C                      //
//     DemoAnalyzer e;                                //
//     e.Loop();                                      //
////////////////////////////////////////////////////////


#define DemoAnalyzer_cxx

R__LOAD_LIBRARY(libPhysics)
//R__LOAD_LIBRARY(/home/aya/programs/Delphes-3.4.2/libDelphes.so)
R__ADD_LIBRARY_PATH(/home/aya/programs/Delphes-3.4.2)
R__LOAD_LIBRARY(libDelphes)

//#include "/HEP_DATA/aya/DemoAnalyzer.h"
#include "/home/aya/DemoAnalyzer.h"
#include "/home/aya/newtree.h"
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
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TFormula.h"
#include <TMath.h>
#include "THStack.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>



void DemoAnalyzer::Loop()
{

   if (fChain == 0) return;
   
   
   // Array act as a counter stores how many events passed a certain selection, 
   Int_t NEvents[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    // size = nb of selection cuts, each element referes to certain selection  
   
   
   //================================================================================================//
   //                                   DATASET SIGNAL H->hh->bb4Mu  NoPU                            //
   //================================================================================================//  
   
   TFile *in_sig = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/gghhbb4Mu_final.root", "READ");
   
   
   //================================================================================================//
   //                                   DATASET Background 14 TeV                                    //
   //================================================================================================// 
   
   //TFile* in_ttbar    = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ttbar.root"    , "READ");
   //TFile *in_ZZ4Mu    = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ZZ4Mu.root"    , "READ");
   //TFile *in_DY       = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/DY.root"       , "READ");   // haven't produced yet
   //TFile *in_ZZbb2Mu  = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ZZbb2Mu.root"  , "READ");
   //TFile *in_ZWpbbMuNL= new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ZWpbbMuNL.root", "READ");
   //TFile *in_ZZZbb4Mu = new TFile("/media/aya/PACKUP/Aya/final_analysis/inputs/ZZZbb4Mu.root" , "READ");   // haven't produced yet
   
   //================================================================================================//
   //                                         Output Root files                                      //
   //================================================================================================//
   
   TFile* out_sig    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ggHhhbb4M_final.root"   , "RECREATE");
   //TFile* ttbar    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ttbar.root"             , "RECREATE");
   //TFile* out_ZZ4Mu= new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZZ4Mu.root"             , "RECREATE");
 /*TFile* DY         = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_DY.root"                , "RECREATE");
   TFile* ZZbb2Mu    = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZZbb2Mu.root"           , "RECREATE");
   TFile* ZZZbb4Mu   = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZZZbb4Mu.root"          , "RECREATE");
   TFile* ZZ4b       = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZZ4b.root"              , "RECREATE");
   TFile* ZWpbbMuNL  = new TFile("/media/aya/PACKUP/Aya/final_analysis/outputs/output_demo_ZWpbbMuNL.root"         , "RECREATE");
   TFile* ZWp2MuMuNL = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_ZWpbbMuNL.root"    , "RECREATE");
   TFile* ZWmbbMuNL  = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_ZWmbbMuNL.root"    , "RECREATE");
   TFile* ZWm2MuMuNL = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_ZWmbbMuNL.root"    , "RECREATE"); 
   TFile* WW2Mu2NL   = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_WW2Mu2NL.root"     , "RECREATE"); 
   TFile* smHZZ4Mu   = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_smHZZ4Mu.root"     , "RECREATE"); 
   TFile* smHbb      = new TFile("/media/aya/LinuxSpace/Pheno_Work_2/final_analysis/output_demo_smHbb.root"        , "RECREATE");  */
   
   
   
   //------------------------WEIGHT Calculation---------------------------
  
   float Lumi_data = 3.e+03;    // in 1/fb
	
   //-----------------------------------//
   //  Lumi_mc = nEvents/xsection(fb);  //
   //-----------------------------------//	
	
   //float Lumi_mc = 1.e+06/898225.;          // DY
   //float Lumi_mc = 1.e+06/7559.17;          // ttbar
   //float Lumi_mc = 1.e+06/ 12.16;           // ZZ4M  
   float Lumi_mc = 1.e+05/ 3.62e-05;          // gghhbb4M_offshellsyntax_alldiagrams
   //float Lumi_mc = 1.e+06/106.6518;         // ZZbb2Mu
   //float Lumi_mc = 1.e+06/280.1;            // ZWpbbMuNL
   float wt = Lumi_data/Lumi_mc;
   //float wt = 1.;    // examine plots with unweighted events
   
   
   //------------------------END WEIGHT CALC----------------------------
   
   
  
   //---------------------START Defining Tree Branches--------------------
   
  
   b_Gen_Muon           = new_tree->Branch("Gen_Muons"           ,   &genMuons  );
   b_Gen_Bjets          = new_tree->Branch("Gen_BJets"           ,   &gen_bjet  );
   b_LooseMuons         = new_tree->Branch("LooseMuons"          ,   &loose_mu  );
   b_Jets               = new_tree->Branch("Jets"                ,   &jets      );
   b_MET                = new_tree->Branch("MET"                 ,   &met       );
   b_4Muons_beforeCut   = new_tree->Branch("Four_Muons_beforeCut",   &four_muons_beforeCut);
   b_deltaR_muons       = new_tree->Branch("DeltaR_Muons"        ,   &drMuons   );
   b_Zmass_comb         = new_tree->Branch("Zmass_Combinations"  ,   &zmass     );
   b_Za_signal          = new_tree->Branch("ZaOFSignal"          ,   &za_sig    );
   b_Zb_signal          = new_tree->Branch("ZbOFSignal"          ,   &zb_sig    );
   b_b1_signal          = new_tree->Branch("bjet_1_OFSignal"     ,   &b1_sig    );
   b_b2_signal          = new_tree->Branch("bjet_2_OFSignal"     ,   &b2_sig    );
   b_1st_Higgs_Ofsignal = new_tree->Branch("1st_Higgs_OFSignal"  ,   &smhiggs1  );
   b_2nd_Higgs_Ofsignal = new_tree->Branch("2nd_Higgs_OFSignal"  ,   &smhiggs2  );
   b_BSM_Higgs_Ofsignal = new_tree->Branch("Heavy_Higgs_OFSignal",   &heavyHiggs); 
  
   
   //---------------------END Defining Tree Branches--------------------
   
   
   /*===================================================================================*/  
  /*------------------------------Looping over ALL Events-----------------------------*/
  /*==================================================================================*/ 
   
   
   Long64_t nentries = fChain->GetEntriesFast();
   cout << "nEvents = " << nentries << endl;

   // define some vectors 

   vector<Int_t> v_muon_idx;            // saves muon index if fullfill Object Selection 
   vector<Int_t> bjet_indx ;            // saves tight b-jet index 
   vector<Int_t> bjet_indx_AfterSel; 	// saves index of tight b-jets after applying series of event selections on them 
								 
   //Long64_t nentries = 10;
   //Long64_t nentries = 1000;
   
   //Long64_t nSelectedEvents = 0;
   
   //double Efficiency; 
    
   Long64_t nbytes = 0, nb = 0;
   
   for (Long64_t jentry = 0; jentry < nentries; jentry++)  
   {
      cout << "******START EVENT LOOP!******    ,    Event nb = " << jentry << endl; 
	  
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      cout << "nbytes for current entry = " << nbytes << endl;
      // if (Cut(ientry) < 0) continue;
       
     
      //============Initialize vectors=============
     
      genMuons.gen_muon_pt.clear();
      genMuons.gen_muon_eta.clear();
      genMuons.gen_muon_phi.clear();
      gen_bjet.gen_bjet_pt.clear();
      gen_bjet.gen_bjet_eta.clear();
      gen_bjet.gen_bjet_phi.clear();
      loose_mu.lmu_pt.clear();
      loose_mu.lmu_eta.clear();
      loose_mu.lmu_phi.clear();
      jets.jet_pt.clear();
      jets.jet_eta.clear();
      jets.jet_phi.clear();
      met.met_MET.clear();
      met.met_eta.clear();
      met.met_phi.clear();
      four_muons_beforeCut.muon1_pt.clear();
      four_muons_beforeCut.muon2_pt.clear();
      four_muons_beforeCut.muon3_pt.clear();
      four_muons_beforeCut.muon4_pt.clear();
      four_muons_beforeCut.muon1_eta.clear();
      four_muons_beforeCut.muon2_eta.clear();
      four_muons_beforeCut.muon3_eta.clear();
      four_muons_beforeCut.muon4_eta.clear();
      four_muons_beforeCut.muon1_phi.clear();
      four_muons_beforeCut.muon2_phi.clear();
      four_muons_beforeCut.muon3_phi.clear();
      four_muons_beforeCut.muon4_phi.clear();
      drMuons.v_DR_mu1mu2.clear();
      drMuons.v_DR_mu3mu4.clear();
      drMuons.v_DR_mu1mu3.clear();
      drMuons.v_DR_mu2mu4.clear();
      drMuons.v_DR_mu1mu4.clear();
      drMuons.v_DR_mu2mu3.clear();
      zmass.Z12_mass.clear();
      zmass.Z34_mass.clear();
      zmass.Z13_mass.clear();
      zmass.Z24_mass.clear();
      zmass.Z14_mass.clear();
      zmass.Z23_mass.clear();
      za_sig.Za_mass.clear();
      za_sig.Za_pt.clear();
      za_sig.Za_eta.clear(); 
      zb_sig.Zb_mass.clear();
      zb_sig.Zb_pt.clear();
      zb_sig.Zb_eta.clear();
      b1_sig.b1_mass.clear();
      b1_sig.b1_pt.clear();
      b1_sig.b1_eta.clear();
      b2_sig.b2_mass.clear();
      b2_sig.b2_pt.clear();
      b2_sig.b2_eta.clear();
      smhiggs1.v_h1_mass.clear();
      smhiggs1.v_h1_pt.clear();
      smhiggs1.v_h1_eta.clear();
      smhiggs2.v_h2_mass.clear();
      smhiggs2.v_h2_pt.clear();
      smhiggs2.v_h2_eta.clear();
      
      
      
      
       //------------------------GEN PARTICLES------------------------------ 	 
		 
	   // Get Thresholds for pT, eta for generated objects 
	 
	  // Looping overall Gen_Particles 	 
	  for (Int_t i = 0; i < Particle_size; i++){
		 
	        // Get particle pdg id 
	        int p_id = Particle_PID[i];
		 
	        // Check for Muons
	        if ( p_id == 13 ) { // pdg_id = 13 for Muon

		  Float_t gen_mu_pt = Particle_PT[i];
	          Float_t gen_mu_eta = Particle_Eta[i];
		  Float_t gen_mu_phi = Particle_Phi[i];
		       
		  genMuons.gen_muon_pt.push_back(gen_mu_pt);
		  genMuons.gen_muon_eta.push_back(gen_mu_eta);
		  genMuons.gen_muon_phi.push_back(gen_mu_phi);
		       
		      
	        } // end if on id 13 	
	     
	        // Check for b quarks
	        if ( p_id == 5 ) { // pdg_id = 5 for b quark 
			 
		  Float_t gen_b_pt = Particle_PT[i];
		  Float_t gen_b_eta = Particle_Eta[i];
		  Float_t gen_b_phi = Particle_Phi[i]; 
			
		  gen_bjet.gen_bjet_pt.push_back(gen_b_pt);
		  gen_bjet.gen_bjet_eta.push_back(gen_b_eta);
		  gen_bjet.gen_bjet_phi.push_back(gen_b_phi);
			 
			
	         } // end if on id 5
		 
           } // end loop over gen particles    
           
       
      
      //------------Loop over Loose Muons--------------
      for (Int_t i = 0; i < MuonLoose_size; i++){
	
	Float_t muon_pt = MuonLoose_PT[i];
	Float_t muon_eta = MuonLoose_Eta[i];
	Float_t muon_phi = MuonLoose_Phi[i];
			
	loose_mu.lmu_pt.push_back(muon_pt);
	loose_mu.lmu_eta.push_back(muon_eta);
	loose_mu.lmu_phi.push_back(muon_phi);
			
		
       } 
      
      //----------------Loop over Jets-----------------
       for ( Int_t i = 0; i < Jet_size; i++ ){
		   
	 Float_t j_pt = MuonLoose_PT[i];
	 Float_t j_eta = MuonLoose_Eta[i];
	 Float_t j_phi = MuonLoose_Phi[i];
		   
	 jets.jet_pt.push_back(j_pt);
	 jets.jet_eta.push_back(j_eta);
	 jets.jet_phi.push_back(j_phi);
		   
       }
      
      //---------------Loop over MET----------------
       
      for (Int_t i = 0; i < MissingET_size; i++){
		  
	Float_t missingET = MissingET_MET[i];
	Float_t MET_eta   = MissingET_Eta[i];
        Float_t MET_phi   = MissingET_Phi[i];
	        
	met.met_MET.push_back(missingET);
	met.met_eta.push_back(MET_eta);
	met.met_phi.push_back(MET_phi);
	       
       }
       
        
          
       //=====================================================================//
       //                        Start 4Muons , ZZ Selections                 //
       //=====================================================================//
	      
	    
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
         

       // 2nd Selection: on number of Muons per event

       if ( v_muon_idx.size() > 3 ){ // having at least 4 Muons per event 
			
      	 NEvents[0]++;
			  
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
      
         // In Delphes tree Muon_PT are sorted from highest to least one
         mu1.SetPtEtaPhiM(MuonLoose_PT[mu1_idx], MuonLoose_Eta[mu1_idx], MuonLoose_Phi[mu1_idx], muon_mass); 
         mu2.SetPtEtaPhiM(MuonLoose_PT[mu2_idx], MuonLoose_Eta[mu2_idx], MuonLoose_Phi[mu2_idx], muon_mass);
         mu3.SetPtEtaPhiM(MuonLoose_PT[mu3_idx], MuonLoose_Eta[mu3_idx], MuonLoose_Phi[mu3_idx], muon_mass);
         mu4.SetPtEtaPhiM(MuonLoose_PT[mu4_idx], MuonLoose_Eta[mu4_idx], MuonLoose_Phi[mu4_idx], muon_mass);  
      
         // Leading Muon pT > 20 GeV (Muon with highest pT)
         pt_mu1 = mu1.Pt();
         eta_mu1 = mu1.Eta();
         phi_mu1 = mu1.Phi();
                
         // subleading Muon pT > 10 GeV (Muon with second-highest pT)    
         pt_mu2 = mu2.Pt();
         eta_mu2 = mu2.Eta();
         phi_mu2 = mu2.Phi();
                
         pt_mu3 = mu3.Pt();
         eta_mu3 = mu3.Eta();
         phi_mu3 = mu3.Phi();
                
         pt_mu4 = mu2.Pt();
         eta_mu4 = mu2.Eta();
         phi_mu4 = mu2.Phi();
            
         four_muons_beforeCut.muon1_pt.push_back(pt_mu1);
         four_muons_beforeCut.muon2_pt.push_back(pt_mu2);
         four_muons_beforeCut.muon3_pt.push_back(pt_mu3);
         four_muons_beforeCut.muon4_pt.push_back(pt_mu4);
         four_muons_beforeCut.muon1_eta.push_back(eta_mu1);
         four_muons_beforeCut.muon2.push_back(eta_mu2);
         four_muons_beforeCut.muon3_eta.push_back(eta_mu3);
         four_muons_beforeCut.muon4_eta.push_back(eta_mu4);
         four_muons_beforeCut.muon1_phi.push_back(phi_mu1);
         four_muons_beforeCut.muon2_phi.push_back(phi_mu2);
         four_muons_beforeCut.muon3_phi.push_back(phi_mu3);
         four_muons_beforeCut.muon4_phi.push_back(phi_mu4);    
                
              
         // Calculate DR bet each 2 muons for all possible combinations 1234, 1324, 1423
         Float_t DR_mu12, DR_mu34, DR_mu13, DR_mu24, DR_mu14, DR_mu23; 
      
         // Initialize variables 
         DR_mu12 = -9999.; DR_mu34 = -9999.; DR_mu13 = -9999.; DR_mu24 = -9999.; DR_mu14 = -9999.; DR_mu23 = -9999.; 
      
         DR_mu12 = TMath::Sqrt(TMath::Power((eta_mu1 - eta_mu2), 2) + TMath::Power((phi_mu1 - phi_mu2), 2));
         DR_mu34 = TMath::Sqrt(TMath::Power((eta_mu3 - eta_mu4), 2) + TMath::Power((phi_mu3 - phi_mu4), 2));
         DR_mu13 = TMath::Sqrt(TMath::Power((eta_mu1 - eta_mu3), 2) + TMath::Power((phi_mu1 - phi_mu3), 2));
         DR_mu24 = TMath::Sqrt(TMath::Power((eta_mu2 - eta_mu4), 2) + TMath::Power((phi_mu2 - phi_mu4), 2));
         DR_mu14 = TMath::Sqrt(TMath::Power((eta_mu1 - eta_mu4), 2) + TMath::Power((phi_mu1 - phi_mu4), 2));
         DR_mu23 = TMath::Sqrt(TMath::Power((eta_mu2 - eta_mu3), 2) + TMath::Power((phi_mu2 - phi_mu3), 2));
      
         drMuons.v_DR_mu1mu2.push_back(DR_mu12);
         drMuons.v_DR_mu3mu4.push_back(DR_mu34);
         drMuons.v_DR_mu1mu3.push_back(DR_mu13);
         drMuons.v_DR_mu2mu4.push_back(DR_mu24);
         drMuons.v_DR_mu1mu4.push_back(DR_mu14);
         drMuons.v_DR_mu2mu3.push_back(DR_mu23); 
       
       
         /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           ^                                         ^ 
           ^            Determine Za, Zb             ^
           ^                                         ^
           ^              for h -> Z Z               ^
           ^                                         ^ 
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
      
      
          // 4 Muons various combination pairs 1234, 1324, 1423
          float mZ = 91.1876; // in GeV 
          Float_t mZ12, mZ34, mZ13, mZ24, mZ14, mZ23, dZc1, dZc2, dZc3;
          Float_t pt_Z12, pt_Z34, pt_Z13, pt_Z24, pt_Z14, pt_Z23;
          Float_t eta_Z12, eta_Z34, eta_Z13, eta_Z24, eta_Z14, eta_Z23;
          Float_t phi_Z12, phi_Z34, phi_Z13, phi_Z24, phi_Z14, phi_Z23;
          Float_t dZ12, dZ23, dZ34, dZ13, dZ14, dZ24;
          Float_t mZa, mZb, ptZa, ptZb, etaZa, etaZb, phiZa, phiZb;
      
          // Initialize variables
          mZ12 = -9999.; mZ34 = -9999.; mZ13 = -9999.; mZ24 = -9999.; mZ14 = -9999.; mZ23 = -9999.; dZc1 = -9999.; dZc2 = -9999.; dZc3 = -9999.;
          pt_Z12 = -9999.; pt_Z34 = -9999.; pt_Z13 = -9999.; pt_Z24 = -9999.; pt_Z14 = -9999.; pt_Z23 = -9999.;
          eta_Z12 = -9999.; eta_Z34 = -9999.; eta_Z13 = -9999.; eta_Z24 = -9999.; eta_Z14 = -9999.; eta_Z23 = -9999.;
          phi_Z12 = -9999.; phi_Z34 = -9999.; phi_Z13 = -9999.; phi_Z24 = -9999.; phi_Z14 = -9999.; phi_Z23 = -9999.;
          dZ12 = -9999.; dZ23 = -9999.; dZ34 = -9999.; dZ13 = -9999.; dZ14 = -9999.; dZ24 = -9999.;
          mZa = -9999.; mZb = -9999.; ptZa = -9999.; ptZb = -9999.; etaZa = -9999.; etaZb = -9999.; phiZa = -9999.; phiZb = -9999.;
      
      
          // 3rd Selection: on 4Muons charge
	    
          if ( MuonLoose_Charge[mu1_idx] + MuonLoose_Charge[mu2_idx] + MuonLoose_Charge[mu3_idx] + MuonLoose_Charge[mu4_idx] == 0){ //Sure that muon pairs are of opposite signs 
       
               NEvents[1]++; 
               
               // 4th Selection: on Leading (highest pT) and Subleading (second highest pT) Muons
	       
               if ( ( pt_mu1 > 20. ) && ( pt_mu2 > 10. ) ){         

                  NEvents[2]++;

		  if ( ( DR_mu12 > 0.02) && ( DR_mu34 > 0.02) ){

		    if ( ( DR_mu13 > 0.02) && ( DR_mu24 > 0.02) ){

		      if ( ( DR_mu14 > 0.02) && ( DR_mu23 > 0.02) ){

			 NEvents[4]++;
                   
                  /////////////////////////////////// 
		  //                               //
		  // Start 4 muon combination 1234 //
                  //                               //
		  ///////////////////////////////////
		  
		  // 5th Selection: on charge of each muon pair

		  if ( MuonLoose_Charge[mu1_idx] + MuonLoose_Charge[mu2_idx] == 0){  // mu1, mu2 
				   
                     if ( MuonLoose_Charge[mu3_idx] + MuonLoose_Charge[mu4_idx] == 0){  // mu3, mu4 
						  
			 NEvents[3]++;

                         // invariant mass for each muon pair
			 double m12 = (mu1 + mu2).M();
			 double m34 = (mu3 + mu4).M(); 
						  
			 // if (  ( DR_mu12 > 0.02) && ( DR_mu34 > 0.02) ){
						 
			     // NEvents[4]++;

			     // mMuMu> 4 GeV for the opposite charge lepton muon pair
                             
			     if (  ( m12 > 4.) &&  ( m34 > 4. ) ) {
			      
			       NEvents[5]++;
			       Z12 = mu1 + mu2;
                               mZ12 = Z12.M();
                               pt_Z12 = Z12.Pt();
                               eta_Z12 = Z12.Eta();
                               phi_Z12 = Z12.Phi();
                             
			       Z34 = mu3 + mu4;
                               mZ34 = Z34.M();
                               pt_Z34 = Z34.Pt();
                               eta_Z34 = Z34.Eta();
                               phi_Z34 = Z34.Phi();
			          
			     } // end if m12,m34
			// } // DR

			 if (mZ12 > 0.) zmass.Z12_mass.push_back(mZ12);
			 if (mZ34 > 0.) zmass.Z34_mass.push_back(mZ34);
			         
		      } // end if on mu3,4 charge
		    } // end if on mu1,2 charge
		  
		    dZ12 = fabs(mZ12 - mZ);
		    dZ34 = fabs(mZ34 - mZ);
		    
		    // condition ? result_if_true : result_if_false  -> syntax for using ? conditional operator 
		    dZc1 = ( dZ12 < dZ34 ) ? dZ12 : dZ34;
		          
                    ///////////////////////////////////
		    //                               //
		    // Start 4 muon combination 1324 //
		    //                               //
		    ///////////////////////////////////
		    
		    if ( MuonLoose_Charge[mu1_idx] + MuonLoose_Charge[mu3_idx] == 0){  // mu1, mu3
					
		       if ( MuonLoose_Charge[mu2_idx] + MuonLoose_Charge[mu4_idx] == 0){ // mu2, mu4
					   
			  NEvents[6]++;

                          // invariant mass for each muon pair
			  double m13 = (mu1 + mu3).M();
			  double m24 = (mu2 + mu4).M();
						 
			  //  if ( ( DR_mu13 > 0.02) && ( DR_mu24 > 0.02) ){
						 
			  //  NEvents[7]++;
						   
                              // mMuMu> 4 GeV for the opposite charge lepton muon pair
			      
			      if ( ( m13 > 4.) && ( m24 > 4. ) ) {  
							   
				 NEvents[8]++;
			         Z13 = mu1 + mu3;
			         mZ13 = Z13.M();
			         pt_Z13 = Z13.Pt();
                                 eta_Z13 = Z13.Eta();
                                 phi_Z13 = Z13.Phi();
                               
			         Z24 = mu2 + mu4;
			         mZ24 = Z24.M();
			         pt_Z24 = Z24.Pt();
                                 eta_Z24 = Z24.Eta();
                                 phi_Z24 = Z24.Phi();
                               
			      } // end if m13, m24

			 // } // DR
					    

		           if (mZ13 > 0.) zmass.Z13_mass.push_back(mZ13);
			   if (mZ24 > 0.) zmass.Z24_mass.push_back(mZ24);
			              
			} // end if on mu1,3 charge
	              } // end if on mu2,4 charge
		    
		      dZ13 = fabs(mZ13 - mZ);
		      dZ24 = fabs(mZ24 - mZ);
		
		      dZc2 = ( dZ13 < dZ24 ) ? dZ13 : dZ24; 
		      
                     ///////////////////////////////////
		     //                               //
		     // Start 4 muon combination 1423 //
		     //                               //
		     ///////////////////////////////////

		     if ( MuonLoose_Charge[0] + MuonLoose_Charge[3] == 0){  // mu1, mu4
				
		        if ( MuonLoose_Charge[1] + MuonLoose_Charge[2] == 0){ // mu2, mu3
						
			    NEvents[9]++;

                            // invariant mass for each muon pair
			    double m14 = (mu1 + mu4).M();
			    double m23 = (mu2 + mu3).M();
						
			    //  if ( ( DR_mu14 > 0.02) && ( DR_mu23 > 0.02) ){
						 
			    //  NEvents[10]++;
						   
                            // mMuMu> 4 GeV for the opposite charge lepton muon pair
			    if ( ( m14 > 4. ) && ( m23 > 4. ) ) {
							   
			       NEvents[11]++;
			       Z14 = mu1 + mu4;
			       mZ14 = Z14.M();
			       pt_Z14 = Z14.Pt();
			       eta_Z14 = Z14.Eta();
			       phi_Z14 = Z14.Phi();
			                   
			       Z23 = mu2 + mu3;
			       mZ23 = Z23.M();
			       pt_Z23 = Z23.Pt();
		               eta_Z23 = Z23.Eta();
			       phi_Z23 = Z23.Phi();
			                   
			     } // end if m14, m23
			 // } // DR

			    if (mZ14 > 0.) zmass.Z14_mass.push_back(mZ14);
			    if (mZ23 > 0.) zmass.Z23_mass.push_back(mZ23);
			              
			 } // end if on mu1,4 charge
		     } // end if on mu2,3 charge
		    
		     dZ14 = fabs(mZ14 - mZ);
		     dZ23 = fabs(mZ23 - mZ);
		
		     dZc3 = ( dZ14 < dZ23 ) ? dZ14 : dZ23; 
		      
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
			
                  
                    // 6th Selection: on Za and Zb masses, Za > 40 GeV (closest to nominal Z mass) & Zb > 12 GeV
                    
                    if ( mZa > 40. && mZa < 120.){
					  
		       if ( mZb > 12. && mZb < 120. ){
						 
			  NEvents[12]++;
						 
			  za_sig.Za_mass.push_back(mZa);
			  za_sig.Za_pt.push_back(ptZa);
			  za_sig.Za_eta.push_back(etaZa);
					     
			  zb_sig.Zb_mass.push_back(mZb);
			  zb_sig.Zb_pt.push_back(ptZb);
			  zb_sig.Zb_eta.push_back(etaZb);
		            
		            
		          //============================
		          // Reconstruct h1 from Za, Zb
		          //============================
            
                          Float_t mh1_ZaZb, pt_h1_ZaZb, eta_h1_ZaZb, phi_h1_ZaZb;
                   
                          // Initialize Variables
                          mh1_ZaZb = -9999.; pt_h1_ZaZb = -9999.; eta_h1_ZaZb = -9999.; phi_h1_ZaZb = -9999.;
                   
                          h1 = Za + Zb;
                          mh1_ZaZb = h1.M(); // get invariant mass of 4Muons
			  pt_h1_ZaZb = h1.Pt();
                          eta_h1_ZaZb = h1.Eta();
                          phi_h1_ZaZb = h1.Phi();
                         
                          // 7th Selection: on invariant mass of 4Muons, in Signal Region 115 ≤ m4l ≤ 135 GeV  
                          // if ( ( mh1_ZaZb >= 115. ) && ( mh1_ZaZb <= 135. ) ) {
							 
		          // NEvents[13]++;
                          // if ( ( mh1_ZaZb >= 115. ) && ( mh1_ZaZb <= 135. ) ) continue;
                         
                          // 8th Selection: on control region or side bands  m4l < 115 or m4l > 135 GeV
			  // if ( ( mh1_ZaZb < 115. ) || ( mh1_ZaZb > 135. ) ) continue;	
						  


			  /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                            ^                                         ^ 
                            ^            Determine b1, b2             ^
                            ^                                         ^
                            ^              for h -> b1 b2             ^
                            ^                                         ^ 
                            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
       
			  // bool found_bjet = false;

			    int nbjets = 0;            // total nb of b-jets found per event
                           
                            bjet_indx.clear();
            
                           
                            // Loop overall Jets
                            
                            for (Int_t i = 0; i < Jet_size; i++){
	 
	                        int bj_indx = i;
	                            
			        // Jet_BTag stores the bit value for each jet

                                //------------------------------------
				//   WP       BitNumber    Bit Values
                                //------------------------------------
				// Loose          0          1,3,5,7
				// Medium         1          2,3,6,7
				// Tight          2          4,5,6,7
				//------------------------------------

				UInt_t jet_bTag = Jet_BTag[i];
			        UInt_t jet_bTag_Algo = Jet_BTagAlgo[i];
			        UInt_t jet_bTag_Phys = Jet_BTagPhys[i];
                               
                                // Check whether a jet has passed the b-tagging criteria defined by the BitNumber,


				// Loose
                               /* Bool_t BtagOk = ( jet_bTag & (1 << 0) );    
                                
                                  cout << "Jet [" << i << "] ,  BTag bit value = " << jet_bTag 
                                       << ",   Is it considered as a Loose WP b-jet ?   Answer: " 
                                       << BtagOk << endl;  */
                 

                                // Tight 
				Bool_t BtagOk = ( jet_bTag & (1 << 2) );     
                                
                                cout << "Jet [" << i << "] ,  BTag bit value = " << jet_bTag 
                                     << ",   Is it considered as a Tight WP b-jet ?   Answer: " 
                                     << BtagOk << ",  " << "BTag_Algo = " << jet_bTag_Algo 
                                     << ",   BTag_Physics = " << jet_bTag_Phys << endl;
                                     
                                
                                // Medium      
                               /* Bool_t BtagOk = ( jet_bTag & (1 << 1) );   
                                
                                   cout << "Jet [" << i << "] ,  BTag bit value = " << jet_bTag 
                                        << ",   Is it considered as a Med WP b-jet ?   Answer: " 
                                        << BtagOk << endl;  */
			

				if ( BtagOk == 1){  
									
				   bjet_indx.push_back(bj_indx);
				   nbjets++;
									
				}  // end if BtagOk 
	                     }  // end loop overall jets
      
					 
		       	     // 2nd Selection: on having at least 2 b-jets/event 

			     if ( nbjets > 1 ){ 
								 
				NEvents[14]++;
								 
				cout << "========================================" << endl;
		 	        cout << " number of b-jets/event =  " << nbjets    << endl;
		                cout << "========================================" << endl;
							 		 
		        	bjet_indx_AfterSel.clear();

				
			        // Loop overall Tight b-jets/event 

				for ( Int_t i = 0; i < bjet_indx.size(); i++ ) { 
									 
				    int bj_indx_AfterSelec = bjet_indx[i];

				    double bjet_pt = Jet_PT[bj_indx_AfterSelec];
				    double bjet_eta = Jet_Eta[bj_indx_AfterSelec];
				    double abs_bjet_eta = fabs(Jet_Eta[bj_indx_AfterSelec]);
				    double bjet_phi = Jet_Phi[bj_indx_AfterSelec]; 


				    // get DR between b-jet & each of 4 muons
				    
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
										 
				       NEvents[15]++;

				       cout << "......Selection on DR between (bjet, lepton) DONE......"  << endl;


				       // 4th Selection: on b-jet pT > 20 GeV

				       if ( bjet_pt > 20.) {   
											 
					  NEvents[16]++;

					  cout << "......Selection on bjet PT DONE......"  << endl;


					  // 5th Selection on b-jet abs_eta < 2.4

					  if ( abs_bjet_eta < 2.4) {
												 
					     NEvents[17]++;

					     // save b-jet index for further selection 

					     bjet_indx_AfterSel.push_back(bj_indx_AfterSelec);

					     cout << "......Selection on bjet abs Eta DONE......"  << endl;

					  } // end absEta 
											
				        } // end if bjet_pt 	 

				      } // end if DR
				    
				 } // end loop over vector elements 
                                   
                                
				 // check if  we still have at least 2 bjets after above selections
				
				 if ( bjet_indx_AfterSel.size() > 1 ) {  
									    
				    NEvents[18]++;

				    cout << "bjet_indx_AfterSel size is : " <<  bjet_indx_AfterSel.size() << endl; 
                                      
											
				    // Select the signal 2 b-jets as the two with highest pT, no scores in Delphes 

				    int signal_bjet_1_indx = bjet_indx_AfterSel[0];
                                    int signal_bjet_2_indx = bjet_indx_AfterSel[1];

				    Float_t signal_bjet_1_pt = Jet_PT[signal_bjet_1_indx];
				    Float_t signal_bjet_2_pt = Jet_PT[signal_bjet_2_indx];
										    
				    cout << "2 B-jets of signal are ( " << signal_bjet_1_indx << ", " << signal_bjet_2_indx 
					 << " ) with Highest PT values " << signal_bjet_1_pt << ", " << signal_bjet_2_pt
					 << " GeV respectively!" << endl;
							          
							          
			            // Set TLorentzVectors for selected 2 b-jets of signal

				    b1.SetPtEtaPhiM( Jet_PT[signal_bjet_1_indx], 
						     Jet_Eta[signal_bjet_1_indx], 
						     Jet_Phi[signal_bjet_1_indx], 
						     Jet_Mass[signal_bjet_1_indx] );
							                     
				    b2.SetPtEtaPhiM( Jet_PT[signal_bjet_2_indx], 
						     Jet_Eta[signal_bjet_2_indx], 
						     Jet_Phi[signal_bjet_2_indx], 
						     Jet_Mass[signal_bjet_2_indx] );
							     
							     
				    Float_t bj1_mass =  b1.M();
			            Float_t bj1_pt   =  b1.Pt();
				    Float_t bj1_eta  =  b1.Eta();
				    Float_t bj1_phi  =  b1.Phi();
							     
				    Float_t bj2_mass =  b1.M();
				    Float_t bj2_pt   =  b1.Pt();
				    Float_t bj2_eta  =  b1.Eta();
				    Float_t bj2_phi  =  b1.Phi();
							     
				    Float_t DeltaEta_b1b2_sqr = TMath::Power( ( bj1_eta - bj2_eta ), 2); 
				    Float_t DeltaPhi_b1b2_sqr = TMath::Power( ( bj1_phi - bj2_phi ), 2); 
				    Float_t DR_b1b2 = TMath::Sqrt( DeltaEta_b1b2_sqr + DeltaPhi_b1b2_sqr );
							     
				    b1_sig.b1_mass.push_back(bj1_mass);
				    b1_sig.b1_pt.push_back(bj1_pt);
				    b1_sig.b1_eta.push_back(bj1_eta);
							            
				    b2_sig.b2_mass.push_back(bj2_mass);
				    b2_sig.b2_pt.push_back(bj2_pt);
				    b2_sig.b2_eta.push_back(bj2_eta);
							            
				    smhiggs1.v_h1_mass.push_back(mh1_ZaZb);
				    smhiggs1.v_h1_pt.push_back(pt_h1_ZaZb);
				    smhiggs1.v_h1_eta.push_back(eta_h1_ZaZb);
                                 

				    //========================================//
		                    //       Reconstruct h2 from b1, b2       //
		                    //========================================//
							     
				    // TLorentzVector for 2nd SM higgs
				    h2 = b1 + b2;
							     
				    Float_t h2_b1b2_mass =  h2.M();
				    Float_t h2_b1b2_pt   =  h2.Pt();
				    Float_t h2_b1b2_eta  =  h2.Eta();
				    Float_t h2_b1b2_phi  =  h2.Phi();
							     
				    smhiggs2.v_h2_mass.push_back(h2_b1b2_mass);
				    smhiggs2.v_h2_pt.push_back(h2_b1b2_pt);
				    smhiggs2.v_h2_eta.push_back(h2_b1b2_eta); 
							     
							     
                                    //======================================//
		                    //                                      //
		                    //      Reconstruct BSM H from SM h     // 
		                    //               H -> h h               //
		                    //                                      //
		                    //======================================//
                
                                    // TLorentzVector for BSM Heavy Higgs 
                                    H = h1 + h2;
                                 
                                    Float_t H_bb4Mu_mass =  H.M();
                                    Float_t H_bb4Mu_pt   =  H.Pt();
                                    Float_t H_bb4Mu_eta  =  H.Eta();
                                    Float_t H_bb4Mu_phi  =  H.Phi(); 
                                 
                                    heavyHiggs.v_H_mass.push_back(H_bb4Mu_mass); 
				    heavyHiggs.v_H_pt.push_back(H_bb4Mu_pt); 
				    heavyHiggs.v_H_eta.push_back(H_bb4Mu_eta); 
								       
								       
		                 } // end if bjet_indx_AfterSel.size() > 1 
			      } // end if nbjets > 1      
		               //  } // end if on m4l invariant mass selection  
		            } // end if mZb  
                          } // end if mZa
		      } // end DR mu1423
		    } // end DR mu1324
		  } // end DR mu1234

	       } // end if on leading and subleading muons pT 

	  } // end if on 4 muons charge 

       } // end if MuonLoose_size > 3
		     		  
   
    new_tree->SetWeight(wt);
    
    cout << "Filling tree" << endl;
    new_tree->Fill();


   } // end loop overall events
   
   
   // Efficiency = nSelectedEvents/nentries;
   
   cout << "***Analysis Loop Ends!***" << endl; 
   
   cout << " Total Number of Events is : " << nentries << " Events" << endl;
   cout << "     " << endl;  
   
   cout << "================================================================================"    << endl;
   cout << "                         Cut                             |  NEvents PASS Cut    |"   << endl;
   cout << "================================================================================"    << endl;
   cout << " number of muons/event > 3                               |   " <<  NEvents[0]        << endl;
   cout << "================================================================================"    << endl;
   cout << " 4 muon charge/event = 0                                 |   " <<  NEvents[1]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " leading, subleading muon pT > 20 & 10 GeV               |   " <<  NEvents[2]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " mu1mu2 , mu3mu4 pairs charge/event = 0 for comb. 1234   |   " <<  NEvents[3]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " DeltaR(mu1mu2, mu3mu4) > 0.02 for comb. 1234            |   " <<  NEvents[4]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " mass of mu1mu2 & mu3mu4/event > 4 GeV for comb. 1234    |   " <<  NEvents[5]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " mu1mu3 , mu2mu4 pairs charge/event = 0 for comb. 1324   |   " <<  NEvents[6]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " DeltaR(mu1mu3, mu2mu4) > 0.02 for comb. 1324            |   " <<  NEvents[7]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " mass of mu1mu3 & mu2mu4/event > 4 GeV for comb. 1324    |   " <<  NEvents[8]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " mu1mu4 , mu2mu3 pairs charge/event = 0 for comb. 1423   |   " <<  NEvents[9]        << endl; 
   cout << "================================================================================"    << endl;
   cout << " DeltaR(mu1mu4, mu2mu3) > 0.02 for comb. 1423            |   " <<  NEvents[10]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " mass of mu1mu4 & mu2mu3/event > 4 GeV for comb. 1423    |   " <<  NEvents[11]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " Za, Zb masses, Za > 40 GeV & Zb > 12 GeV                |   " <<  NEvents[12]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " invariant mass of 4Muons, in SR 115 ≤ m4l ≤ 135 GeV     |   " <<  NEvents[13]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " at least 2 b-jets/event  -  no cuts on b-jets           |   " <<  NEvents[14]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " DeltaR(b-jet,lepton) of ZZ candidates > 0.3             |   " <<  NEvents[15]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " b-jet pT > 20 GeV                                       |   " <<  NEvents[16]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " absolute b-jet Eta < 2.4                                |   " <<  NEvents[17]       << endl; 
   cout << "================================================================================"    << endl;
   cout << " at least 2 b-jets/event after applying above selections |   " <<  NEvents[18]       << endl; 
   cout << "================================================================================"    << endl;   
   
   
  /* cout << " Total Number of Events = " << nentries 
        << ", Number of Selected Events = " << nSelectedEvents
        << ", Efficiency = Nsel/Ngen = " << Efficiency << endl;
  */    
   
   cout << "Writing tree!" << endl;  
   new_tree->Write();
   
   cout << "Writing output file! " << endl;
   out_sig->Write();
  // out_ttbar->Write();
   //out_ZZ4Mu->Write();
   /* out_DY->Write(); 
    ZZbb2Mu->Write();  
    ZZ4b->Write();       
    ZWpbbMuNL->Write(); 
    ZWp2MuMuNL->Write(); 
    ZWmbbMuNL->Write();  
    ZWm2MuMuNL->Write(); 
    WW2Mu2NL->Write(); 
    smHZZ4Mu->Write();  
    smHbb->Write();  */
    
   cout << "saving..." << endl;    
   out_sig->Close();   
  // out_ttbar->Close();
   //out_ZZ4Mu->Close();
   /* out_DY->Close(); 
    ZZbb2Mu->Close();  
    ZZ4b->Close();       
    ZWpbbMuNL->Close(); 
    ZWp2MuMuNL->Close(); 
    ZWmbbMuNL->Close();  
    ZWm2MuMuNL->Close(); 
    WW2Mu2NL->Close(); 
    smHZZ4Mu->Close();  
    smHbb->Close();  */
   
   cout << "---DONE---" << endl;
   cout << "ROOT file: " << out_sig << " has been created sucessfully!" << endl;
   
   
}
