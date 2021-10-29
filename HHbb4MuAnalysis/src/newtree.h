//============================================
//  Header file defining Tree and branches
//           for Hhhbb4M
//============================================

#ifndef newtree_h
#define newtree_h

#include <TROOT.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>


struct GeneratedMuons{
	
	
   vector<Float_t> gen_muon_pt;
   vector<Float_t> gen_muon_eta;
   vector<Float_t> gen_muon_phi;
   

};

GeneratedMuons genMuons;
/////////////////////////////////////////////////////

struct GeneratedBHadrons{
	
   vector<Float_t> gen_bjet_pt;
   vector<Float_t> gen_bjet_eta;
   vector<Float_t> gen_bjet_phi;
	
	
};
GeneratedBHadrons gen_bjet;

//////////////////////////////////////////////////////

struct LooseMuons{
	
   vector<Float_t> lmu_pt;
   vector<Float_t> lmu_eta;
   vector<Float_t> lmu_phi;
	

};
LooseMuons loose_mu;

//////////////////////////////////////////////////////

struct Jets{
	
   vector<Float_t> jet_pt;
   vector<Float_t> jet_eta;
   vector<Float_t> jet_phi;
	

};
Jets jets;

///////////////////////////////////////////////////////

struct MET{
	
   vector<Float_t> met_MET;
   vector<Float_t> met_eta;
   vector<Float_t> met_phi;

		
};
MET met;

//////////////////////////////////////////////////////

struct DeltaR_Muons{
	
   vector<Float_t> v_DR_mu1mu2;
   vector<Float_t> v_DR_mu3mu4;
   vector<Float_t> v_DR_mu1mu3;
   vector<Float_t> v_DR_mu2mu4;
   vector<Float_t> v_DR_mu1mu4;	
   vector<Float_t> v_DR_mu2mu3;	
	
};
DeltaR_Muons drMuons;

///////////////////////////////////////////////////////

struct FourMuons{
	
   vector<Float_t> muon1_pt;
   vector<Float_t> muon2_pt;
   vector<Float_t> muon3_pt;
   vector<Float_t> muon4_pt;
   vector<Float_t> muon1_eta;	
   vector<Float_t> muon2_eta;
   vector<Float_t> muon3_eta;	
   vector<Float_t> muon4_eta;	
   vector<Float_t> muon1_phi;
   vector<Float_t> muon2_phi;
   vector<Float_t> muon3_phi;
   vector<Float_t> muon4_phi;	
	
};
FourMuons four_muons;


///////////////////////////////////////////////////////

struct ZmassCombinations{
	
   vector<Float_t> Z12_mass;
   vector<Float_t> Z34_mass;
   vector<Float_t> Z13_mass;
   vector<Float_t> Z24_mass;
   vector<Float_t> Z14_mass;
   vector<Float_t> Z23_mass;	
	
		
};
ZmassCombinations zmass;

///////////////////////////////////////////////////////

struct Za_OfSignal{
	
   vector<Float_t> Za_mass;
   vector<Float_t> Za_pt;
   vector<Float_t> Za_eta;
   	
	
};
Za_OfSignal za_sig;

////////////////////////////////////////////////////////

struct Zb_OfSignal{
	
   vector<Float_t> Zb_mass;
   vector<Float_t> Zb_pt;
   vector<Float_t> Zb_eta;
   	
	
};
Zb_OfSignal zb_sig;

////////////////////////////////////////////////////////

struct b1_OfSignal{
	
   vector<Float_t> b1_mass;
   vector<Float_t> b1_pt;
   vector<Float_t> b1_eta;
   	
	
};
b1_OfSignal b1_sig;

////////////////////////////////////////////////////////

struct b2_OfSignal{
	
   vector<Float_t> b2_mass;
   vector<Float_t> b2_pt;
   vector<Float_t> b2_eta;
   	
	
};
b2_OfSignal b2_sig;

//////////////////////////////////////////////////////////

struct FirstHiggsOfSignal{
	
   vector<Float_t> v_h1_mass;
   vector<Float_t> v_h1_pt;
   vector<Float_t> v_h1_eta;
		
}; 
FirstHiggsOfSignal smhiggs1;

///////////////////////////////////////////////////////////

struct SecondHiggsOfSignal{
	
   vector<Float_t> v_h2_mass;
   vector<Float_t> v_h2_pt;
   vector<Float_t> v_h2_eta;
		
}; 
SecondHiggsOfSignal smhiggs2;

//////////////////////////////////////////////////////////

struct BSM_HiggsOfSignal{
	
   vector<Float_t> v_H_mass;
   vector<Float_t> v_H_pt;
   vector<Float_t> v_H_eta;
		
}; 
BSM_HiggsOfSignal heavyHiggs;

//////////////////////////////////////////////////////////


TTree* new_tree  = new TTree("output_demo", "output_demo");


TBranch*   b_Gen_Muon;
TBranch*   b_Gen_Bjets;
TBranch*   b_LooseMuons;
TBranch*   b_Jets;
TBranch*   b_MET;
TBranch*   b_4Muons;
TBranch*   b_deltaR_muons;
TBranch*   b_Zmass_comb;
TBranch*   b_Za_signal;
TBranch*   b_Zb_signal;
TBranch*   b_b1_signal;
TBranch*   b_b2_signal;
TBranch*   b_1st_Higgs_Ofsignal;
TBranch*   b_2nd_Higgs_Ofsignal;
TBranch*   b_BSM_Higgs_Ofsignal;


#endif
