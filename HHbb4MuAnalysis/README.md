Hhhbb4Mu Analysis: 

1- Studying expectation for BSM second heavy CP even neutral Higgs production and decay into two SM like Higgs in the final state 2b-jets+4Muons at tree level via ggF for CMS Phase2 upgrade with Lumimosity 3000 1/fb . The Heavy Higgs mass to be investigated is ~ 400 GeV, due to a parameter space scan study performed on Z2 symmetry breaking parameter M12 of 2HDM-typeI model, which has obvious effect on heavy Higgs mass. Such that for high values of M12 we got lower masses for heavy Higgs. 


2- Signal and Background samples are generated using MG5+Pythia8+Delphes.
 
 
BTagging in CMS Phase2 Delphes Simulation Cards:

1- There three Working Points WP which are Loose, Medium and Tight for b-jets, 
   defined in Delphes card, each WP is stored in a certain BitNumber position,
   for e.g; BitNumber=0 (default) means stored in 1st position, BitNumber=1 
   stored in 2nd position, and so on.....
   
2- Each BitNumber stores certain bit values, specifically defined in CMS Phase2
   Card as :
   
   
       WP       BitNumber    Bit Values
   
     Loose          0          1,3,5,7
     Medium         1          2,3,6,7
     Tight          2          4,5,6,7
     
     
    Bit Value       0    1    2    3    4    5    6    7
   
    Binary form    000  001  010  011  100  101  110  111 
   


3- Each bit value is flagged to a certain WP or many.

      Bit_Value       WP_Loose        WP_Medium         WP_Tight
      
         0              NOT              NOT               NOT
         1            FLAGGED            NOT               NOT 
         2              NOT            FLAGGED             NOT
         3            FLAGGED          FLAGGED             NOT
         4              NOT              NOT             FLAGGED
         5            FLAGGED            NOT             FLAGGED
         6              NOT            FLAGGED           FLAGGED
         7            FLAGGED          FLAGGED           FLAGGED
         
         
4- Bit value stored in a variable called Jet_BTag[] for each jet in 
   Jet container in Delphes tree.  
   
   
5- Selection of certain WP B-jets for doing analysis mainly depends on BitNUmber value selection.
