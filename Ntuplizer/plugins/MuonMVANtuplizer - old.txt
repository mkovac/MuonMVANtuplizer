// -*- C++ -*-
//
// Package:    Analyzer/Test
// Class:      Test
//
/**\class Test Test.cc Analyzer/Test/plugins/Test.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marko Kovac
//         Created:  Wed, 23 Jan 2019 13:58:45 GMT
//
//


// User include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "RecoEgamma/EgammaTools/interface/MVAVariableManager.h" // Specify variables using .txt file
#include "RecoEgamma/EgammaTools/interface/MultiToken.h"


#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "MuonMVAReader/Reader/interface/MuonGBRForestReader.hpp"

#include "TTree.h"
#include "TFile.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

int YEAR = 2018;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace reco;
using namespace edm;
using namespace std;


class MuonMVANtuplizer:public edm::one::EDAnalyzer<edm::one::SharedResources> {
   
   public:
      
      explicit MuonMVANtuplizer(const edm::ParameterSet&);
      ~MuonMVANtuplizer() override;
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
   
      void analyze(const edm::Event &, const edm::EventSetup &) override;

      // Method called once each job just before starting event loop
      void beginJob() override {};
      
      // Method called once each job just after ending the event loop
      void endJob() override {};

      // MC matching
      template<class T, class V>
      int matchToTruth(const T &el, const V &genParticles, int &genIdx);


      // Other
      TTree* tree_;

//      MVAVariableManager<reco::GsfElectron> mvaVarMngr_;
//      std::vector<float> vars_;
//      int nVars_;

      // Global variables
      int nEvent_, nRun_, nLumi_;
      int genNpu_;
      int vtxN_;
      int n_muons_;

      // Muon variables
      double mu_Q_;
      double mu_pT_;
      double mu_eta_;
      double mu_phi_;
      double mu_E_;
      double Z_mass_;

      int matchedToGenMu_;
      int matchedToGenMu_POG_codes_;
      int matchedToGenMuPOG_;
      int matchedGenIdx_;
      int mu_mother_;
      
      bool is_global_mu_;
      bool is_tracker_mu_;
      bool is_pf_mu_;
      
      int mu_N_hits_;
      int mu_N_matched_stations_;
      int mu_N_pixel_hits_;
      int mu_N_tracker_hits_;
      
      double mu_chi_square_;
      
      // Vertexing
      double mu_ip_;
      double mu_ip_error_;
      double mu_sip_;
      double mu_dxy_;
      double mu_dz_;
      
      double mu_pf_charged_had_iso_;
      double mu_pf_neutral_had_iso_;
      double mu_pf_photon_iso_;
      double mu_pu_charged_had_iso_;
      double mu_rho_;
      
      double mu_dR_;
      float score_;


      // ID decisions and categories
      std::vector<int> mvaPasses_;
      std::vector<float> mvaValues_;
      std::vector<int> mvaCats_;

      // Config
      const bool is_MC_;
      const double deltaR_;
      const double pT_threshold_;
      int setup_;
      
      TLorentzVector *mu_p_ = new TLorentzVector();
      TLorentzVector *mu_m_ = new TLorentzVector();


      // ID decisions objects
//      const std::vector< std::string > eleMapTags_;
//      std::vector< edm::EDGetTokenT< edm::ValueMap<bool> > > eleMapTokens_;
//      const std::vector< std::string > eleMapBranchNames_;
//      const size_t nEleMaps_;

      // MVA values and categories (optional)
//      const std::vector< std::string > valMapTags_;
//      std::vector< edm::EDGetTokenT<edm::ValueMap<float> > > valMapTokens_;
//      const std::vector< std::string > valMapBranchNames_;
//      const size_t nValMaps_;

//      const std::vector< std::string > mvaCatTags_;
//      std::vector< edm::EDGetTokenT<edm::ValueMap<int> > > mvaCatTokens_;
//      const std::vector< std::string > mvaCatBranchNames_;
//      const size_t nCats_;

      // Tokens for AOD and MiniAOD case
      MultiTokenT<edm::View<pat::Muon>>           src_;
      MultiTokenT<std::vector<reco::Vertex>>      vertices_;
      MultiTokenT<std::vector<PileupSummaryInfo>> pileup_;
//      MultiTokenT<std::double>                    rho_;
      MultiTokenT<edm::View<reco::GenParticle>>   genParticles_;
      edm::EDGetTokenT<double> rhoToken;
      
      // MVA reader
      MuonGBRForestReader *r;

};


enum MuonMatchType{
   UNMATCHED,
   TRUE_PROMPT_MUON,
   TRUE_MUON_FROM_TAU,
   TRUE_NON_PROMPT_MUON,
   }; // The last does not include tau parents


// Constructors and destructor
MuonMVANtuplizer::MuonMVANtuplizer(const edm::ParameterSet& iConfig)
//  : mvaVarMngr_            (iConfig.getParameter<std::string>("variableDefinition"))
  : is_MC_                 (iConfig.getParameter<bool>("isMC"))
  , deltaR_                (iConfig.getParameter<double>("deltaR"))
  , pT_threshold_          (iConfig.getParameter<double>("ptThreshold"))
//  , setup_                 (iConfig.getParameter<int>("setup"))
//  , eleMapTags_            (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVAs"))
//  , eleMapBranchNames_     (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVALabels"))
//  , nEleMaps_              (eleMapBranchNames_.size())
//  , valMapTags_            (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVAValMaps"))
//  , valMapBranchNames_     (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVAValMapLabels"))
//  , nValMaps_              (valMapBranchNames_.size())
//  , mvaCatTags_            (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVACats"))
//  , mvaCatBranchNames_     (iConfig.getUntrackedParameter<std::vector<std::string>>("eleMVACatLabels"))
//  , nCats_                 (mvaCatBranchNames_.size())
  , src_                   (consumesCollector(), iConfig, "srcMiniAOD")
  , vertices_              (src_, consumesCollector(), iConfig, "verticesMiniAOD")
  , pileup_                (src_, consumesCollector(), iConfig, "pileupMiniAOD")
//  , rho_                   (src_, consumesCollector(), iConfig, "rhoMiniAOD")
  , genParticles_          (src_, consumesCollector(), iConfig, "genParticlesMiniAOD")
{


   rhoToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));


/*
   // Muon maps
   for (size_t k = 0; k < nEleMaps_; ++k)
   {
   	  eleMapTokens_.push_back(consumes<edm::ValueMap<bool> >(edm::InputTag(eleMapTags_[k])));
   	  
   	  // Initialize vectors for holding ID decisions
      mvaPasses_.push_back(0);
   }

    // valMaps
    for (size_t k = 0; k < nValMaps_; ++k) {
        valMapTokens_.push_back(consumes<edm::ValueMap<float> >(edm::InputTag(valMapTags_[k])));

        // Initialize vectors for holding MVA values
        mvaValues_.push_back(0.0);
    }

    // categories
    for (size_t k = 0; k < nCats_; ++k) {
        mvaCatTokens_.push_back(consumes<edm::ValueMap<int> >(edm::InputTag(mvaCatTags_[k])));

        // Initialize vectors for holding MVA values
        mvaCats_.push_back(0);
    }
*/
    
    
//==========
// Book tree
//==========
   
   usesResource(TFileService::kSharedResource);
   edm::Service<TFileService> fs;
   tree_ = fs->make<TTree>("tree", "tree");

//   nVars_ = mvaVarMngr_.getNVars();

   tree_->Branch("nEvent",  &nEvent_);
   tree_->Branch("nRun",    &nRun_);
   tree_->Branch("nLumi",   &nLumi_);
   tree_->Branch("n_muons", &n_muons_);
   
   if (is_MC_) tree_->Branch("genNpu", &genNpu_);
   
   tree_->Branch("vtxN",   &vtxN_);

   tree_->Branch("mu_Q",   &mu_Q_);
   tree_->Branch("mu_pT",  &mu_pT_);
   tree_->Branch("mu_eta", &mu_eta_);
   tree_->Branch("mu_phi", &mu_phi_);
   tree_->Branch("mu_E",   &mu_E_);
   tree_->Branch("Z_mass", &Z_mass_);
   
   tree_->Branch("is_global_mu",  &is_global_mu_);
   tree_->Branch("is_tracker_mu", &is_tracker_mu_);
   tree_->Branch("is_pf_mu",      &is_pf_mu_);
   
   tree_->Branch("mu_N_hits",             &mu_N_hits_);
   tree_->Branch("mu_N_matched_stations", &mu_N_matched_stations_);
   tree_->Branch("mu_N_pixel_hits",       &mu_N_pixel_hits_);
   tree_->Branch("mu_N_tracker_hits",     &mu_N_tracker_hits_);
   tree_->Branch("mu_chi_square",         &mu_chi_square_);
   
   tree_->Branch("mu_ip",       &mu_ip_);
   tree_->Branch("mu_ip_error", &mu_ip_error_);
   tree_->Branch("mu_sip",      &mu_sip_);
   tree_->Branch("mu_dxy",      &mu_dxy_);
   tree_->Branch("mu_dz",       &mu_dz_);
   
   tree_->Branch("mu_pf_charged_had_iso", &mu_pf_charged_had_iso_);
   tree_->Branch("mu_pf_neutral_had_iso", &mu_pf_neutral_had_iso_);
   tree_->Branch("mu_pf_photon_iso",      &mu_pf_photon_iso_);
   tree_->Branch("mu_pu_charged_had_iso", &mu_pu_charged_had_iso_);
   tree_->Branch("mu_rho", &mu_rho_);
   tree_->Branch("mu_dR", &mu_dR_);
   
   tree_->Branch("score", &score_);


   if (is_MC_)
   {
       tree_->Branch("matchedToGenMu", &matchedToGenMu_);
       tree_->Branch("matchedToGenMuPOG", &matchedToGenMuPOG_);
       tree_->Branch("matchedToGenMu_POG_codes", &matchedToGenMu_POG_codes_);
   }

   // Has to be in two different loops
//   for (int i = 0; i < nVars_; ++i) {
//       vars_.push_back(0.0);
//   }
//   for (int i = 0; i < nVars_; ++i) {
//       tree_->Branch(mvaVarMngr_.getName(i).c_str(), &vars_[i]);
//   }

   // IDs
//   for (size_t k = 0; k < nValMaps_; ++k) {
//       tree_->Branch(valMapBranchNames_[k].c_str() ,  &mvaValues_[k]);
//   }
//
//   for (size_t k = 0; k < nEleMaps_; ++k) {
//       tree_->Branch(eleMapBranchNames_[k].c_str() ,  &mvaPasses_[k]);
//   }
//
//   for (size_t k = 0; k < nCats_; ++k) {
//       tree_->Branch(mvaCatBranchNames_[k].c_str() ,  &mvaCats_[k]);
//   }

   // All tokens for event content needed by this MVA
   // Tags from the variable helper
//   mvaVarMngr_.setConsumes(consumesCollector());

   r = new MuonGBRForestReader(YEAR);

}


// Destructor
//===================================
MuonMVANtuplizer::~MuonMVANtuplizer()
{
}
//===================================


//=====================================================================================
void MuonMVANtuplizer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
   
   // Fill global event info
   nEvent_ = iEvent.id().event();
   nRun_   = iEvent.id().run();
   nLumi_  = iEvent.luminosityBlock();

   // Get Handles
   auto src          = src_.getValidHandle(iEvent);
   auto vertices     = vertices_.getValidHandle(iEvent);

   // Get MC only Handles, which are allowed to be non-valid
   auto genParticles = genParticles_.getHandle(iEvent);
   auto pileup       = pileup_.getHandle(iEvent);
//   auto rho          = rho_.getHandle(iEvent);

   // Rho (for isolation)
   edm::Handle<double> rhoHandle;
   iEvent.getByToken(rhoToken, rhoHandle);
   mu_rho_ = *rhoHandle;
     
     
   vtxN_ = vertices->size();

   // Fill with true number of pileup
   if (is_MC_)
   {
      for (const auto &pu : *pileup)
      {
         int bx = pu.getBunchCrossing();
         if (bx == 0)
         {
            genNpu_ = pu.getPU_NumInteractions();
            break;
         }
      }
   }

    // Get MVA decisions
//    edm::Handle<edm::ValueMap<bool> > decisions[nEleMaps_];
//    for (size_t k = 0; k < nEleMaps_; ++k) {
//        iEvent.getByToken(eleMapTokens_[k],decisions[k]);
//    }

    // Get MVA values
//    edm::Handle<edm::ValueMap<float> > values[nValMaps_];
//    for (size_t k = 0; k < nValMaps_; ++k) {
//        iEvent.getByToken(valMapTokens_[k],values[k]);
//    }

    // Get MVA categories
//    edm::Handle<edm::ValueMap<int> > mvaCats[nCats_];
//    for (size_t k = 0; k < nCats_; ++k) {
//        iEvent.getByToken(mvaCatTokens_[k],mvaCats[k]);
//    }

//	cout << "nEvent = " << nEvent_ << endl;
//	cout << src->size() << endl;
	

   for (size_t iMu = 0; iMu < src->size(); ++iMu)
   {
   	
//		if ( src->size() != 2 ) continue;
   	
      const auto mu = src->ptrAt(iMu);
      
      mu_Q_   = mu->charge();
      mu_pT_  = mu->pt();
      mu_eta_ = mu->eta();
      mu_phi_ = mu->phi();
      mu_E_   = mu->energy();
      
      n_muons_ = src->size();
      
      

      if (mu->pt() < pT_threshold_) continue;

// Fill variables from .txt file
//        for (int iVar = 0; iVar < nVars_; ++iVar) {
//            vars_[iVar] = mvaVarMngr_.getValue(iVar, ele, iEvent);
//        }

      if (is_MC_)
      {
         matchedToGenMu_ = matchToTruth( mu, genParticles, matchedGenIdx_);
         
         if ( matchedToGenMu_ == UNMATCHED )
         {
         	matchedToGenMu_POG_codes_ = 0;
        	}
        	else if( matchedToGenMu_ == TRUE_PROMPT_MUON )
        	{
        		matchedToGenMu_POG_codes_ = 4;
        	}
        	else if( matchedToGenMu_ == TRUE_MUON_FROM_TAU )
        	{
        		matchedToGenMu_POG_codes_ = 3;
        	}
        	else if( matchedToGenMu_ == TRUE_NON_PROMPT_MUON )
        	{
        		matchedToGenMu_POG_codes_ = 2;
        	}
         
         matchedToGenMuPOG_ = mu->simType();
         
//         cout << "matchedToGenMu_ = " << matchedToGenMu_ << endl;
//         cout << "matchedToGenMu_POG_codes_ = " << matchedToGenMu_POG_codes_ << endl;
//         cout << "Muon POG = " << matchedToGenMuPOG_ << endl;
      }
      
      is_global_mu_  = mu->isGlobalMuon();
      is_tracker_mu_ = mu->isTrackerMuon();
      is_pf_mu_      = mu->isPFMuon(); // Particle-Flow muon id

      if (is_global_mu_)
      {
         // Number of muon chamber hits included in the the global muon track fit
         mu_N_hits_ = (mu->globalTrack()->hitPattern().numberOfValidMuonHits());
         
         // Chi2 of the global track fit
         mu_chi_square_ = (mu->globalTrack()->normalizedChi2());
      }
      else
      {
         mu_N_hits_     = -1;
         mu_chi_square_ = -1;
      }
      
      // Muon segments in at least two muon stations. This implies that the muon is also an arbitrated tracker muon.
      mu_N_matched_stations_ = mu->numberOfMatchedStations();
      
      // Number of hits in the pixel detector
      bool valid_KF = false;
      reco::TrackRef myTrackRef = mu->innerTrack();
      valid_KF = (myTrackRef.isAvailable());
      valid_KF = (myTrackRef.isNonnull());


      if (valid_KF)
      {
         // Number of pixel hits
         mu_N_pixel_hits_ = mu->innerTrack()->hitPattern().numberOfValidPixelHits();
      
         // Number of hits in the tracker layers
         mu_N_tracker_hits_ = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
      }
      else
      {
         mu_N_pixel_hits_ = -1;
         mu_N_tracker_hits_ = -1;
      }
    
      // Vertexing
      double IP       = fabs(mu->dB(pat::Muon::PV3D));
      double IP_error = mu->edB(pat::Muon::PV3D);
      double SIP      = -999;
      
      if (IP_error != 0)
      {
         SIP =  IP/IP_error;
      }
      
      mu_ip_ = IP;
      mu_ip_error_ = IP_error;
      mu_sip_ = SIP;
      

//      cout << "IP = " << mu_ip_ << "   " << "SIP = " << mu_sip_ << "   " << "dR = " << mu_dR_ << endl;

    
      double dxy = -999.;
      double dz  = -999.;
    
      const Vertex *vertex = 0;
      
      if (vertices->size() > 0)
      {
         vertex = &(vertices->front());
         dxy    = fabs(mu->muonBestTrack()->dxy(vertex->position()));
         dz     = fabs(mu->muonBestTrack()->dz(vertex->position()));
      }
    
      mu_dxy_ = dxy; 
      mu_dz_ = dz; 

      // Isolation   
      mu_pf_charged_had_iso_ = mu->pfIsolationR03().sumChargedHadronPt;
      mu_pf_neutral_had_iso_ = mu->pfIsolationR03().sumNeutralHadronEt;
      mu_pf_photon_iso_      = mu->pfIsolationR03().sumPhotonEt;
      mu_pu_charged_had_iso_ = mu->pfIsolationR03().sumPUPt;


// Look up and save the ID decisions
//      for (size_t k = 0; k < nEleMaps_; ++k)
//      {
//         mvaPasses_[k] = (int)(*decisions[k])[ele];
//      }
//
//      for (size_t k = 0; k < nValMaps_; ++k)
//      {
//         mvaValues_[k] = (*values[k])[ele];
//      }
//
//      for (size_t k = 0; k < nCats_; ++k)
//      {
//         mvaCats_[k] = (*mvaCats[k])[ele];
//      }


		if ( mu_Q_ == 1 )  mu_p_->SetPtEtaPhiE(mu_pT_, mu_eta_, mu_phi_, mu_E_);
		if ( mu_Q_ == -1 ) mu_m_->SetPtEtaPhiE(mu_pT_, mu_eta_, mu_phi_, mu_E_);		

      
      score_ = r->Get_MVA_value(mu_pT_, mu_eta_, mu_N_hits_, mu_N_pixel_hits_, mu_N_tracker_hits_, mu_chi_square_, mu_pf_photon_iso_,
                          mu_pf_charged_had_iso_, mu_pf_neutral_had_iso_, mu_rho_, mu_sip_, mu_dxy_, mu_dz_);
                          
      tree_->Fill();
 
   }
   
//   auto Z_ = (*mu_p_) + (*mu_m_);
//   cout << "Z mass = " << Z_.M() << endl;

   
}


template<class T, class V>
int MuonMVANtuplizer::matchToTruth(const T &mu, const V &genParticles, int &genIdx)
{

   // Find the closest status 1 gen muon to the reco muon
   double dR = 999;
   genIdx = -1;
  
   for (size_t i = 0; i < genParticles->size(); i++)
   {
      const auto particle = genParticles->ptrAt(i);
      
      // Drop everything that is not muon or not status 1
      if (abs(particle->pdgId()) != 13 || particle->status() != 1 ) continue;
   
      double dRtmp = ROOT::Math::VectorUtil::DeltaR( mu->p4(), particle->p4() );
      
//      cout << "dR = " << dRtmp << endl;
      
      if (dRtmp < dR)
      {
         dR = dRtmp;
         genIdx = i;
      }
   }
  
  	mu_dR_ = dR;
  
   // See if the closest muon is close enough. If not, no match found.
   if (genIdx == -1 || dR >= deltaR_)
   {
      return UNMATCHED;
   }

   const auto closestMuon = genParticles->ptrAt(genIdx);

//   if (closestMuon->fromHardProcessFinalState())
   if (closestMuon->isPromptFinalState())
   {
   	return TRUE_PROMPT_MUON;
   }
  
   if (closestMuon->isDirectHardProcessTauDecayProductFinalState())
   {
      return TRUE_MUON_FROM_TAU;
   }
  
   // What remains is true non-prompt muons
   return TRUE_NON_PROMPT_MUON;
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonMVANtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

    edm::ParameterSetDescription desc;
//    desc.add<edm::InputTag>("src");
//    desc.add<edm::InputTag>("vertices");
//    desc.add<edm::InputTag>("pileup");
//    desc.add<edm::InputTag>("genParticles");
    desc.add<edm::InputTag>("srcMiniAOD");
    desc.add<edm::InputTag>("verticesMiniAOD");
    desc.add<edm::InputTag>("pileupMiniAOD");
    desc.add<edm::InputTag>("genParticlesMiniAOD");
//    desc.add<std::string>("variableDefinition");
    desc.add<bool>("isMC");
    desc.add<double>("deltaR", 0.1);
    desc.add<double>("ptThreshold", 5.0);
//    desc.addUntracked<std::vector<std::string>>("eleMVAs");
//    desc.addUntracked<std::vector<std::string>>("eleMVALabels");
//    desc.addUntracked<std::vector<std::string>>("eleMVAValMaps");
//    desc.addUntracked<std::vector<std::string>>("eleMVAValMapLabels");
//    desc.addUntracked<std::vector<std::string>>("eleMVACats");
//    desc.addUntracked<std::vector<std::string>>("eleMVACatLabels");
    descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonMVANtuplizer);
