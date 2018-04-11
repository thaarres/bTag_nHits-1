
// system include files
#include <memory>
#include <string> 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// For event setup
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

// ROOT includes
#include "TLorentzVector.h" 
#include "TMath.h" 
#include "TTree.h" 

// For output
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// For btag
#include "DataFormats/BTauReco/interface/JetTag.h"


// For jet
#include "DataFormats/JetReco/interface/PFJet.h" 
#include "DataFormats/JetReco/interface/PFJetCollection.h" 

// #include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
// #include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
// #include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"

// For pixel clusters
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"

// Pixel topology
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#define NEW_ID
#ifdef NEW_ID
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#else 
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#endif 


#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"


//  For gen particle
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// For vertices
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"



class HitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit HitAnalyzer(const edm::ParameterSet& conf);
  ~HitAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  //const std::vector<reco::Candidate> findDaughters(const reco::GenParticle *particle); 

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void reset( void );
  const reco::GenParticle* findMother(const reco::GenParticle *particle);
      
  edm::ParameterSet conf_;
  edm::InputTag src_;
  bool printLocal;
  bool phase1_;
       
  edm::EDGetTokenT< reco::GenParticleCollection>          genPtoken;
  edm::EDGetTokenT< reco::PFJetCollection >               ak4CHStoken;
  edm::EDGetTokenT< reco::PFJetCollection >               ak8CHStoken;
  edm::EDGetTokenT< edmNew::DetSetVector<SiPixelCluster>> clusterToken;
  edm::EDGetTokenT< reco::VertexCollection >              svToken;
  edm::EDGetTokenT< reco::VertexCollection >              pvToken;
  edm::EDGetTokenT< reco::TrackCollection >               trackToken;
  edm::EDGetTokenT< reco::JetTagCollection >              csv2Token;
     
  //
  // vector<reco::Vertex>                  "inclusiveSecondaryVertices"   ""                "RECO"
  //   vector<reco::TrackExtra>              "generalTracks"             ""                "RECO"
  //     vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >    "pfGhostTrackVertexTagInfos"   ""                "RECO"
  //     vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >    "pfInclusiveSecondaryVertexFinderCvsLTagInfos"   ""                "RECO"
  //     vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >    "pfInclusiveSecondaryVertexFinderTagInfos"   ""                "RECO"
  //     vector<reco::TemplatedSecondaryVertexTagInfo<reco::IPTagInfo<vector<edm::Ptr<reco::Candidate> >,reco::JetTagInfo>,reco::VertexCompositePtrCandidate> >    "pfSecondaryVertexTagInfos"   ""                "RECO"
  //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfSimpleInclusiveSecondaryVertexHighEffBJetTags"   ""                "RECO"
  //           edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>    "pfSimpleSecondaryVertexHighEffBJetTags"   ""                "RECO"

  // ----------member data ---------------------------
      
  int                  nJets;
  std::vector<double>  jet_pt;
  std::vector<double>  jet_eta;
  std::vector<double>  jet_phi;
  std::vector<double>  jet_mass;
  std::vector<int>     jet_pdgId;
  std::vector<double>  jet_bTag;
  // std::vector<double>  jet_vx_x   ;
  // std::vector<double>  jet_vx_y   ;
  // std::vector<double>  jet_vx_z   ;
  // std::vector<double>  jet_vx_eta ;
  // std::vector<double>  jet_vx_phi ;
  // std::vector<double>  jet_vx_r   ;
  
  
  
  std::vector<double>  dr_jetGen;
       
  int                  nGenParticles;
  std::vector<double>  genParticle_pt;
  std::vector<double>  genParticle_eta;
  std::vector<double>  genParticle_phi;
  std::vector<double>  genParticle_mass;
  std::vector<int>     genParticle_pdgId;
  std::vector<int>     genParticle_mother_pdgId;
  std::vector<int>     genParticle_status;
  std::vector<double>  genParticle_vx_x   ;
  std::vector<double>  genParticle_vx_y   ;
  std::vector<double>  genParticle_vx_z   ;
  std::vector<double>  genParticle_vx_eta ;
  std::vector<double>  genParticle_vx_phi ;
  std::vector<double>  genParticle_vx_r   ;
  std::vector<double>  genParticle_decayvx_x   ;
  std::vector<double>  genParticle_decayvx_y   ;
  std::vector<double>  genParticle_decayvx_z   ;
  std::vector<double>  genParticle_decayvx_eta ;
  std::vector<double>  genParticle_decayvx_phi ;
  std::vector<double>  genParticle_decayvx_r   ;
  

  int                  nTracks;
  std::vector<double>  track_pt;
  std::vector<double>  track_eta;
  std::vector<double>  track_phi;
  //std::vector<double>  track_mass;
  std::vector<double>  track_charge   ;
  std::vector<double>  track_px   ;
  std::vector<double>  track_py   ;
  std::vector<double>  track_pz ;
  std::vector<double>  track_vx_x   ;
  std::vector<double>  track_vx_y   ;
  std::vector<double>  track_vx_z   ;
  //std::vector<double>  track_vx_eta ;
  std::vector<double>  track_vx_phi ;
  std::vector<double>  track_vx_r   ;
  
  
       
  int                          nDetUnits;
  // std::vector<unsigned int>    detUnit_detType;
  std::vector<unsigned int>    detUnit_subdetId;
       
  // barrel ids
  // std::vector< unsigned int>  detUnit_layerC;
 //  std::vector< unsigned int>  detUnit_ladderC;
 //  std::vector< unsigned int>  detUnit_zindex;
 //  std::vector< int>           detUnit_shell ; // shell id // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
 //  std::vector< int>           detUnit_sector; // 1-8
  // std::vector< int>           detUnit_ladder; // 1-22
  std::vector< int>           detUnit_layer ; // 1-3
  // std::vector< int>           detUnit_module; // 1-4
  // std::vector< bool>          detUnit_half ;

  // Endcap ids
  std::vector<unsigned int>   detUnit_disk;   //1,2,3
  // std::vector<unsigned int>   detUnit_blade;  //1-24
  // std::vector<unsigned int>   detUnit_moduleF; // plaquette 1,2,3,4
  std::vector<unsigned int>   detUnit_side;   //size=1 for -z, 2 for +z
  // std::vector<unsigned int>   detUnit_panel;  //panel=1
  
   std::vector< double >           detUnit_X     ;
   std::vector< double >           detUnit_Y     ;
   std::vector< double >           detUnit_Z     ;
   std::vector< double >           detUnit_R     ;
   std::vector< double >           detUnit_Phi   ;
   
   
   
   
       
  std::vector<int>                  nClusters;
  std::vector<std::vector<double>>  cluster_sizex;
  std::vector<std::vector<double>>  cluster_sizey;
  std::vector<std::vector<double>>  cluster_globalz;
  std::vector<std::vector<double>>  cluster_globalx;
  std::vector<std::vector<double>>  cluster_globaly;
  std::vector<std::vector<double>>  cluster_globalPhi;
  std::vector<std::vector<double>>  cluster_globalR;
  std::vector<std::vector<double>>  cluster_charge;     
       
  TTree *tree;
       
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HitAnalyzer::HitAnalyzer(const edm::ParameterSet& conf)
: conf_(conf), src_(conf.getParameter<edm::InputTag>( "src" )) { 
  
  printLocal = conf.getUntrackedParameter<bool>("Verbosity",false);
  phase1_ = conf.getUntrackedParameter<bool>("phase1",false);
        
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>( "tree", "tree" );
  
  tree->Branch( "nJets"             , &nJets );
  tree->Branch( "jet_pt"            , &jet_pt );
  tree->Branch( "jet_eta"           , &jet_eta );
  tree->Branch( "jet_phi"           , &jet_phi );
  tree->Branch( "jet_mass"          , &jet_mass );
  // tree->Branch( "jet_pdgId"         , &jet_pdgId );
  tree->Branch( "jet_bTag"          , &jet_bTag );
  // tree->Branch( "jet_vx_x"  , &jet_vx_x   );
  // tree->Branch( "jet_vx_y"  , &jet_vx_y   );
  // tree->Branch( "jet_vx_z"  , &jet_vx_z   );
  // tree->Branch( "jet_vx_eta", &jet_vx_eta );
  // tree->Branch( "jet_vx_phi", &jet_vx_phi );
  // tree->Branch( "jet_vx_r"  , &jet_vx_r   );
  
  
  
  // tree->Branch( "dr_jetGen"         , &dr_jetGen );
  
  tree->Branch( "nGenParticles"     , &nGenParticles );
  tree->Branch( "genParticle_pt"    , &genParticle_pt );
  tree->Branch( "genParticle_eta"   , &genParticle_eta );
  tree->Branch( "genParticle_phi"   , &genParticle_phi );
  tree->Branch( "genParticle_mass"  , &genParticle_mass );
  tree->Branch( "genParticle_pdgId" , &genParticle_pdgId );
  tree->Branch( "genParticle_mother_pdgId" , &genParticle_mother_pdgId );
  tree->Branch( "genParticle_status", &genParticle_status );
  tree->Branch( "genParticle_vx_x"   , &genParticle_vx_x   );
  tree->Branch( "genParticle_vx_y"   , &genParticle_vx_y   );
  tree->Branch( "genParticle_vx_z"   , &genParticle_vx_z   );
  tree->Branch( "genParticle_vx_eta" , &genParticle_vx_eta );
  tree->Branch( "genParticle_vx_phi" , &genParticle_vx_phi );
  tree->Branch( "genParticle_vx_r"   , &genParticle_vx_r   );
  tree->Branch( "genParticle_decayvx_x"   , &genParticle_decayvx_x   );
  tree->Branch( "genParticle_decayvx_y"   , &genParticle_decayvx_y   );
  tree->Branch( "genParticle_decayvx_z"   , &genParticle_decayvx_z   );
  tree->Branch( "genParticle_decayvx_eta" , &genParticle_decayvx_eta );
  tree->Branch( "genParticle_decayvx_phi" , &genParticle_decayvx_phi );
  tree->Branch( "genParticle_decayvx_r"   , &genParticle_decayvx_r   );
  
  tree->Branch( "nTracks"      , &nTracks );
  tree->Branch( "track_pt"     , &track_pt );
  tree->Branch( "track_eta"    , &track_eta );
  tree->Branch( "track_phi"    , &track_phi );
  //tree->Branch( "track_mass"   , &track_mass );
  tree->Branch( "track_charge" , &track_charge)   ;
  tree->Branch( "track_px"     , &track_px)   ;
  tree->Branch( "track_py"     , &track_py)   ;
  tree->Branch( "track_pz"     , &track_pz) ;
  tree->Branch( "track_vx_x"   , &track_vx_x   );
  tree->Branch( "track_vx_y"   , &track_vx_y   );
  tree->Branch( "track_vx_z"   , &track_vx_z   );
  //tree->Branch( "track_vx_eta" , &track_vx_eta );
  tree->Branch( "track_vx_phi" , &track_vx_phi );
  tree->Branch( "track_vx_r"   , &track_vx_r   );
  
  
  
  
  
  
  
  
  
  tree->Branch( "nDetUnits"         ,&nDetUnits           );
  // tree->Branch( "detUnit_detType"   ,&detUnit_detType     );
  tree->Branch( "detUnit_subdetId"  ,&detUnit_subdetId    );
  // tree->Branch( "detUnit_layerC"    ,&detUnit_layerC      );
  // tree->Branch( "detUnit_ladderC"   ,&detUnit_ladderC     );
  // tree->Branch( "detUnit_zindex"    ,&detUnit_zindex      );
  // tree->Branch( "detUnit_shell"     ,&detUnit_shell       );
  // tree->Branch( "detUnit_sector"    ,&detUnit_sector      );
  // tree->Branch( "detUnit_ladder"    ,&detUnit_ladder      );
  tree->Branch( "detUnit_layer"     ,&detUnit_layer       ); 
  // tree->Branch( "detUnit_module"    ,&detUnit_module      );
  // tree->Branch( "detUnit_half"      ,&detUnit_half        );
  tree->Branch( "detUnit_disk"      ,&detUnit_disk        );   
  // tree->Branch( "detUnit_blade"     ,&detUnit_blade       );
  // tree->Branch( "detUnit_moduleF"   ,&detUnit_moduleF     );
  tree->Branch( "detUnit_side"      ,&detUnit_side        );
  // tree->Branch( "detUnit_panel"     ,&detUnit_panel       );
  
  tree->Branch( "detUnit_X"   ,&detUnit_X                   );
  tree->Branch( "detUnit_Y"   ,&detUnit_Y                   );
  tree->Branch( "detUnit_Z"   ,&detUnit_Z                   );
  tree->Branch( "detUnit_R"   ,&detUnit_R                   );
  tree->Branch( "detUnit_Phi" ,&detUnit_Phi                 );
  
  
  
  tree->Branch( "nClusters"         ,&nClusters           );
  tree->Branch( "cluster_sizex"     ,&cluster_sizex           );
  tree->Branch( "cluster_sizey"     ,&cluster_sizey           );
  tree->Branch( "cluster_globalz"   ,&cluster_globalz     );
  tree->Branch( "cluster_globalx"   ,&cluster_globalx     );
  tree->Branch( "cluster_globaly"   ,&cluster_globaly     );
  tree->Branch( "cluster_globalPhi" ,&cluster_globalPhi   );
  tree->Branch( "cluster_globalR"   ,&cluster_globalR     );
  tree->Branch( "cluster_charge"    ,&cluster_charge	  );
  
    
   
  std::string labelgenP("genParticles");
  std::string labelAK8s("ak8PFJetsCHS");
  std::string labelAK4s("ak4PFJetsCHS");
  std::string labelClusters("siPixelClusters");
  std::string labelSVs("inclusiveSecondaryVertices");
  std::string labelPVs("offlinePrimaryVertices");
  std::string labelTracks("generalTracks");
  std::string labelCSV("pfCombinedSecondaryVertexV2BJetTags");
  genPtoken      = consumes<reco::GenParticleCollection         > (edm::InputTag(labelgenP));
  ak8CHStoken    = consumes<reco::PFJetCollection               > (edm::InputTag(labelAK8s));
  ak4CHStoken    = consumes<reco::PFJetCollection               > (edm::InputTag(labelAK4s));
  clusterToken   = consumes<edmNew::DetSetVector<SiPixelCluster>> (src_);
  svToken        = consumes<reco::VertexCollection              > (edm::InputTag(labelSVs));
  pvToken        = consumes<reco::VertexCollection              > (edm::InputTag(labelPVs));
  csv2Token      = consumes<reco::JetTagCollection              > (edm::InputTag(labelCSV));
  trackToken     = consumes< reco::TrackCollection               >(edm::InputTag(labelTracks));

}


HitAnalyzer::~HitAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
  HitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  using namespace edm;
   
  // Make sure vectors are clear
  reset();
   
  // // Get event setup
  edm::ESHandle<TrackerGeometry> geom;
  iSetup.get<TrackerDigiGeometryRecord>().get( geom );
  const TrackerGeometry& theTracker(*geom);
  
#ifdef NEW_ID
  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopoH;
  iSetup.get<TrackerTopologyRcd>().get(tTopoH);
  const TrackerTopology *tTopo=tTopoH.product();
#endif
  
   
  // Get handles
  Handle<edmNew::DetSetVector<SiPixelCluster> > clusters ; iEvent.getByToken( clusterToken , clusters );
  Handle<reco::PFJetCollection                > ak8CHS   ; iEvent.getByToken( ak8CHStoken  , ak8CHS   );
  Handle<reco::PFJetCollection                > ak4CHS   ; iEvent.getByToken( ak4CHStoken  , ak4CHS   );
  Handle<reco::GenParticleCollection          > genPs    ; iEvent.getByToken( genPtoken    , genPs    );
  Handle<reco::JetTagCollection               > CSVs     ; iEvent.getByToken( csv2Token    , CSVs     );
  Handle<reco::VertexCollection               > SVs      ; iEvent.getByToken( svToken      , SVs      );
  Handle<reco::TrackCollection                > tracks   ; iEvent.getByToken( trackToken   , tracks   );
  const reco::JetTagCollection & bTags = *(CSVs.product()); 
  
 
  // Loop over jets
  for ( reco::PFJetCollection::const_iterator jet = ak4CHS->begin(); jet != ak4CHS->end(); ++jet ) {
    if (jet->pt()<200.) continue;
    TLorentzVector TVjet;
    TVjet.SetPtEtaPhiM(jet->pt(),jet->eta(),jet->phi(),jet->mass());
    jet_pt.push_back(jet->pt());
    jet_eta.push_back(jet->eta());
    jet_phi.push_back(jet->phi());
    jet_mass.push_back(jet->mass());
    
    // b-tag infos
    double match = 0.4;
    double csv2 = -99.;
    for (unsigned int i = 0; i != bTags.size(); ++i) {
      if (bTags[i].first->pt()<170.) continue;
      TLorentzVector bTagJet;
      bTagJet.SetPtEtaPhiM(bTags[i].first->pt(),bTags[i].first->eta(),bTags[i].first->phi(),bTags[i].first->mass());
      float dR = TVjet.DeltaR(bTagJet);
      if (dR > match ) continue;
      match = dR;
      csv2 = bTags[i].second; 
    }
    jet_bTag.push_back(csv2);
  }
  nJets         = jet_pt.size();
  // jet_pdgId.push_back(pdgId);
  // jet_vx_x   .push_back(jet->vertex().Coordinates().X());
  // jet_vx_y   .push_back(jet->vertex().Coordinates().Y());
  // jet_vx_z   .push_back(jet->vertex().Coordinates().Z());
  // jet_vx_eta .push_back(jet->vertex().Coordinates().Eta());
  // jet_vx_phi .push_back(jet->vertex().Coordinates().Phi());
  // jet_vx_r   .push_back(jet->vertex().Coordinates().R());
    

    
    // Loop opver gen particles
    TLorentzVector selectedGenP ;
    std::cout<< " " << std::endl;
    for ( reco::GenParticleCollection::const_iterator gp = genPs->begin(); gp != genPs->end(); ++gp ) {
      bool BMother = false;
      int motherID = false;
      for( unsigned int m=0; m<gp->numberOfMothers(); ++m ){
        motherID =  gp->mother(m)->pdgId();
        if ( (fabs(gp->mother(m)->pdgId())>500 && fabs(gp->mother(m)->pdgId())<600) || (fabs(gp->mother(m)->pdgId())>5000 && fabs(gp->mother(m)->pdgId())<6000) ) {
          BMother = true;
	  if (BMother) break;
        }
      }
      
      //if (BMother && (fabs(gp->pdgId())>500 && fabs(gp->pdgId())<600)) continue;

      std::string idString = std::to_string(fabs(gp->pdgId()));
      /*
      if ( BMother == true || (gp->status()<30 && gp->status()>20) || (fabs(gp->pdgId())>500 && fabs(gp->pdgId())<600) 
      || (idString.find("511") != std::string::npos)  || (idString.find("521") != std::string::npos)
      || (idString.find("513") != std::string::npos)  || (idString.find("523") != std::string::npos)
      || (idString.find("531") != std::string::npos)  || (idString.find("541") != std::string::npos)
      || (idString.find("533") != std::string::npos)  || (idString.find("543") != std::string::npos) 
      || (idString.find("551") != std::string::npos)  || (idString.find("553") != std::string::npos) 
      || (idString.find("555") != std::string::npos)  || (idString.find("557") != std::string::npos) ){ */ 
      if (true){
	std::vector<const reco::Candidate*> genParticle_Daughters;
	double decayvx_x = 0;
	double decayvx_y = 0;
	double decayvx_z = 0;
	double decayvx_r = 0;
	double decayvx_phi = 0;
	double decayvx_eta = 0;

	bool IdenticalDaughter = false;
	unsigned int dau_identical = 0;
	for (unsigned int dau=0; dau < gp->numberOfDaughters(); ++dau) {	//check if daughter identical to particle
		if ((gp->daughter(dau)->pdgId() == gp->pdgId()) || ((fabs(gp->pdgId())>500 && fabs(gp->pdgId())<600) && (fabs(gp->daughter(dau)->pdgId())>500 && fabs(gp->daughter(dau)->pdgId())<600)) || ((fabs(gp->pdgId())>5000 && fabs(gp->pdgId())<6000) && (fabs(gp->daughter(dau)->pdgId())>5000 && fabs(gp->daughter(dau)->pdgId())<6000))) {
			IdenticalDaughter = true;
			dau_identical = dau;
			break;
			}
		}
	if (IdenticalDaughter){			//iterate over generations until daughter differs from particle
//		std::cout<<std::endl<<"particle id = "<<gp->pdgId()<<std::endl;
//		std::cout<<"Generation = "<<0<<std::endl;
//		std::cout<<"identical daughter, id = "<<gp->daughter(0)->pdgId()<<std::endl; 	
		int GenerationCounter = 1;
		const reco::Candidate *tmp; 
		tmp = gp->daughter(dau_identical);
		LOOP:
//		std::cout<<"Generation = "<<GenerationCounter<<std::endl;
		IdenticalDaughter = false;
		for (unsigned int dau=0; dau < tmp->numberOfDaughters(); ++dau) {
			if ((tmp->pdgId() == tmp->daughter(dau)->pdgId()) || ((fabs(tmp->pdgId())>500 && fabs(tmp->pdgId())<600) && (fabs(tmp->daughter(dau)->pdgId())>500 && fabs(tmp->daughter(dau)->pdgId())<600)) || ((fabs(tmp->pdgId())>5000 && fabs(tmp->pdgId())<6000) && (fabs(tmp->daughter(dau)->pdgId())>5000 && fabs(tmp->daughter(dau)->pdgId())<6000))) {
				tmp = tmp->daughter(dau);
				dau = 0;
				IdenticalDaughter = true;
				++GenerationCounter;
//				std::cout<<"identical daughter, id = "<<tmp->pdgId()<<std::endl;
				break;
				}		
			}
		if (IdenticalDaughter){
			goto LOOP;
			}
		else {
//			std::cout<<"Final daugthers:"<<std::endl;
			for (unsigned int dau=0; dau < tmp->numberOfDaughters(); ++dau) {
				genParticle_Daughters.push_back(tmp->daughter(dau));
//				std::cout<<"d"<<dau<<"_id = "<<tmp->daughter(dau)->pdgId()<<std::endl;	
				}
			}
		}
	else {
		for (unsigned int dau=0; dau < gp->numberOfDaughters(); ++dau) {
			genParticle_Daughters.push_back(gp->daughter(dau));			
			}
		}
	for (unsigned int dau=0; dau < genParticle_Daughters.size(); ++dau) {	//loop over the actual daughters
		if (dau == 0){
		decayvx_x = genParticle_Daughters[dau]->vertex().x();
		decayvx_y = genParticle_Daughters[dau]->vertex().y();
		decayvx_z = genParticle_Daughters[dau]->vertex().z();
		decayvx_r = genParticle_Daughters[dau]->vertex().r();
		decayvx_phi = genParticle_Daughters[dau]->vertex().phi();
		decayvx_eta = genParticle_Daughters[dau]->vertex().eta();
		}
		break;								//break after first daughter since their vertices are the same
		}
	

        TLorentzVector genP;
        genP.SetPtEtaPhiM(gp->pt(),gp->eta(),gp->phi(),gp->mass());
        int  pdgId = gp->pdgId();
        int  status = gp->status();
        double vx_x   = gp->vertex().Coordinates().X(); 
        double vx_y   = gp->vertex().Coordinates().Y(); 
        double vx_z   = gp->vertex().Coordinates().Z(); 
        double vx_eta = gp->vertex().Coordinates().Eta(); 
        double vx_phi = gp->vertex().Coordinates().Phi(); 
        double vx_r   = gp->vertex().Coordinates().R(); 

        /*
	if ((fabs(gp->pdgId())>500 && fabs(gp->pdgId())<600) 
      || (idString.find("511") != std::string::npos)  || (idString.find("521") != std::string::npos)
      || (idString.find("513") != std::string::npos)  || (idString.find("523") != std::string::npos)
      || (idString.find("531") != std::string::npos)  || (idString.find("541") != std::string::npos)
      || (idString.find("533") != std::string::npos)  || (idString.find("543") != std::string::npos) 
      || (idString.find("551") != std::string::npos)  || (idString.find("553") != std::string::npos) 
      || (idString.find("555") != std::string::npos)  || (idString.find("557") != std::string::npos)){
	const double distance = std::sqrt((decayvx_x-vx_x)*(decayvx_x-vx_x)+(decayvx_y-vx_y)*(decayvx_y-vx_y)+(decayvx_z-vx_z)*(decayvx_z-vx_z));
	std::cout<<std::endl<<"particle pdgId() = "<<gp->pdgId()<<"\t pt = "<<gp->pt()<<std::endl<<"particle_vx = ("<<vx_x<<", "<<vx_y<<", "<<vx_z<<") \t decay_vx = ("<<decayvx_x<<", "<<decayvx_y<<", "<<decayvx_z<<")"<<std::endl<<"distance travelled = "<<distance<<std::endl;
	std::cout<<"daughters: ";
	for (unsigned int dau=0; dau < genParticle_Daughters.size(); ++dau) {
	std::cout<<"   d"<<dau<<"_pdgId() = "<<genParticle_Daughters[dau]->pdgId();
	}
	std::cout<<std::endl;
	}
	*/
        genParticle_pt  .push_back(genP.Pt());
        genParticle_eta .push_back(genP.Eta());
        genParticle_phi .push_back(genP.Phi());
        genParticle_mass.push_back(genP.M());
        genParticle_pdgId.push_back(pdgId);
        genParticle_mother_pdgId.push_back(motherID);
        genParticle_status.push_back(status);
        genParticle_vx_x   .push_back(vx_x   );
        genParticle_vx_y   .push_back(vx_y   );
        genParticle_vx_z   .push_back(vx_z   );
        genParticle_vx_eta .push_back(vx_eta );
        genParticle_vx_phi .push_back(vx_phi );
        genParticle_vx_r   .push_back(vx_r   );
	genParticle_decayvx_x   .push_back(decayvx_x   );
        genParticle_decayvx_y   .push_back(decayvx_y   );
        genParticle_decayvx_z   .push_back(decayvx_z   );
 	genParticle_decayvx_r   .push_back(decayvx_r   );
        genParticle_decayvx_eta   .push_back(decayvx_eta   );
        genParticle_decayvx_phi   .push_back(decayvx_phi   );
 

    }
    }
    

  nGenParticles = genParticle_pt.size();
  
  // Loop over tracks
  for ( reco::TrackCollection::const_iterator track = tracks->begin(); track != tracks->end(); ++track ) {
    TLorentzVector TVtrack;
    // TVtrack.SetPtEtaPhiM(track->pt(),track->eta(),track->phi(),track->mass());
    track_pt .push_back(track->pt());
    track_eta.push_back(track->eta());
    track_phi.push_back(track->phi()); 
    track_px.push_back(track->px());
    track_py.push_back(track->py());
    track_pz.push_back(track->pz());
    track_vx_x   .push_back(track->vertex().x() );
    track_vx_y   .push_back(track->vertex().y() );
    track_vx_z   .push_back(track->vertex().z() );
    //track_vx_eta .push_back(track->vertex().eta());
    track_vx_phi .push_back(track->vertex().phi());
    track_vx_r   .push_back(track->vertex().r() );
    track_charge .push_back(track->charge());
    // track_mass.push_back(track->mass());
  }
  nTracks = track_pt.size();

  // Get vector of detunit ids and loop
  const edmNew::DetSetVector<SiPixelCluster>& input = *clusters;
  
  int numberOfDetUnits = 0;
  
  for ( edmNew::DetSetVector<SiPixelCluster>::const_iterator detUnit = input.begin(); detUnit != input.end(); ++detUnit ) {
    if (detUnit->size()<1) continue;
    unsigned int detid = detUnit->detId();
    DetId detId = DetId(detid);       // Get the Detid object
    unsigned int detType=detId.det(); // det type, pixel=1
    if(detType!=1) continue; // look only at pixels
    unsigned int subid=detId.subdetId(); //subdetector type, pix barrel=1, forward=2
    // Subdet id, pix barrel=1, forward=2
    
    // Get the geom-detector
    const PixelGeomDetUnit * theGeomDet = dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDet(detId) );
    const PixelTopology * topol = &(theGeomDet->specificTopology());
    double detX = theGeomDet->surface().position().x();
    double detY = theGeomDet->surface().position().y();
    double detZ = theGeomDet->surface().position().z();
    double detR = theGeomDet->surface().position().perp();
    double detPhi = theGeomDet->surface().position().phi();
    // barrel ids
    // unsigned int layerC=0;
    // unsigned int ladderC=0;
    // unsigned int zindex=0;
    // int shell  = 0; // shell id // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
    // int sector = 0; // 1-8
    // int ladder = 0; // 1-22
    int layer  = 0; // 1-3
    // int module = 0; // 1-4
    // bool half  = false; //

    // Endcap ids
    unsigned int disk=0; //1,2,3
    // unsigned int blade=0; //1-24
    // unsigned int moduleF=0; // plaquette 1,2,3,4
    unsigned int side=0; //size=1 for -z, 2 for +z
    // unsigned int panel=0; //panel=1
   
   
    if(subid==2) {  // forward
#ifdef NEW_ID
      disk=tTopo->pxfDisk(detid); //1,2,3
      // blade=tTopo->pxfBlade(detid); //1-24
      // zindex=tTopo->pxfModule(detid); //
      side=tTopo->pxfSide(detid); //size=1 for -z, 2 for +z
      // panel=tTopo->pxfPanel(detid); //panel=1
      PixelEndcapName pen(detid,tTopo,phase1_);
#else 
      PXFDetId pdetId = PXFDetId(detid); 
      disk=pdetId.disk();      //1,2,3
      // blade=pdetId.blade();    //1-24
      // moduleF=pdetId.module(); // plaquette
      side=pdetId.side();      //size=1 for -z, 2 for +z
      // panel=pdetId.panel();    //panel=1
#endif
    }
    else if (subid==1) {  // barrel
#ifdef NEW_ID
      // layerC=tTopo->pxbLayer(detid);
      // ladderC=tTopo->pxbLadder(detid);
      // zindex=tTopo->pxbModule(detid);
      PixelBarrelName pbn(detid,tTopo,phase1_);
#else      
      PXBDetId pdetId = PXBDetId(detid);
     
      // layerC=pdetId.layer();       // Barell layer = 1,2,3
      // ladderC=pdetId.ladder();     // Barrel ladder id 1-20,32,44.
      // zindex=pdetId.module();      // Barrel Z-index=1,8
      PixelBarrelName pbn(pdetId); // Convert to online     
#endif

      // PixelBarrelName::Shell sh = pbn.shell(); //enum
      // sector = pbn.sectorName();
      // ladder = pbn.ladderName();
      layer  = pbn.layerName();
      // module = pbn.moduleName();
      // half  = pbn.isHalfModule();
      // shell = int(sh);
    
      // if(shell==1 || shell==2) module = -module; // change the module sign for z<0
      // if(shell==1 || shell==3) ladder = -ladder; // change ladeer sign for Outer )x<0)
    }
    

    
    std::vector<double>  _cluster_sizex;
    std::vector<double>  _cluster_sizey;
    std::vector<double>  _cluster_globalz;
    std::vector<double>  _cluster_globalx;
    std::vector<double>  _cluster_globaly;
    std::vector<double>  _cluster_globalPhi;
    std::vector<double>  _cluster_globalR;
    std::vector<double>  _cluster_charge;
    int numberOfClusters = 0;
    for ( edmNew::DetSet<SiPixelCluster>::const_iterator clustIt = detUnit->begin(); clustIt != detUnit->end(); ++clustIt ) {
      numberOfClusters++;

      // get global position of the cluster
      int sizeX = clustIt->sizeX(); //x=row=rfi, 
      int sizeY = clustIt->sizeY(); //y=col=z_global
      float x = clustIt->x(); // row, cluster position in pitch units, as float (int+0.5);
      float y = clustIt->y(); // column, analog average
      LocalPoint lp = topol->localPosition(MeasurementPoint(x,y));
      // float lx = lp.x(); // local cluster position in cm
      // float ly = lp.y();

      double charge = clustIt->charge();
	
      GlobalPoint clustgp = theGeomDet->surface().toGlobal( lp );
      double gZ = clustgp.z();  // global z
      double gX = clustgp.x();
      double gY = clustgp.y();
      TVector3 v(gX,gY,gZ);
      float gPhi = v.Phi(); // phi of the hit
      float gR = v.Perp(); // r of the hit
      
      _cluster_sizex.push_back(sizeX);
      _cluster_sizey.push_back(sizeY);
      _cluster_globalz.push_back(gZ);
      _cluster_globalx.push_back(gX);
      _cluster_globaly.push_back(gY);
      _cluster_globalPhi.push_back(gPhi);
      _cluster_globalR.push_back(gR);
      _cluster_charge.push_back(charge);
    }
    if( numberOfClusters < 1.) continue;
    
    // detUnit_detType.push_back(detType);
    detUnit_subdetId.push_back(subid);  
    // detUnit_layerC.push_back( layerC);
//     detUnit_ladderC.push_back(ladderC);
//     detUnit_zindex.push_back( zindex);
//     detUnit_shell .push_back( shell );
//     detUnit_sector.push_back( sector);
//     detUnit_ladder.push_back( ladder);
    detUnit_layer .push_back( layer ); 
    // detUnit_module.push_back( module);
  //   detUnit_half .push_back( half );
    detUnit_disk.push_back(disk);   
    // detUnit_blade.push_back(blade);
//     detUnit_moduleF.push_back(moduleF);
    detUnit_side.push_back(side);   
    // detUnit_panel.push_back(panel);
    detUnit_X         .push_back(detX     );
    detUnit_Y         .push_back(detY     );
    detUnit_Z         .push_back(detZ     );
    detUnit_R         .push_back(detR     );
    detUnit_Phi       .push_back(detPhi   );
    
    nClusters.push_back(numberOfClusters);
    cluster_sizex    .push_back(_cluster_sizex);
    cluster_sizey    .push_back(_cluster_sizey);
    cluster_globalz  .push_back(_cluster_globalz);
    cluster_globalx  .push_back(_cluster_globalx);
    cluster_globaly  .push_back(_cluster_globaly);
    cluster_globalPhi.push_back(_cluster_globalPhi);
    cluster_globalR  .push_back(_cluster_globalR);
    cluster_charge   .push_back(_cluster_charge);

    numberOfDetUnits++;
  }
  nDetUnits = numberOfDetUnits;
  
  tree->Fill();
}

// Private methods
void HitAnalyzer::reset( void ){

  nJets = 0;
  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_mass.clear();
  jet_pdgId.clear();
  jet_bTag.clear();
  // jet_vx_x  .clear();
  // jet_vx_y  .clear();
  // jet_vx_z  .clear();
  // jet_vx_eta.clear();
  // jet_vx_phi.clear();
  // jet_vx_r  .clear();
  dr_jetGen.clear();
  
  nGenParticles = 0;
  genParticle_pt.clear();
  genParticle_eta.clear();
  genParticle_phi.clear();
  genParticle_mass.clear();
  genParticle_pdgId.clear();
  genParticle_mother_pdgId.clear();
  genParticle_status.clear();
  genParticle_vx_x  .clear(); 
  genParticle_vx_y  .clear(); 
  genParticle_vx_z  .clear(); 
  genParticle_vx_eta.clear(); 
  genParticle_vx_phi.clear(); 
  genParticle_vx_r  .clear(); 
  genParticle_decayvx_x  .clear(); 
  genParticle_decayvx_y  .clear(); 
  genParticle_decayvx_z  .clear(); 
  genParticle_decayvx_eta.clear(); 
  genParticle_decayvx_phi.clear(); 
  genParticle_decayvx_r  .clear(); 
  
  nClusters.clear();
  cluster_sizex.clear();
  cluster_sizey.clear();
  cluster_globalz.clear();
  cluster_globalx.clear();
  cluster_globaly.clear();
  cluster_globalPhi.clear();
  cluster_globalR.clear();
  cluster_charge.clear();
  
  nDetUnits = 0.;
  // detUnit_detType.clear();
  detUnit_subdetId.clear();
  // detUnit_layerC.clear();
  // detUnit_ladderC.clear();
  // detUnit_zindex.clear();
  // detUnit_shell .clear();
  // detUnit_sector.clear();
  // detUnit_ladder.clear();
  detUnit_layer .clear();
  // detUnit_module.clear();
  // detUnit_half .clear();
  detUnit_disk.clear();
  // detUnit_blade.clear();
  // detUnit_moduleF.clear();
  detUnit_side.clear();
  // detUnit_panel.clear();
  detUnit_layer .clear();
  detUnit_X     .clear();
  detUnit_Y     .clear();
  detUnit_Z     .clear();
  detUnit_R     .clear();
  detUnit_Phi   .clear();
  
  
  nTracks = 0.;
  track_pt     .clear();
  track_eta    .clear();
  track_phi    .clear();
  //track_mass   .clear();
  track_px     .clear();
  track_py     .clear();
  track_pz     .clear();
  track_charge .clear();
  track_vx_x   .clear();
  track_vx_y   .clear();
  track_vx_z   .clear();
  //track_vx_eta .clear();
  track_vx_phi .clear();
  track_vx_r   .clear();
  
  
}

// Find mother
const reco::GenParticle *HitAnalyzer::findMother(const reco::GenParticle *particle) {
      const reco::GenParticle *tmp = particle;
      int pdgId = particle->pdgId();
      while(const reco::GenParticle *mother = dynamic_cast<const reco::GenParticle *>(tmp->mother())) {
        if(mother->pdgId() == pdgId) {
          tmp = mother;
          continue;
        }
        return mother;
      }
      return 0;
    }

// Find Daughters
/*
const std::vector<reco::Candidate> HitAnalyzer::findDaughters(const reco::GenParticle *particle) {
	std::vector<reco::Candidate> daughters;
	int pdgId = particle->pdgId();
	bool identical = false;
	unsigned int tmp;
	if (particle->numberOfDaughters() == 0) {
		return daughters;
		}
	for (unsigned int dau=0; dau < particle->numberOfDaughters(); ++dau ){
		if (particle->daughter(dau)->pdgId() == pdgId) {
			identical = true;
			tmp = dau;
			break;
			}
		}	
	if (identical) {
		return HitAnalyzer::findDaughters(particle->daughter(tmp));
		}	
	else {
		for( unsigned int dau=0; dau < particle->numberOfDaughters(); ++dau ){
			daughters.push_back(particle->daughter(dau));
			}
		return daughters;
		}
	}
*/
              
// ------------ method called once each job just before starting event loop  ------------
void 
  HitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
  HitAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HitAnalyzer);
