/////////////////////////////////////////////////////////////////////////
//                                                                     //
//  Analyzer for making mini-ntuple for HiggstoTauTau GEN level events //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

///////////////////////
// DATA FORMATS HEADERS
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "DataFormats/Common/interface/Ref.h"

//#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
//#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
//#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
//#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
//#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
//#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
//#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
//#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
//#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
//#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticle.h"
//#include "DataFormats/L1TrackTrigger/interface/L1TkEmParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//#include "SimDataFormats/SLHC/interface/L1EGCrystalCluster.h"
#include "DataFormats/Math/interface/LorentzVector.h"

////////////////
// PHYSICS TOOLS
#include "CommonTools/UtilAlgos/interface/TFileService.h"

///////////////
// ROOT HEADERS
#include <TROOT.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

//////////////
// STD HEADERS
#include <memory>
#include <string>
#include <iostream>

//////////////
// NAMESPACES
using namespace std;
using namespace edm;
using namespace reco;
using namespace math;

//////////////////////////////
//                          //
//     CLASS DEFINITION     //
//                          //
//////////////////////////////

class GeneratorNtupleMaker : public edm::EDAnalyzer
{
public:

  // Constructor/destructor
  explicit GeneratorNtupleMaker(const edm::ParameterSet& iConfig);
  virtual ~GeneratorNtupleMaker();

  // Mandatory methods
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void getStableDaughters(const reco::Candidate & p,std::vector<const reco::Candidate *>& stabledaughters);
  virtual void printDaughters(const reco::Candidate & p, int & index, int mother_index, bool details);
  virtual void getTauDecayVertex(const reco::Candidate & tau, double& decayvx, double& decayvy, double& decayvz);

protected:
  
private:
  
  // Containers of parameters passed by python configuration file
  edm::ParameterSet config; 
  
  int MyProcess;

  //-----------------------------------------------------------------------------------------------
  // tree & branches for mini-ntuple

  TTree* eventTree;

  // GEN Objects
  std::vector<float>* m_higgs_pt;
  std::vector<float>* m_higgs_eta;
  std::vector<float>* m_higgs_phi;
  //std::vector<float>* m_higgs_vx;
  //std::vector<float>* m_higgs_vy;
  std::vector<float>* m_higgs_vr;
  std::vector<float>* m_higgs_vz; 
  std::vector<float>* m_higgs_mass;
  std::vector<int>*   m_higgs_id;

  std::vector<float>* m_tau_pt;
  std::vector<float>* m_tau_eta;
  std::vector<float>* m_tau_phi;
  //std::vector<float>* m_tau_vx; 
  //std::vector<float>* m_tau_vy;
  std::vector<float>* m_tau_vr;
  std::vector<float>* m_tau_vz;
  std::vector<float>* m_tau_mass;
  std::vector<int>*   m_tau_id;

  std::vector<int>*   m_tau_nprongs;
  std::vector<int>*   m_tau_nprongsK;

  //std::vector<int>*   m_tau_ndaug;

  std::vector<int>*   m_taudecay_index;
  //std::vector<float>* m_taudecay_vx;
  //std::vector<float>* m_taudecay_vy;
  std::vector<float>* m_taudecay_vr;
  std::vector<float>* m_taudecay_vz;
  std::vector<int>*   m_taudecayprod_id;

  std::vector<int>*   m_stabledaug_index;  //come from which tau
  std::vector<int>*   m_tau_nvisstabledaug;
  std::vector<float>* m_tau_vispt;
  std::vector<float>* m_tau_visE;
  std::vector<float>* m_tau_visP;
  std::vector<float>* m_stabledaug_pt;
  std::vector<float>* m_stabledaug_eta;
  std::vector<float>* m_stabledaug_phi;
  std::vector<float>* m_stabledaug_mass;
  std::vector<int>  * m_stabledaug_id;

  std::vector<float>* m_recohiggs_mass;
};


//////////////////////////////////
//                              //
//     CLASS IMPLEMENTATION     //
//                              //
//////////////////////////////////

//////////////
// CONSTRUCTOR
GeneratorNtupleMaker::GeneratorNtupleMaker(edm::ParameterSet const& iConfig) : 
  config(iConfig)
{

  MyProcess = iConfig.getParameter< int >("MyProcess");
  
}

/////////////
// DESTRUCTOR
GeneratorNtupleMaker::~GeneratorNtupleMaker()
{
}  

//////////
// END JOB
void GeneratorNtupleMaker::endJob()
{
  // things to be done at the exit of the event Loop
  cerr << "GeneratorNtupleMaker::endJob" << endl;

}

////////////
// BEGIN JOB
void GeneratorNtupleMaker::beginJob()
{

  // things to be done before entering the event Loop
  cerr << "GeneratorNtupleMaker::beginJob" << endl;


  //-----------------------------------------------------------------------------------------------
  // book histograms / make ntuple
  edm::Service<TFileService> fs;


  // initilize
  m_higgs_pt   = new std::vector<float>;
  m_higgs_eta  = new std::vector<float>;
  m_higgs_phi  = new std::vector<float>;
  //m_higgs_vx   = new std::vector<float>;
  //m_higgs_vy   = new std::vector<float>;
  m_higgs_vr   = new std::vector<float>;
  m_higgs_vz   = new std::vector<float>;
  m_higgs_mass = new std::vector<float>;
  m_higgs_id   = new std::vector<int>;

  m_tau_pt   = new std::vector<float>;
  m_tau_eta  = new std::vector<float>;
  m_tau_phi  = new std::vector<float>;
  //m_tau_vx   = new std::vector<float>;
  //m_tau_vy   = new std::vector<float>;
  m_tau_vr   = new std::vector<float>;
  m_tau_vz   = new std::vector<float>;
  m_tau_mass = new std::vector<float>;
  m_tau_id   = new std::vector<int>;

  m_tau_nprongs = new std::vector<int>;
  m_tau_nprongsK = new std::vector<int>;
 
  //m_tau_ndaug = new std::vector<int>;
  m_taudecay_index = new std::vector<int>;
  //m_taudecay_vx = new std::vector<float>;
  //m_taudecay_vy = new std::vector<float>;
  m_taudecay_vr = new std::vector<float>;
  m_taudecay_vz = new std::vector<float>;
  m_taudecayprod_id = new std::vector<int>;

  m_stabledaug_index = new std::vector<int>;
  m_tau_nvisstabledaug = new std::vector<int>;
  m_tau_vispt = new std::vector<float>;
  m_tau_visE  = new std::vector<float>;
  m_tau_visP  = new std::vector<float>;
  m_stabledaug_pt = new std::vector<float>;
  m_stabledaug_eta = new std::vector<float>;
  m_stabledaug_phi = new std::vector<float>;
  m_stabledaug_mass = new std::vector<float>;
  m_stabledaug_id  = new std::vector<int>;

  m_recohiggs_mass = new std::vector<float>;
  
  // ntuple
  eventTree = fs->make<TTree>("H2TauTauTree", "Event tree");

  eventTree->Branch("higgs_pt",      &m_higgs_pt);
  eventTree->Branch("higgs_eta",     &m_higgs_eta);
  eventTree->Branch("higgs_phi",     &m_higgs_phi);
  //eventTree->Branch("higgs_vx",      &m_higgs_vx);
  //eventTree->Branch("higgs_vy",      &m_higgs_vy);
  eventTree->Branch("higgs_vr",      &m_higgs_vr);
  eventTree->Branch("higgs_vz",      &m_higgs_vz);
  eventTree->Branch("higgs_mass",    &m_higgs_mass);
  eventTree->Branch("higgs_id",      &m_higgs_id);

  eventTree->Branch("tau_pt",      &m_tau_pt);
  eventTree->Branch("tau_eta",     &m_tau_eta);
  eventTree->Branch("tau_phi",     &m_tau_phi);
  //eventTree->Branch("tau_vx",       &m_tau_vx);
  //eventTree->Branch("tau_vy",       &m_tau_vy);
  eventTree->Branch("tau_vr",      &m_tau_vr);
  eventTree->Branch("tau_vz",      &m_tau_vz);
  eventTree->Branch("tau_mass",    &m_tau_mass);
  eventTree->Branch("tau_id",      &m_tau_id);

  eventTree->Branch("tau_nprongs", &m_tau_nprongs);
  eventTree->Branch("tau_nprongsK",&m_tau_nprongsK);
  
  eventTree->Branch("taudecay_index",   &m_taudecay_index);
  //eventTree->Branch("taudecay_vx",   &m_taudecay_vx);
  //eventTree->Branch("taudecay_vy",   &m_taudecay_vy);
  eventTree->Branch("taudecay_vr",   &m_taudecay_vr);
  eventTree->Branch("taudecay_vz",   &m_taudecay_vz);
  eventTree->Branch("taudecayprod_id",  &m_taudecayprod_id);

  eventTree->Branch("stabledaug_index", &m_stabledaug_index);
  eventTree->Branch("tau_nvisstabledaug",  &m_tau_nvisstabledaug);
  eventTree->Branch("tau_vispt",        &m_tau_vispt);
  eventTree->Branch("tau_visE",         &m_tau_visE);
  eventTree->Branch("tau_visP",         &m_tau_visP);
  eventTree->Branch("stabledaug_pt",    &m_stabledaug_pt);
  eventTree->Branch("stabledaug_eta",   &m_stabledaug_eta);
  eventTree->Branch("stabledaug_phi",   &m_stabledaug_phi);
  eventTree->Branch("stabledaug_mass",  &m_stabledaug_mass);
  eventTree->Branch("stabledaug_id",    &m_stabledaug_id);

  eventTree->Branch("recohiggs_mass",   &m_recohiggs_mass);
}


//
// member functions
//
/*
int GeneratorNtupleMaker::tauClass(std::vector<const reco::Candidate *>& stabledaughters, double& maxpt) {

  // -999 means not classified
  // 1 means electron
  // 2 means muon
  // 3 means had 1 prong
  // 4 means had 3 prong
  // 5 means had 5 prong
  
  maxpt=-999.9;

  if (stabledaughters[1]->pdgId()==11||stabledaughters[1]->pdgId()==-11) return 1;
  if (stabledaughters[1]->pdgId()==13||stabledaughters[1]->pdgId()==-13) return 2;
    
  int nprong=0;  

  for (unsigned int i=0;i<stabledaughters.size();i++) {
    if (stabledaughters[i]->pdgId()==211||stabledaughters[i]->pdgId()==-211) {
      nprong++;
      if (stabledaughters[i]->pt()>maxpt) maxpt=stabledaughters[i]->pt();
    }
  }

  if (nprong==1) return 3;
  if (nprong==3) return 4;
  if (nprong==5) return 5;

  return -999;

}*/


void GeneratorNtupleMaker::getStableDaughters(const reco::Candidate & p,std::vector<const reco::Candidate *>& stabledaughters){

  int ndaug=p.numberOfDaughters();
  
  for(int j=0;j<ndaug;j++){
    const reco::Candidate * daug = p.daughter(j);
    if (daug->status()==1) {
      stabledaughters.push_back(daug);
    }
    else {
      getStableDaughters(*daug,stabledaughters);
    }
  }  
}

void GeneratorNtupleMaker::printDaughters(const reco::Candidate & p, int & index, int mother_index, bool details = false){

  int ndaug=p.numberOfDaughters();

  for(int j=0;j<ndaug;j++){
    const reco::Candidate * daug = p.daughter(j);
    index++;
    cout << index <<"\t"<<daug->pdgId() <<"\t"<< daug->status() <<"\t"<< mother_index <<"\t"<< daug->numberOfDaughters();
    if (details) {
      cout <<"\t"<<daug->pt()<<"\t"<<daug->eta()<<"\t"<<daug->phi()<<"\t"<<daug->vx()<<"\t"<<daug->vy()<<"\t"<<daug->vz()<<endl;
    }
    else cout<<endl;
    
    if (daug->status()!=1) printDaughters(*daug,index,index,details);
  }
  
}

void GeneratorNtupleMaker::getTauDecayVertex(const reco::Candidate & tau, double& decayvx, double& decayvy, double& decayvz) {  //Find the tau that does decay into other particles（instead of tau->tau->tau...). This tau should have status=2, and return its decay vertex
  //check if the input particle is tau
  if (fabs(tau.pdgId())!=15){
    cout << "not tau! return -999..."<< endl;
    decayvx = -999;
    decayvy = -999;
    decayvz = -999;
  }
  else {
    if (tau.status()==2) {
      const reco::Candidate * taudaug = tau.daughter(0);
      decayvx=taudaug->vx();
      decayvy=taudaug->vy();
      decayvz=taudaug->vz();
    } 
    else { //if not the one decay into other particles, loop over its daughters and find the tau it decays into
      int ndaug = tau.numberOfDaughters();
      for (int j=0;j<ndaug;j++){
	const reco::Candidate * nexttau = tau.daughter(j);
	if (fabs(nexttau->pdgId())!=15) continue;
	getTauDecayVertex(*nexttau,decayvx,decayvy,decayvz);
      }
    }

  }

}

//////////
// ANALYZE
void GeneratorNtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // clear variables
  m_higgs_pt->clear();
  m_higgs_eta->clear();
  m_higgs_phi->clear();
  //m_higgs_vx->clear();
  //m_higgs_vy->clear();
  m_higgs_vr->clear();
  m_higgs_vz->clear();
  m_higgs_mass->clear();
  m_higgs_id->clear();

  m_tau_pt->clear();
  m_tau_eta->clear();
  m_tau_phi->clear();
  //m_tau_vx->clear();
  //m_tau_vy->clear();
  m_tau_vr->clear();
  m_tau_vz->clear();
  m_tau_mass->clear();
  m_tau_id->clear();

  m_tau_nprongs->clear();
  m_tau_nprongsK->clear();

  //m_tau_ndaug->clear();
  m_taudecay_index->clear();
  //m_taudecay_vx->clear();
  //m_taudecay_vy->clear();
  m_taudecay_vr->clear();
  m_taudecay_vz->clear();
  m_taudecayprod_id->clear();

  m_stabledaug_index->clear();
  m_tau_nvisstabledaug->clear();
  m_tau_vispt->clear();
  m_tau_visE ->clear();
  m_tau_visP ->clear();
  m_stabledaug_pt->clear();
  m_stabledaug_eta->clear();
  m_stabledaug_phi->clear();
  m_stabledaug_mass->clear();
  m_stabledaug_id->clear();

  m_recohiggs_mass->clear();

  //-----------------------------------------------------------------------------------------------
  // retrieve various containers
  //-----------------------------------------------------------------------------------------------
  
  // MC truth association maps
  //edm::Handle< TTTrackAssociationMap< Ref_PixelDigi_ > > MCTruthTTTrackHandle;
  //iEvent.getByLabel("TTTrackAssociatorFromPixelDigis", "Level1TTTracks", MCTruthTTTrackHandle);
  //edm::Handle< TTClusterAssociationMap< Ref_PixelDigi_ > > MCTruthTTClusterHandle;
  //iEvent.getByLabel("TTClusterAssociatorFromPixelDigis", "ClusterAccepted", MCTruthTTClusterHandle);
  //edm::Handle< TTStubAssociationMap< Ref_PixelDigi_ > > MCTruthTTStubHandle;
  //iEvent.getByLabel("TTStubAssociatorFromPixelDigis", "StubAccepted", MCTruthTTStubHandle);

  // tracking particles
  //edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  //edm::Handle< std::vector< TrackingVertex > > TrackingVertexHandle;
  //iEvent.getByLabel("mix", "MergedTrackTruth", TrackingParticleHandle);
  //iEvent.getByLabel("mix", "MergedTrackTruth", TrackingVertexHandle);


  //Get generator MC truth

  Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);

  int nhiggs=0;

  //cout << "gen size: " << genParticles->size() << endl;

  for(size_t i = 0; i < genParticles->size(); ++ i) {
    const GenParticle & p = (*genParticles)[i];
    int id = p.pdgId();
    //int st = p.status();
    if (id!=25) continue;
    double tmp_eta=p.p4().eta();
    double tmp_phi=p.p4().phi();
    double tmp_pt=p.p4().pt();
    double tmp_vx = p.vx();
    double tmp_vy = p.vy();
    double tmp_vr = sqrt(tmp_vx*tmp_vx+tmp_vy*tmp_vy);
    double tmp_vz = p.vz();
    double tmp_mass = p.mass();

    if (p.daughter(0)->pdgId()==25) continue;  //skip if daughter is higgs
    m_higgs_pt  -> push_back(tmp_pt);
    m_higgs_eta -> push_back(tmp_eta);
    m_higgs_phi -> push_back(tmp_phi);
    m_higgs_vr  -> push_back(tmp_vr);
    m_higgs_vz  -> push_back(tmp_vz);
    m_higgs_mass -> push_back(tmp_mass);
    m_higgs_id  -> push_back(id);


    int ndaug = p.numberOfDaughters();

    PtEtaPhiMLorentzVector p4Tau[2];

    for(int idaug=0;idaug<ndaug;idaug++) {
      const reco::Candidate * daug = p.daughter(idaug);
      if (daug->pdgId()==25) continue;
      m_tau_pt    -> push_back(daug->pt());
      m_tau_eta   -> push_back(daug->eta());
      m_tau_phi   -> push_back(daug->phi());
      m_tau_id    -> push_back(daug->pdgId());
      m_tau_mass  -> push_back(daug->mass());
      double tmp_tauvx = daug->vx();
      double tmp_tauvy = daug->vy();
      m_tau_vr    -> push_back(sqrt(tmp_tauvx*tmp_tauvx+tmp_tauvy*tmp_tauvy));
      m_tau_vz    -> push_back(daug->vz());

      //loop over daugthers of Tau

      //Decay vertex
      double tmp_vx;
      double tmp_vy;
      double tmp_vz;
      getTauDecayVertex(*daug,tmp_vx,tmp_vy,tmp_vz);

      m_taudecay_vr -> push_back(sqrt(tmp_vx*tmp_vx+tmp_vy*tmp_vy));
      m_taudecay_vz -> push_back(tmp_vz);
      m_taudecay_index  -> push_back(idaug);  //from which tau


      //To check if the taus have the same kinematics and decay vertices. Also print out their status
      bool tau2tau = false;
      bool printTauDaugs = false;  //turn on to print list of tau daughters for events containing tau->tau 
      int ntaudaug = daug->numberOfDaughters();
      for (int itaudaug=0;itaudaug<ntaudaug;itaudaug++) {
	const reco::Candidate * grddaug = daug->daughter(itaudaug);
	int grd_id = grddaug->pdgId();
	if (grd_id == 15 || grd_id == -15) {
	  tau2tau = true;
	  continue;
	}

      }
      
      if (tau2tau && printTauDaugs) {
	cout << "**********************************************************" << endl;
	cout << "Find Tau -> Tau  event!" << endl;
	cout << "List of Tau daughters: " << endl;
	cout << "__________________________________________________________" << endl;
	cout << "index"<<"\t"<<"pdgid"<<"\t"<<"status"<<"\t"<<"mother"<<"\t"<<"ndaughters"<<"\t"<<"pt"<<"\t"<<"eta"<<"\t"<<"phi"<<"\t"<<"vertex"<<endl;
	int index = 0;
	printDaughters(*daug, index, 0, true);
	
	cout << "**********************************************************" << endl;
      }


      //Tau stable daughters
      std::vector<const reco::Candidate *> stabledaughters;

      getStableDaughters(*daug,stabledaughters);

      //m_tau_nvisstabledaug -> push_back(stabledaughters.size());
      int nvisdaug = 0;
      int nprongs = 0;
      int nprongsK = 0;
      bool foundPi0 = false;
      bool printPi0Evts = true; //turn on to print the events with Pi0 final states

      for (unsigned int j=0; j<stabledaughters.size();j++) {
	m_stabledaug_index -> push_back(idaug);  //from which tau
	m_stabledaug_pt    -> push_back(stabledaughters[j]->pt());
	m_stabledaug_eta   -> push_back(stabledaughters[j]->eta());
	m_stabledaug_phi   -> push_back(stabledaughters[j]->phi());
	m_stabledaug_mass  -> push_back(stabledaughters[j]->mass());
	m_stabledaug_id    -> push_back(stabledaughters[j]->pdgId());

	//number of stable daughters and tau's visible pt
	if (fabs(stabledaughters[j]->pdgId())!=12 && fabs(stabledaughters[j]->pdgId())!=14 && fabs(stabledaughters[j]->pdgId())!=16) {
	  p4Tau[idaug] = p4Tau[idaug]+PtEtaPhiMLorentzVector(stabledaughters[j]->pt(),stabledaughters[j]->eta(),stabledaughters[j]->phi(),stabledaughters[j]->mass());
	  nvisdaug++;
	}
	//number of prongs
	if (stabledaughters[j]->pdgId()==211 || stabledaughters[j]->pdgId()==-211) nprongs++;
	if (fabs(stabledaughters[j]->pdgId())==211 || fabs(stabledaughters[j]->pdgId())==321 ) nprongsK++;
	//pion0
	if (stabledaughters[j]->pdgId()==111) foundPi0 = true;
      }
      
      m_tau_nprongs -> push_back(nprongs);
      m_tau_nprongsK -> push_back(nprongsK);
      m_tau_nvisstabledaug -> push_back(nvisdaug);
      m_tau_vispt -> push_back(p4Tau[idaug].Pt());
      m_tau_visE  -> push_back(p4Tau[idaug].E());
      m_tau_visP  -> push_back(p4Tau[idaug].P());
      
      if (foundPi0 && printPi0Evts) {
	cout << "**********************************************************" << endl;
	cout << "Find an pi0 event!" << endl;
	cout << "List of Tau daughters: " << endl;
	cout << "__________________________________________________________" << endl;
	cout << "index"<<"\t"<<"pdgid"<<"\t"<<"status"<<"\t"<<"mother"<<"\t"<<"ndaughters"<<endl;
	int index = 0;
	printDaughters(*daug, index, 0, false);
	
	cout << "**********************************************************" << endl;
      }
      
    }  //end of tau loop

    m_recohiggs_mass -> push_back((p4Tau[0]+p4Tau[1]).M());
    //cout << "number of daughters = " << ndaug << endl;

    nhiggs++;
  } //end of gen higgs loop

  eventTree->Fill();

} //end of analyze()


///////////////////////////
// DEFINE THIS AS A PLUG-IN
DEFINE_FWK_MODULE(GeneratorNtupleMaker);
