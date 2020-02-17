/*-*- c++ -*-*/

#ifndef MC_TRUTH_MATCH_SVC_H
#define MC_TRUTH_MATCH_SVC_H

#include "GaudiKernel/Service.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include <map>
#include <vector>
#include "McTruth/McParticle.h"
#include "McTruth/MdcMcHit.h"
#include "MdcRecEvent/RecMdcHit.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "MCTruthMatchSvc/IMCTruthMatchSvc.h"
#include "EventNavigator/EventNavigator.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "McDecayModeSvc/McDecayModeSvc.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
using std::vector;
using Event::McParticle;
using Event::McParticleCol;
using Event::MdcMcHitCol;
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif

class MC2TrkMap;
class McPar2ParMap;

typedef std::vector<const RecMdcHit*> RecMdcHitVector;

template <class TYPE> class CnvFactory;

bool comp(const SmartRef<Event::McParticle> , const SmartRef<Event::McParticle> );

class MCTruthMatchSvc : public Service, virtual public IMCTruthMatchSvc{
    
    friend class CnvFactory<MCTruthMatchSvc>;

    public :

    MCTruthMatchSvc(const std::string& name, ISvcLocator* svcLoc);
    virtual ~MCTruthMatchSvc();
    virtual StatusCode initialize();
    virtual StatusCode finalize();
    virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvIF);
    bool  matchTrack(const EvtRecTrack* recTrk, Event::McParticle* mcTrk);
    bool  matchTrack(const EvtRecTrack* recTrk,int pid, Event::McParticle* mcTrk);
    bool  match(vector<const EvtRecTrack*> recTracks,vector<int> Pids,vector<int> pdgid,vector<int> dechain);
    bool match(vector<const EvtRecTrack*> recTracks, vector<int> Pids, const McParticle *aEta);
    bool match(vector<const EvtRecTrack*> recTracks, vector<int> Pids, int motherPID);
    bool match(const EvtRecTrack* aTrack, int pid, int motherPID);
    private :
        bool m_init;
        bool m_OnlyCompareDirForGamma;
        double m_minValue ;
        IDataProviderSvc*  eventSvc_;
        McDecayModeSvc* m_svc;
        mutable EventNavigator* m_navigator;
        vector<McParticle*> m_pars;
        vector<int> m_noTracing, _essenPart;
        HepLorentzVector  getp4(const EvtRecTrack* aTrk, int pid);
        void FillFinalParticle();
        double findMin(HepLorentzVector);
        double findMinDir(HepLorentzVector);
        RecMdcHitVector getMdcRecHits(McParticle* mcTrk);
        bool irrational(HepLorentzVector p4);
        int m_run, m_event;
        void GetFinalChildren(McParticle *aEta, vector<McParticle*> &mcParticlaList);
};

#endif
