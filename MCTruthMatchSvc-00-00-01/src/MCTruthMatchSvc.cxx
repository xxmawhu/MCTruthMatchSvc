#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "EventModel/EventHeader.h"
#include "MCTruthMatchSvc/MCTruthMatchSvc.h"
#include <algorithm>

using CLHEP::HepLorentzVector;
using namespace Event;
using std::cout;
using std::endl;
MCTruthMatchSvc::MCTruthMatchSvc(const std::string& name, ISvcLocator* svcLoc)
    : Service(name, svcLoc) {
    declareProperty("NoTracingList", m_noTracing);
    declareProperty("OnlyCompareDirForGamma", m_OnlyCompareDirForGamma = true);
    declareProperty("MinValue", m_minValue = 0.9);
    m_run = 0;
    m_event = 0;
    _essenPart.push_back(211);
    _essenPart.push_back(321);
    _essenPart.push_back(2212);
    _essenPart.push_back(11);
    _essenPart.push_back(13);
    _essenPart.push_back(22);
}
void MCTruthMatchSvc::FillFinalParticle() {
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc_,
                                                 "/Event/EventHeader");
    if (!eventHeader) {
        return;
    }
    if (m_run == eventHeader->runNumber() &&
        m_event == eventHeader->eventNumber()) {
        return;
    } else {
        m_run = eventHeader->runNumber();
        m_event = eventHeader->eventNumber();
    }
    m_pars.erase(m_pars.begin(), m_pars.end());

    SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc_,
                                                     "/Event/MC/McParticleCol");
    if (!mcParticleCol) {
        cout << "Could not retireve McparticleCol" << endl;
        return;
    }
    McParticleCol::iterator iter = mcParticleCol->begin();
    for (; iter != mcParticleCol->end(); iter++) {
        if (!(*iter)->decayFromGenerator()) continue;
        int m_pid = abs((*iter)->particleProperty());
        if (std::find(_essenPart.begin(), _essenPart.end(), m_pid) ==
            _essenPart.end()) {
            continue;
        }
        HepLorentzVector p4Truth = (*iter)->initialFourMomentum();
        if (irrational(p4Truth)) {
            continue;
        }
        m_pars.push_back(const_cast<McParticle*>(*iter));
    }
}
bool MCTruthMatchSvc::irrational(HepLorentzVector p4) {
    for (int i = 0; i < 4; i++) {
        if (isnan(p4[i])) return true;
    }
    if (p4.rho() > 10) return true;
    if (p4.rho() < 0.0001) return true;
    if (p4[3] < 0.0001) return true;
    if (p4[3] > 10) return true;
    return false;
}
bool MCTruthMatchSvc::matchTrack(const EvtRecTrack* recTrk, int pid,
                                 McParticle* mcTrk) {
    this->FillFinalParticle();

    int m_pid = fabs(mcTrk->particleProperty());

    if (m_pid != pid) {
        return false;
    }

    if (m_pid == 22) {
        if (!const_cast<EvtRecTrack*>(recTrk)->isEmcShowerValid()) {
            return false;
        }
    } else if (!const_cast<EvtRecTrack*>(recTrk)->isMdcKalTrackValid()) {
        return false;
    }

    HepLorentzVector p4Rec = getp4(recTrk, m_pid);
    HepLorentzVector p4Truth = mcTrk->initialFourMomentum();

    // for same events is truth information is not aviable,
    // and the p4 is not a number
    if (irrational(p4Truth)) {
        return false;
    }

    double value = p4Rec.howNear(p4Truth);
    double valueInEvent = this->findMin(p4Rec);
    if (m_OnlyCompareDirForGamma && m_pid == 22) {
        HepLorentzVector p4RecUnit = p4Rec / p4Rec.rho();
        HepLorentzVector p4TruthUnit = p4Truth / p4Truth.rho();
        value = p4RecUnit.howNear(p4TruthUnit);
        valueInEvent = findMinDir(p4Rec);
    }

    // for same track the hownear return 0.0, but the match is false
    if (fabs(value - valueInEvent) > 0.00001 || value > m_minValue) {
        return false;
    }
    return true;
    SmartDataPtr<EventNavigator> navigator(eventSvc_, "/Event/Navigator");
    if (!navigator) {
        cout << "Can't retrieve the Navigator!" << endl;
        return true;
    }
    // m_navigator = navigator;
    RecMdcTrack* recMdcTrk = const_cast<EvtRecTrack*>(recTrk)->mdcTrack();
    int nRecHits = recMdcTrk->getVecHits().size();
    cout << "MCTruthMatchSvc:: nRecHits " << nRecHits << endl;
    return true;
    std::vector<const McParticle*>& mcVec =
        navigator->getMcParticles(recMdcTrk);
    cout << " McParticleCol size: " << mcVec.size() << endl;
    if (mcVec.empty()) {
        cout << "MCTruthMatchSvc:: mcVec is empty!!!" << endl;
        // return true;
    }
    std::vector<const McParticle*>::iterator mcIt = mcVec.begin();
    int mcIndex = (mcVec.size() != 0) ? 0 : 1;
    cout << " mcTrk trackId: " << mcTrk->trackIndex() << endl;
    cout << " mcTrk pid: " << mcTrk->particleProperty() << endl;
    //   while (mcIt != mcVec.end()) {
    //       cout<<"ID: "<< (*mcIt)->trackIndex()<<endl;
    //       if((*mcIt)->trackIndex() == mcTrk->trackIndex()) break;
    //       if ((*mcIt) == mcTrk) break;
    //       mcIt++;
    //       mcIndex++;
    //   }
    //   if(mcIt == mcVec.end()){
    //       cout<<"MCTruthMatchSvc:: can't find the MC track"
    //           <<endl;
    //       return false;
    //   }

    int recTrkId = recTrk->trackId();

    int nMatchedHits = 0;
    std::vector<const RecMdcHit*> hits = getMdcRecHits(mcTrk);
    cout << "hits.size( ) = " << hits.size() << endl;
    for (unsigned int i = 0; i < hits.size(); i++) {
        if (hits[i]->getTrkId() == recTrkId) nMatchedHits++;
    }
    cout << "nMatchedHits  " << nMatchedHits << endl;
    cout << "nRecHits " << nRecHits << endl;

    return true;
}

bool MCTruthMatchSvc::matchTrack(const EvtRecTrack* recTrk, McParticle* mcTrk) {
    int pid = abs(mcTrk->particleProperty());
    if (std::find(_essenPart.begin(), _essenPart.end(), pid) ==
        _essenPart.end()) {
        return false;
    }
    return this->matchTrack(recTrk, pid, mcTrk);
}

double MCTruthMatchSvc::findMin(HepLorentzVector p4Rec) {
    double value = 10.;

    for (int i = 0; i < m_pars.size(); i++) {
        HepLorentzVector p4Truth = m_pars[i]->initialFourMomentum();
        double tmp = p4Rec.howNear(p4Truth);
        if (tmp < value && tmp > 0.00001) value = tmp;
    }
    return value;
}
double MCTruthMatchSvc::findMinDir(HepLorentzVector p4Rec) {
    double value = 10.;

    for (int i = 0; i < m_pars.size(); i++) {
        HepLorentzVector p4Truth = m_pars[i]->initialFourMomentum();

        HepLorentzVector p4TruthUnit = p4Truth / p4Truth.rho();
        double tmp = (p4Rec / p4Rec.rho()).howNear(p4TruthUnit);
        if (tmp < value && tmp > 0.00001) value = tmp;
    }
    return value;
}

HepLorentzVector MCTruthMatchSvc::getp4(const EvtRecTrack* aTrk, int pid) {
    int ID = fabs(pid);

    if (ID == 22) {
        RecEmcShower* emcShw = const_cast<EvtRecTrack*>(aTrk)->emcShower();
        double energy = emcShw->energy();
        double theta = emcShw->theta();
        double phi = emcShw->phi();
        HepLorentzVector p4Rec(energy * sin(theta) * cos(phi),
                               energy * sin(theta) * sin(phi),
                               energy * cos(theta), energy);
        return p4Rec;
    } else {
        RecMdcKalTrack* mdcKalTrack =
            const_cast<EvtRecTrack*>(aTrk)->mdcKalTrack();
        HepVector zhelix;
        double mass = 0;
        if (ID == 11) {
            zhelix = mdcKalTrack->getZHelixE();
            mass = 0.000511;
        }
        if (ID == 13) {
            zhelix = mdcKalTrack->getZHelixMu();
            mass = 0.105658;
        }
        if (ID == 211) {
            zhelix = mdcKalTrack->getZHelix();
            mass = 0.139570;
        }
        if (ID == 321) {
            zhelix = mdcKalTrack->getZHelixK();
            mass = 0.493677;
        }
        if (ID == 2212) {
            mass = 0.938272;
            zhelix = mdcKalTrack->getZHelixP();
        }

        double dr(0), phi0(0), kappa(0), dz(0), tanl(0);
        dr = zhelix[0];
        phi0 = zhelix[1];
        kappa = zhelix[2];
        dz = zhelix[3];
        tanl = zhelix[4];

        double pxy = 0;
        if (kappa != 0) pxy = 1.0 / fabs(kappa);
        double px = pxy * (-sin(phi0));
        double py = pxy * cos(phi0);
        double pz = pxy * tanl;
        double e = sqrt(pxy * pxy + pz * pz + mass * mass);
        return HepLorentzVector(px, py, pz, e);
    }
}

MCTruthMatchSvc::~MCTruthMatchSvc() {
    m_init = false;
    m_pars.erase(m_pars.begin(), m_pars.end());
}

StatusCode MCTruthMatchSvc::initialize() {
    MsgStream log(messageService(), name());
    log << MSG::INFO << "@initialize()" << endreq;

    StatusCode sc = Service::initialize();
    sc = serviceLocator()->service("EventDataSvc", eventSvc_, true);
    IMcDecayModeSvc* i_McDecayModeSvc;
    StatusCode sc_McDecayModeSvc = service("McDecayModeSvc", i_McDecayModeSvc);

    if (sc_McDecayModeSvc.isFailure()) {
        log << MSG::FATAL << "could not load McDecayModeSvc" << endreq;
        return sc_McDecayModeSvc;
    }
    m_svc = dynamic_cast<McDecayModeSvc*>(i_McDecayModeSvc);

    m_pars.erase(m_pars.begin(), m_pars.end());
    m_init = false;
    return sc;
}

StatusCode MCTruthMatchSvc::finalize() {
    MsgStream log(messageService(), name());
    log << MSG::INFO << "@initialize()" << endreq;

    StatusCode sc = Service::finalize();

    m_pars.erase(m_pars.begin(), m_pars.end());
    m_init = false;
    return sc;
}

StatusCode MCTruthMatchSvc::queryInterface(const InterfaceID& riid,
                                           void** ppvIF) {
    if (IMCTruthMatchSvc::interfaceID().versionMatch(riid)) {
        *ppvIF = dynamic_cast<IMCTruthMatchSvc*>(this);
    } else {
        return Service::queryInterface(riid, ppvIF);
    }
    addRef();
    // cout<<"MCTruthMatchSvc::Inf:queryInterface"<<endl;
    return StatusCode::SUCCESS;
}

bool MCTruthMatchSvc::match(vector<const EvtRecTrack*> recTracks,
                            vector<int> Pids, vector<int> pdgid,
                            vector<int> dechain) {
    if (dechain.empty()) {
        return false;
    }
    if (Pids.size() != recTracks.size()) {
        cout << "The PID list's size is not equl the number of tracks, please "
                "check it" << endl;
        return false;
    }

    if (pdgid.size() != dechain.size()) {
        cout << "Please check the decay chain!!!" << endl;
        return false;
    }

    dechain[0] = -1;

    SmartDataPtr<McParticleCol> mcParticleCol(eventSvc_,
                                              "/Event/MC/McParticleCol");

    map<int, McParticle*> parts;
    map<int, int> index;
    map<int, int> mother;
    int len = pdgid.size();

    // initial
    for (int i = 0; i < len; i++) {
        parts[i] = NULL;
        index[i] = 0;
    }

    for (int i = 0; i < len; i++) {
        for (int j = 0; j <= i; j++) {
            if (pdgid[j] == pdgid[i] && dechain[i] == dechain[j]) index[i] += 1;
        }
    }

    for (int i = 0; i < len; i++) {
        int mth = dechain[i];
        int pid = pdgid[i];

        if (mth < 0) {
            int n = index[i];
            McParticleCol::iterator ptr = mcParticleCol->begin();
            for (; ptr != mcParticleCol->end(); ptr++) {
                if ((*ptr)->particleProperty() == pid) {
                    if (n == 1) {
                        parts[i] = const_cast<McParticle*>(*ptr);
                        index[i] = index[i] - 1;
                        break;
                    } else
                        n = n - 1;
                }
            }
            if (parts[i] == NULL) {
                return false;
            }
        }

        if (mth >= 0) {
            int n = index[i];
            const SmartRefVector<McParticle>& decays =
                parts[mth]->daughterList();
            for (int j = 0; j < decays.size(); j++) {
                if (decays[j]->particleProperty() == pid) {
                    if (n == 1) {
                        parts[i] = const_cast<McParticle*>(decays[j].target());
                        index[i] = index[i] - 1;
                        break;
                    } else {
                        n = n - 1;
                    }
                }
            }
            if (parts[i] == NULL) {
                return false;
            }
        }
    }

    for (int k = 0; k < recTracks.size(); k++) {
        bool matched = false;
        const EvtRecTrack* atrack = recTracks[k];
        for (int p = 1; p < len; p++) {
            int pid = pdgid[p];
            if (!(abs(pid) == 11 || abs(pid) == 13 || pid == 22 ||
                  abs(pid) == 211 || abs(pid) == 321 || abs(pid) == 2212)) {
                continue;
            }

            if (matchTrack(atrack, Pids[k], parts[p])) {
                matched = true;
                break;
            }
        }
        if (!matched) return false;
    }
    return true;
}
void MCTruthMatchSvc::GetFinalChildren(McParticle* aEta,
                                       vector<McParticle*>& mcParticlaList) {
    mcParticlaList.erase(mcParticlaList.begin(), mcParticlaList.end());
    vector<int> _pid, _motherId;
    m_svc->extract(aEta, _pid, _motherId);
    if (_pid.empty()) {
        return;
    }
    map<int, McParticle*> parts;
    int* index = new int[_pid.size()];
    for (int i = 0; i < _pid.size(); i++) {
        parts[i] = NULL;
        index[i] = 0;
    }
    parts[0] = aEta;
    // initial
    for (int i = 0; i < _pid.size(); i++) {
        for (int j = 0; j <= i; j++) {
            if (_pid[j] == _pid[i] && _motherId[i] == _motherId[j]) {
                index[i] += 1;
            }
        }
    }
    parts[0] = aEta;
    for (int i = 1; i < _pid.size(); i++) {
        int n = index[i];
        const SmartRefVector<McParticle>& decays =
            parts[_motherId[i]]->daughterList();
        for (int j = 0; j < decays.size(); j++) {
            if (decays[j]->particleProperty() == _pid[i]) {
                if (index[i] == 1) {
                    parts[i] = const_cast<McParticle*>(decays[j].target());
                    break;
                }
                index[i] += -1;
            }
        }
        if (parts[i] == NULL) {
            return;
        }
    }
    for (int i = 1; i < _pid.size(); i++) {
        mcParticlaList.push_back(parts[i]);
    }
    return;
}
bool MCTruthMatchSvc::match(vector<const EvtRecTrack*> recTracks,
                            vector<int> Pids, const McParticle* aEta) {
    if (recTracks.size() != Pids.size() || recTracks.empty()) {
        return false;
    }

    vector<McParticle*> parts;
    this->GetFinalChildren(const_cast<McParticle*>(aEta), parts);

    if (parts.size() < recTracks.size()) {
        return false;
    }
    for (int k = 0; k < recTracks.size(); k++) {
        bool matched = false;
        // cout<<"pid:"<<Pids[k]<<endl;
        const EvtRecTrack* atrack = recTracks[k];
        for (int p = 0; p < parts.size(); p++) {

            if (matchTrack(atrack, Pids[k], parts[p])) {
                matched = true;
                break;
            }
        }
        if (!matched) {
            return false;
        }
    }

    return true;
}
bool MCTruthMatchSvc::match(vector<const EvtRecTrack*> recTracks,
                            vector<int> Pids, int motherPID) {
    SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc_,
                                                     "/Event/MC/McParticleCol");
    if (!mcParticleCol) {
        return false;
    }
    McParticleCol::iterator iter = mcParticleCol->begin();
    bool _ismatch = false;
    for (; iter != mcParticleCol->end(); iter++) {
        int m_pid = (*iter)->particleProperty();
        if (m_pid == motherPID && this->match(recTracks, Pids, (*iter))) {
            _ismatch = true;
            break;
        }
    }
    return _ismatch;
}
bool MCTruthMatchSvc::match(const EvtRecTrack* aTrack, int pid, int motherPID) {
    vector<const EvtRecTrack*> tracks;
    tracks.push_back(aTrack);
    vector<int> pids;
    pids.push_back(pid);
    return this->match(tracks, pids, motherPID);
}

RecMdcHitVector MCTruthMatchSvc::getMdcRecHits(McParticle* mcTrk) {
    RecMdcHitVector hitsCol;
    MdcMcHitCol::const_iterator itMcHits;
    RecMdcHitCol::const_iterator itRecHits;
    SmartDataPtr<MdcMcHitCol> mdcMcHits(eventSvc_, "/Event/MC/MdcMcHitCol");
    SmartDataPtr<RecMdcHitCol> recMdcHits(eventSvc_,
                                          "/Event/Recon/RecMdcHitCol");
    for (itMcHits = mdcMcHits->begin(); itMcHits != mdcMcHits->end();
         itMcHits++) {
        if ((*itMcHits)->getTrackIndex() != mcTrk->trackIndex()) continue;
        int mdcId = (*itMcHits)->identify().get_value();
        for (itRecHits = recMdcHits->begin(); itRecHits != recMdcHits->end();
             itRecHits++) {
            if ((*itRecHits)->getMdcId().get_value() == mdcId) {
                hitsCol.push_back(*itRecHits);
            }
        }
    }
    return hitsCol;
}
