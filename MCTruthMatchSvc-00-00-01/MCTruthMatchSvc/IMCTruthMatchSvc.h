#ifndef IMC_TRUTH_MATCH_SVC_H
#define IMC_TRUTH_MATCH_SVC_H

#include "GaudiKernel/IService.h"
#include "McTruth/McParticle.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "CLHEP/Vector/LorentzVector.h"

/* Decaration of the interface ID */
static const InterfaceID IID_IMCTruthMatchSvc("IMCTruthMatchSvc", 1, 0);

using CLHEP::HepLorentzVector;
class EvtRecMC;
namespace Event {
    class McParticle;
}

class IMCTruthMatchSvc : virtual public IService{
    public :

        virtual ~IMCTruthMatchSvc(){}
        static const InterfaceID& interfaceID(){ 
            return IID_IMCTruthMatchSvc;
        }
        virtual bool matchTrack(const EvtRecTrack* recTrk, Event::McParticle* mcTrk) = 0; 
        virtual bool matchTrack(const EvtRecTrack* recTrk,int pid, Event::McParticle* mcTrk) = 0; 
        virtual bool match(vector<const EvtRecTrack*> recTracks,vector<int> Pids,vector<int> pdgid,vector<int> dechain) = 0;
};

#endif
