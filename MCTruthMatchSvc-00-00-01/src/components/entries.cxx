#include "GaudiKernel/DeclareFactoryEntries.h"

#include "MCTruthMatchSvc/MCTruthMatchSvc.h"

DECLARE_SERVICE_FACTORY( MCTruthMatchSvc )

DECLARE_FACTORY_ENTRIES( MCTruthMatchSvc ) { 
  DECLARE_SERVICE( MCTruthMatchSvc );
}
