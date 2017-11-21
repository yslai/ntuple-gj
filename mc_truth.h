#include <AliMCEvent.h>

namespace {

    bool final_state_primary(AliMCEvent *mc_event, Int_t index)
    {
        if (mc_event->HasSubsidiaries()) {
            AliMCEvent *e = NULL;
            Int_t j = mc_event->FindIndexAndEvent(index, e);
            return final_state_primary(e, j);
        }
        else if (mc_event != NULL) {
            AliStack *s = mc_event->Stack();

            return index < s->GetNprimary() &&
                !(s->Particle(index)->GetStatusCode() > 1);
        }
        else {
            return false;
        }
    }

}
