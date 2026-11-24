# MDR updates

The original MDR design document that Evan Brenner uncouthly deleted is restored 
from Emily Boyer's `85a1860`. An update to the design is proposed with the block
below and the code in this PR.

## MDR update proposal (2026-11-17)

- `vanilla` now selects drug and drug class models, while `mdr` selects MDR models.
- MDR runs only in `mdr` mode. Stratified, cross-test and leave-one-out modes all 
  exclude MDR. The computational cost and model assumptiondon't justify their inclusion.
- MAGIC rescoring and network export will remain limited to vanilla drug/class
  models for now. Native MDR coefficients retain their labels. Further MDR scoring
  and label handling require dedicated PRs!
- Optional parallel CV now changes execution only. Tuning candidates, including
  Ridge and below-path penalties retain previous eligibility behavior. Path notes 
  are diagnostic only. Fitting explicit penalty paths is deferred for later PR.
  
The comparisons in the original document are historical validation because some
guy accidentally deleted it. Its MDR cross-test and leave-one-out rules no longer 
describe supported modes, so the mdr.md document should be updated to reflect
these changes if they are approved.