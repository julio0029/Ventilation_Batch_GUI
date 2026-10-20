# Evidence map: species specific NeuroKit2 ventilation methods

| Claim | Status | Original source and location | Version or date | Verified | Uncertainty or limitation |
|---|---|---|---|---|---|
| Analysis uses species specific regions of interest and image difference motion | CONFIRMED | `ventilation_batch_gui.py`, `extract_motion`, `VentilationGUI.choose_roi` | 21 October 2026 | Source reviewed and syntax checked | Biological validity requires manually labelled data. |
| NeuroKit2 general peak detection is used for the image motion trace | CONFIRMED | `ventilation_batch_gui.py`, `neurokit_peaks`; Makowski et al., 2021, doi:10.3758/s13428-020-01516-y | 21 October 2026 | Clean Python 3.11 synthetic signal test | NeuroKit2 respiratory belt functions are not used because the input is not a respiratory belt signal. |
| Epaulette shark and shovelnose ray labels remain separate in results | CONFIRMED | `ventilation_batch_gui.py`, `safe_metadata` | 21 October 2026 | Synthetic label tests | Scientific names were not inferred from the user's local labels. |
| Epaulette shark respiratory response to progressive hypoxia is a relevant species specific physiological context | CONFIRMED | Routley et al., 2002, doi:10.1016/S1095-6433(01)00484-6 | 2002 | Bibliographic metadata verified | Does not validate the video algorithm. |
| Shovelnose ray hypoxia context is relevant | PRELIMINARY | Speers Roesch et al., 2012, doi:10.1242/jeb.059642 | 2012 | Bibliographic metadata verified | Confirm the scientific identity of the local shovelnose ray before citing species specific physiology. |
| Automated estimates are valid for hypothesis testing | UNKNOWN | No annotated study videos or manual counts available | 21 October 2026 | No | Compare against blinded manual counts stratified by species, oxygen treatment, drug treatment, orientation, and camera geometry. |
