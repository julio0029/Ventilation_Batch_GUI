# Evidence map: ventilation analysis methods draft

| Claim | Primary source | Verification state | Notes |
|---|---|---|---|
| Recursive local video processing and supported formats | `ventilation_batch_gui.py`, lines 20 and batch scanning implementation | CONFIRMED | Current GUI source. |
| Metadata fields and keyword inference | `ventilation_batch_gui.py`, `safe_metadata` | CONFIRMED | Automatic labels require experimental record check. |
| Region of interest handling and frame difference motion signal | `ventilation_batch_gui.py`, `extract_motion` | CONFIRMED | Uses greyscale mean absolute consecutive frame difference. |
| Smoothing, detrending, spectral and peak calculations | `ventilation_batch_gui.py`, `moving_average`, `linear_detrend`, `estimate_rate`, `local_peaks` | CONFIRMED | Defaults are 0.07 s smoothing, 20 to 180 bpm, prominence 0.5. |
| Diagnostic plots and CSV result fields | `ventilation_batch_gui.py`, `save_diagnostic`, `analyse_video` | CONFIRMED | Output reflects the current implementation. |
| Biological validity against manual counts | No completed validation dataset in the project | UNKNOWN | Must be completed and reported before manuscript submission. |
| Camera setup and exposure protocol | No verified experimental record read for this draft | UNKNOWN | Add from laboratory records. |
