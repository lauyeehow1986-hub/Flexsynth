# flexsynth 0.0.0.9000

* Phase 0 scaffold: package skeleton, MIT license, and CI (`R-CMD-check`).
* Public API with settled signatures and input validation: `synth()`,
  `synth_linked()`, `synth_control()`, `dp_control()`.
* Two-track privacy design documented: high-utility default (Track A) and an
  opt-in person-level differentially private track (Track B).
* Synthetic cardiac example-data generator (`data-raw/make_toy_cardiac.R`).
