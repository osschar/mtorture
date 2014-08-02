* intrinsics for non-symmetric matrices
* intrinsics for inversions
* intrinsics for custom fitting functions

* Track fitting test ... sigh. Somewhat unclear whether I should go back to
  mictest or get pieces of it here. Investigate ...
   * This is now partially imported in here.
   * On host: x23.7 (x22.8 with -O2),
   * On MIC:  crashes with -O2 and -O3 (x19 with -O1).

* One of kamlan mults might not have 3x3?
