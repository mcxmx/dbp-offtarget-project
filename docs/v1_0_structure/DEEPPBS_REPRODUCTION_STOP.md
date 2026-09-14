# DeepPBS Reproduction Stop

The four published applicability cases (DBP5, DBP6, DBP9 and DBP35) now have
their design PDB inputs. Before any PBM comparison, the native DeepPBS
preprocessing was probed on a recovered design PDB. It stops at import/runtime:
this Windows workspace lacks the official Linux 3DNA/X3DNA preprocessing
runtime and the PyG compiled extensions (`torch_cluster`, `torch_scatter`, and
related packages). The available historical DBP35/DBP48 predictions are kept
with their original input paths and are not relabeled as v1.0 reruns.

This is an execution-environment limitation. It is not evidence that DeepPBS
cannot process the recovered structures and it does not support a biological
failure claim. A supported Linux/container run is required before making the
requested qualitative motif-reproduction or full DeepPBS coverage claim.
