---
title: Structural Variant Calling
nav_order: 7
---

# Structural Variant Calling

In conjunction with a recent structural variant caller, BISCUIT can be used to make large-scale structural variant
calls. We recommend using either `lumpyexpress` or `manta`, but the standard-compliant BAM output by BISCUIT should work
with other similar callers.

## Recommended Structural Variant Calling Pipeline

  1. Create aligned, sorted, and duplicate marked BAM using the biscuitSifter pipeline. Note, while we would not suggest
  using the `-M` flag otherwise, for structural variant calling, we suggest using the `-M` flag in `biscuit` when
  running the pipeline. Also, if planning to use `lumpy` rather than `lumpyexpress`, a different pipeline must be used, as
  `lumpy` requires split and discordant reads to be separated from the original BAM. However, most other structural
  variant callers extract the split and discordant reads during processing, so this is generally not needed.
  2. Using the [BISCUIT SV calling container available on GitHub](https://github.com/huishenlab/sv_calling_docker), call
  structural variants using your preferred structural variant caller. Out of the box, the container includes Manta,
  Delly, Smoove, and Lumpy, but it can be easily modified to include other callers.

For callers that allow the user to restrict to specified call regions, we suggest (at minimum) restricting to primary
chromosomes (not including the mitochondrial chromosome) and only taking regions not in the ENCODE black list of your
genome of choice. The callable regions can be further restricted depending on your specific needs.

As with most things in biological research, structural variant calling is hard, so this recommendation may not be the
best option for every case. *It is purely a recommendation and a reasonable place to start.*
