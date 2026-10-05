# Invalidated attempt

Slurm job `1708915` was stopped after the first configuration completed.
The original dispatcher submitted one future and immediately waited for it,
so nominal worker counts greater than one were effectively serial. The
partial `default_w1` record is retained as diagnostic evidence only; this run
root must not be used for worker selection or production sizing.

The wrapper was corrected to use the shared bounded dispatcher. The validated
replacement uses a new run root.
