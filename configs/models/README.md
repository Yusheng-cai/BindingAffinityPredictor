# Model Configurations

Each model configuration should eventually record:

- adapter name;
- upstream source URL and pinned revision;
- checkpoint identifier, revision, and checksum;
- environment/lock-file reference;
- native input and output conventions;
- default inference arguments;
- score units and direction;
- sampling, seed, recycling, MSA, template, and cache policy;
- known hardware limitations.

Configurations with `status: proposed_not_installed` may be used to preregister
an installation and smoke test. They are not executable sources of truth until
the installed revision, checkpoint identity, and native output have been
verified; only then should their status change to `verified`.
