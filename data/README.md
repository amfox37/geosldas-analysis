# External Test Data

The large reference datasets that previously lived under `GEOSldas_diagnostics/test_data` stay in that archived repo.  Keep this repo focused on code by:

1. Archiving `/Users/amfox/Desktop/GEOSldas_diagnostics` (or wherever you store the legacy tree).
2. Creating a local symlink that points back to the archived `test_data` directory:
   ```bash
   ln -s /Users/amfox/Desktop/GEOSldas_diagnostics/test_data data/external/test_data
   ```
3. Updating the path above if you relocate the archive.  Every notebook/script now references data via `data/external/test_data`, so the symlink is the only thing that needs to change per machine.

Never commit raw data into this repo—each `data` subfolder should contain only lightweight configs or README files that describe the required inputs.
