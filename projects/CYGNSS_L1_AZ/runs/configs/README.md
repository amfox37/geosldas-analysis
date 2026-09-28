# cygl1_operator_test run configs (fixed-operator experiments)

These are copies of the configs in `/discover/nobackup/projects/land_da/cygl1_operator_test/` (the project directory, which is not a
git repo), as of 2026-09-28. See `../../docs/cygl1_operator_test_project_README.md` for the experiment table and how to use them.

- `exeinp/<EXP_ID>.txt`: the `ldas_setup` exeinp file of each fixed-operator experiment. The four 2020–2022 arms (`*coh040216*`,
  `DA_L3_fixedop`, `DA_SMAP_fixedop`) were extended past their exeinp END_DATE by editing `run/CAP.rc` to `END_DATE: 20230101 000000`.
- `bat_inp/`: batch inputs (account s3208, 24 tasks, limited AZ domain).
- `templates/<name>/`: exeinp templates and `LDASsa_SPECIAL_inputs_ens{upd,prop}.nml`. Each experiment's `NML_INPUT_PATH` points at
  one of these directories.

The paths inside these files are absolute Discover paths (project dir, `scaling_params/`, obs dirs). If a config changes in the
project dir, re-copy it here.
