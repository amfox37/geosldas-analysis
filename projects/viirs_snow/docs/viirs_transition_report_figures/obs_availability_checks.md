# Analysis A checks

- OL: n_obs == n_tiles on all days: True; days x species where an assimilated species has unassimilated obs: 0; where a monitor-only species has assimilated obs: 0
- T_to_J2: n_obs == n_tiles on all days: True; days x species where an assimilated species has unassimilated obs: 0; where a monitor-only species has assimilated obs: 0
- TA_to_J1J2: n_obs == n_tiles on all days: True; days x species where an assimilated species has unassimilated obs: 0; where a monitor-only species has assimilated obs: 0
- T_control: n_obs == n_tiles on all days: True; days x species where an assimilated species has unassimilated obs: 0; where a monitor-only species has assimilated obs: 0
- TA_control: n_obs == n_tiles on all days: True; days x species where an assimilated species has unassimilated obs: 0; where a monitor-only species has assimilated obs: 0

## Zero-obs periods (full days with no obs of the species)

| arm | species | first | last | days | assimilated there |
|---|---|---|---|---:|---|
| OL | Aqua | 20220402 | 20220416 | 15 | no |
| OL | Aqua | 20260102 | 20260731 | 211 | no |
| OL | Terra | 20221011 | 20221022 | 12 | no |
| OL | Terra | 20250511 | 20250514 | 4 | no |
| OL | Terra | 20260101 | 20260731 | 212 | no |
| OL | J2 | 20200802 | 20230209 | 922 | no |
| OL | J2 | 20240228 | 20240228 | 1 | no |
| T_to_J2 | Aqua | 20220402 | 20220416 | 15 | no |
| T_to_J2 | Aqua | 20260102 | 20260731 | 211 | no |
| T_to_J2 | Terra | 20221011 | 20221022 | 12 | yes |
| T_to_J2 | Terra | 20250511 | 20250514 | 4 | no |
| T_to_J2 | Terra | 20260101 | 20260731 | 212 | no |
| T_to_J2 | J2 | 20200802 | 20230209 | 922 | no |
| T_to_J2 | J2 | 20240228 | 20240228 | 1 | yes |
| TA_to_J1J2 | Aqua | 20220402 | 20220416 | 15 | yes |
| TA_to_J1J2 | Aqua | 20260102 | 20260731 | 211 | no |
| TA_to_J1J2 | Terra | 20221011 | 20221022 | 12 | yes |
| TA_to_J1J2 | Terra | 20250511 | 20250514 | 4 | no |
| TA_to_J1J2 | Terra | 20260101 | 20260731 | 212 | no |
| TA_to_J1J2 | J2 | 20200802 | 20230209 | 922 | no |
| TA_to_J1J2 | J2 | 20240228 | 20240228 | 1 | yes |
| T_control | Terra | 20250511 | 20250514 | 4 | yes |
| T_control | J2 | 20240228 | 20240228 | 1 | no |
| TA_control | Terra | 20250511 | 20250514 | 4 | yes |
| TA_control | J2 | 20240228 | 20240228 | 1 | no |
