# `robin list-job-types` and job model

`robin list-job-types` prints the **job types** ROBIN knows about, grouped by **queue**. Understanding this helps you build `-w` / `--workflow` strings for [`robin workflow`](workflow.md).

## Command

```bash
robin list-job-types
```

You must complete the **disclaimer** (`I agree`).

## Queues and job types

The orchestration layer assigns each job type to a queue (simplified names here; internal Ray queue names may differ slightly):

| Queue (concept) | Job types | Role |
|-----------------|-----------|------|
| **Preprocessing** | `preprocessing` | Read BAM headers/metadata; entry point for each new file. |
| **BED conversion** | `bed_conversion` | Prepare inputs for classifiers that need BED-level views. |
| **Analysis** | `mgmt`, `cnv`, `target`, `fusion` | Methylation (MGMT), copy number, targeted variant/fusion panels. |
| **Classification** | `sturgeon`, `nanodx`, `pannanodx` | Methylation / expression classifiers. |
| **Slow** | `random_forest` | Heavier or batched models. |

## Automatic steps

When you use the **simplified** workflow format (`-w mgmt,sturgeon`, …):

1. **`preprocessing`** is prepended if you did not list it.
2. **`bed_conversion`** is inserted when any of **`sturgeon`**, **`nanodx`**, **`pannanodx`**, or **`random_forest`** appear — those jobs expect BED conversion upstream.

You do not need to list `bed_conversion` manually for those classifiers unless you are hand-editing **legacy** queue-prefixed pipelines.

## Valid job type names

The CLI accepts only these **job** identifiers in workflow strings:

`preprocessing`, `bed_conversion`, `mgmt`, `cnv`, `target`, `fusion`, `sturgeon`, `nanodx`, `pannanodx`, `random_forest`

Unknown names produce warnings and are skipped.

## Examples

```bash
# Minimal classifier run (preprocessing + bed_conversion added as needed)
robin workflow /data/bams -w sturgeon --center Demo --target-panel rCNS2 -d /out --reference /ref/hg38.fa

# Full stack (typical)
robin workflow /data/bams \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --center Sherwood \
  --target-panel rCNS2 \
  -d ~/results \
  --reference ~/references/hg38.fa
```

## Related

- [`robin workflow`](workflow.md)  
- [Quickstart](../getting-started/quickstart.md)  
