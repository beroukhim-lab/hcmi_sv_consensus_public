# sv_consensus

# Original code breakdown:
Step 1: run_sv_merge.py
1. Get paths to each type of input (delly, snowman, dranger, brass)
2. Run below code to begin bash script that does actual processing

```
pcawg6_sv_merge_master.sh <run_id> <dranger_symlink> <snowman_symlink> <brass_symlink> <delly_symlink>
```

Step 2: pcawg6_sv_merge_master.sh (now called optimizing_sv_consensus.bash)


# Dependencies
pybedtools (pairToPair): https://daler.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.pair_to_pair.html
- Pretty sure original bedTools doesn't include -slop tag, so we need this version
