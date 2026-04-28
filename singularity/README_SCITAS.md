# Copying data to cluster

Some rsync commands to transfer data going through the Haas because it seems faster.

#### To a cluster for processing

Selective rsync copy of files to cluster: remove `--dry-run` for the real copy
```bash
rsync -aPv --dry-run   --include='*corrected*'   --include='*/'   --exclude='*'   bisi@haas056.rcp.epfl.ch:/mnt/lsens-analysis/Axel_Bisi/data/AB163/   /scratch/bisi/data/AB163/
```

Same as before but for a list of mouse folders, parallelized and with SSH multiplexing (1 pwd prompt):

```bash
printf "%s\n" AB{077..164} MH{001..075} | \
xargs -P 4 -I {} rsync -aPv --dry-run \
  -e "ssh -o ControlMaster=auto -o ControlPersist=10m -o ControlPath=~/.ssh/cm-%r@%h:%p" \
  --include='*corrected*' \
  --include='*/' \
  --exclude='*' \
  bisi@haas056.rcp.epfl.ch:/mnt/lsens-analysis/Axel_Bisi/data/{}/ \
  /scratch/bisi/data/{}/
  ```

Similar, but change `-aPn` into `-aPv` to remove the `--dry-run`:

```bash
printf "%s\n" AB{160..164} MH{001..010} | xargs -P 4 -I {} rsync -aPn   -e "ssh -o ControlMaster=auto -o ControlPersist=10m -o ControlPath=~/.ssh/cm-%r@%h:%p"  --include='*corrected*'   --include='*/'   --exclude='*'   bisi@haas056.rcp.epfl.ch:/mnt/lsens-analysis/Axel_Bisi/data/{}/   /scratch/bisi/data/{}/(modifié)
```

#### From a cluster after processing

Selective rsync copy of folders and associated files from cluster: same remove `--dry-run`

```bash
rsync -aPv \
  --include='*/' \
  --include='dredge/***' \
  --exclude='*' \
  /scratch/bisi/data/AB163/ \
  bisi@haas056.rcp.epfl.ch:/mnt/lsens-analysis/Axel_Bisi/data/AB163/
```

