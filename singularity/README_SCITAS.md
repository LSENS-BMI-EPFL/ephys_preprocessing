# Copying data to cluster

### rsync
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
  --no-perms \
  --no-times \
  --include='*/' \
  --include='dredge_fast/***' \
  --include='kilosort4/***' \
  --include='preprocess/***' \
  --exclude='*' \
  /scratch/bisi/data/AB163/ \
  bisi@haas056.rcp.epfl.ch:/mnt/lsens-analysis/Axel_Bisi/data/AB163/
```

Or for a list of mouse folders, parallelized and with SSH multiplexing (1 pwd prompt):

````bash
printf "%s\n" AB{163..164} | \
xargs -P 4 -I {} rsync -av --no-perms --no-times --dry-run \
  -e "ssh -o ControlMaster=auto -o ControlPersist=10m -o ControlPath=~/.ssh/cm-%r@%h:%p" \
  --include='*/' \
  --include='dredge/***' \
  --exclude='*' \
  /scratch/bisi/data/{}/ \
  bisi@haas056.rcp.epfl.ch:/mnt/lsens-analysis/Axel_Bisi/data/{}/
````

### rclone
Or using rclone directly in a SCITAS cluster:

```
module load rclone
```

Configure: ``rclone config``

```
[NAS-SV]
type = smb
host = sv-nas1.rcp.epfl.ch
user = USER
pass = PASSWORD
domain = INTRANET.EPFL.CH
spn = sv-nas1.rcp.epfl.ch
```

check the remote works:
```
rclone lsd nas-sv:"Petersen-Lab/analysis/Axel_Bisi/data/"
``` 
It should list the content of the folders. Then run:

```bash
rclone copy nas-sv:"Petersen-Lab/analysis/Axel_Bisi/data/" "/scratch/bisi/data/"   --filter "+ */"   --filter "+ **/*corrected*"   --filter "- *"  --transfers 8 --checkers 8 --progress --log-file rclone_transfer.log
```

Or, with more files requires in preprocessing:
````bash
rclone copy \
  nas-sv:"Petersen-Lab/analysis/Axel_Bisi/data/" \
  "/scratch/bisi/data/" \
  --filter "+ */" \
  --filter "+ **/catgt_*/**/*corrected*" \
  --filter "+ **/catgt_*/**/*.meta" \
  --filter "+ **/catgt_*/*.meta" \
  --filter "+ **/catgt_*/*.txt" \
  --filter "+ **/catgt_*/**/*.txt" \
  --filter "+ **/tracks/**" \
  --filter "- *" \
  --transfers 8 \
  --checkers 8 \
  --progress \
  --stats 30s \
  --log-level INFO \
  --log-file rclone_transfer.log
````
Make sure the filters are correct.