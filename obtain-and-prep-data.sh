

# 1. Download the indexed Chinook salmon genome, and put it in a location in scratch
# somewhere. (I am being intentionally vague so that people can put this where they
# want to)

# Most everyone should already have downloaded this indexed genome into the directory
# chinook-genome-idx in the chinook-play homework:

# Here is where I have the genome file stored (and the index is with it):
# /home/eanderson/scratch/chinook-play/chinook-genome-idx/GCA_002872995.1_Otsh_v1.0_genomic.fna.gz


# Note that if you need to download it you can find the directory called 
# pre-indexed-chinook-genome in our shared course
# folder at:  https://drive.google.com/drive/folders/11fLG7b0RV1Uij9CbYh_jt_Xp_dyfLetO?usp=sharing
# that I have shared with everyone in the class already

# If I had needed to download this, I would have used a command like this:
rclone copy --tpslimit 10 --fast-list -P --drive-shared-with-me gdrive-pers:CSU-con-gen-2020/pre-indexed-chinook-genome  pre-indexed-chinook-genome

# You would need to modify that for your own system.


# 2. Download the 2,560 gzipped fastq files.  These are in the shared folder at
#  https://drive.google.com/drive/folders/11fLG7b0RV1Uij9CbYh_jt_Xp_dyfLetO?usp=sharing
# inside the directory fastqs-chr32-160-chinook-8-lanes.  Use rclone to get this
# to scratch on your cluster.  Since there are so many files before sure to
# use the --tps-limit 19 and --fast-list options.

# Here is how I did it:
rclone copy --tpslimit 10 --fast-list -P --drive-shared-with-me gdrive-pers:CSU-con-gen-2020/fastqs-chr32-160-chinook-8-lanes  fastqs-chr32-160-chinook-8-lanes




# 3. Move the folder with the preindexed genome into this repository directory
# and rename it to "genome" (or make a Symbolic link).

# Here is what I did (making a symbolic link)
ln -s /home/eanderson/scratch/chinook-play/chinook-genome-idx genome

# You will have to modify this to reflect the proper location of the
# directory on your system.


# 4. Move (or make a symbolic link), the directory fastqs-chr32-160-chinook-8-lanes containing
# all the gzipped FASTQ files into a directory named "fastqs" in this repository directory.

# Here is what I did (making a symbolic link)
ln -s  /home/eanderson/scratch/course_stuff/fastqs-chr32-160-chinook-8-lanes fastqs

# You will have to modify this to reflect the proper location of the
# directory on your system.




# With all of this in place.  You are ready to get an interactive
# shell on a compute node and make sure that you can make the command
# lines in map-N-files-from-K.sh work.  