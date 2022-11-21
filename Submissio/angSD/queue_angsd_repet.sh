
# qsub -cwd -V -N 'chr_X' -j y -q h11.q -l hostname=hercules11 -pe ompi64h11 24 -b y bash Submissio/queue_angsd.sh X
# qsub -cwd -V -N 'chr_U1' -j y -q h12.q -l hostname=hercules12 -pe ompi128h12 24 -b y bash Submissio/queue_angsd.sh U1
# qsub -cwd -V -N 'chr_1' -hold_jid 1469843 -j y -q h12.q -l hostname=hercules12 -pe ompi128h12 24 -b y bash Submissio/queue_angsd.sh 1

# qsub -cwd -V -N 'chr_U1_repet' -j y -q h12.q -l hostname=hercules12 -pe ompi128h12 24 -b y bash Submissio/queue_angsd.sh Data/'chr'
# qsub -cwd -V -N 'chr_U2_repet' -j y -q h11.q -l hostname=hercules11 -pe ompi64h11 24 -b y bash Submissio/queue_angsd.sh Data/'chr'

#angsd_theta_repet.sh Deng_reads.bam Dsil_genome.fasta 'chr' 'min' 'max'

# Per veure dades dels resultats: less output/chr_U2_mn*/Results_stats/*
# Per veure SFSs: less output/chr_U2_mn*/*.sfs

bash Submissio/angsd_theta_sites.sh Data/Deng_reads.bam Data/Dsil_genome.fasta $1 7 37
bash Submissio/angsd_theta_sites.sh Data/Deng_reads.bam Data/Dsil_genome.fasta $1 5 19
bash Submissio/angsd_theta_sites.sh Data/Deng_reads.bam Data/Dsil_genome.fasta $1 18 28
bash Submissio/angsd_theta_sites.sh Data/Deng_reads.bam Data/Dsil_genome.fasta $1 10 15
