#!/bin/bash
#SBATCH --job-name="post_proc_TBA"
#SBATCH --time=26:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute-p1,compute-p2
#SBATCH --mem-per-cpu=3800MB
#SBATCH --account=Research-AS-BN
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL

module load 2025 openmpi gromacs
module load 2025 py-pandas/2.2.3
module load 2025 python/3.11.9
module load 2025 py-matplotlib/3.9.2
module load 2025 py-scipy

NDX="index_extended.ndx"
TRAJ="traj_center.xtc"
TPR="md.tpr"

# Group map — verify with: echo q | gmx_mpi make_ndx -f md.tpr -n index_extended.ndx
#
#  0  System              37129 atoms
#  1  DNA                  1700 atoms
#  2  Na+                   124 atoms
#  3  K+                      2 atoms
#  4  Cl-                    75 atoms
#  5  Water               35228 atoms
#  6  SOL                 35228 atoms
#  7  non-Water            1901 atoms
#  8  T-tail                288 atoms
#  9  G-plex                490 atoms   (Gquad_O6 in earlier versions)
# 10  Gquad_O6               8 atoms
# 11  Gquad_K                2 atoms
# 12  Tetrad1_top           132 atoms   (even — gmx distance works)
# 13  Tetrad2_bot           133 atoms   (odd  — use mindist)
# 14  Loop1_TT               64 atoms
# 15  Loop2_TGT              97 atoms
# 16  Loop3_TT               64 atoms
# 17  TBA_loops             225 atoms
# 18  TBA_tetrads           265 atoms
# 19  r_1-23 / Immobilized   737 atoms
# 20  r_24-68 / TBA_strand   963 atoms
#
# If a Protein group is present it will appear as group 1, pushing
# all others down by the number of auto-generated protein groups.
# The detection block below handles this automatically.

# Detect whether a Protein group is present in the index file.
# All protein-specific steps are skipped if it is absent.
if grep -q '^\[ Protein \]' "${NDX}"; then
    HAS_PROTEIN=true
    echo ">>> Protein group detected — protein analyses will run."
    # With protein the group numbers from the thrombin index apply:
    #  1  Protein             4485 atoms
    #  3  C-alpha              277 atoms
    #  4  Backbone             831 atoms
    # 12  DNA                 2180 atoms
    # 19  T-tail               288 atoms
    # 20  G-plex               490 atoms
    # 21  Gquad_O6               8 atoms
    # 22  Gquad_K                4 atoms
    # 23  Tetrad1_bot          133 atoms  (odd — mindist)
    # 24  Tetrad2_top          132 atoms  (even — distance)
    # 25  Loop1_TT              64 atoms
    # 26  Loop2_TGT             97 atoms
    # 27  Loop3_TT              64 atoms
    # 28  TBA_loops            225 atoms
    # 29  Tba_tetrads          265 atoms
    # 30  Protein_DNA         6665 atoms
    DNA=12; GPLEX=29; TTAIL=19; GQUAD_O6=21; GQUAD_K=22
    TET1=23; TET2=24; L1=25; L2=26; L3=27; LOOPS=28; TETS=29
    PROT=1; PROT_CA=3; PROT_BB=4; PROT_DNA=30
else
    HAS_PROTEIN=false
    echo ">>> No Protein group found — protein analyses skipped."
    # DNA-only index numbers (as in this script's group map above)
    DNA=1; GPLEX=18; TTAIL=8; GQUAD_O6=10; GQUAD_K=11
    TET1=12; TET2=13; L1=14; L2=15; L3=16; LOOPS=17; TETS=18
fi

echo "================================================"
echo " Post-processing: TBA aptamer"
echo "================================================"


# ============================================================
# STEP A: TRAJECTORY CORRECTION
# ============================================================

echo ">>> STEP A: trajectory correction"

printf '0\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f md.xtc \
    -o traj_whole.xtc \
    -pbc whole

printf '0\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f traj_whole.xtc \
    -o traj_nojump.xtc \
    -pbc nojump
rm traj_whole.xtc

# Center on TBA strand, output all atoms
printf '30\n0\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f traj_nojump.xtc \
    -n ${NDX} \
    -o ${TRAJ} \
    -center -pbc mol -ur compact
rm traj_nojump.xtc

# Subsampled trajectories for quick visual inspection
printf '0\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f ${TRAJ} -o traj_center_10.xtc  -skip 10
printf '0\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f ${TRAJ} -o traj_center_100.xtc -skip 100

# PBC sanity check — minimum image distance must stay above the non-bonded cutoff (~1 nm)
printf "30\n" | srun gmx_mpi mindist \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -pi -od mindist_periodic.xvg -xvg none


# ============================================================
# STEP B: FITTED TRAJECTORY FOR PCA AND RMSF
#
# Fit on G-plex (stable core), output full DNA or Protein_DNA.
# This trajectory is used for ALL RMSF and PCA steps — overall
# tumbling must be removed here or per-residue fluctuations
# will be meaningless.
# ============================================================

echo ">>> STEP B: rotationally fitted trajectory"

if [ "$HAS_PROTEIN" = true ]; then
    printf "${GPLEX}\n${PROT_DNA}\n" | srun gmx_mpi trjconv \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -o traj_fit.xtc -fit rot+trans
else
    printf "${GPLEX}\n${DNA}\n" | srun gmx_mpi trjconv \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -o traj_fit.xtc -fit rot+trans
fi

echo ">>> STEP B2: per-group fitted trajectory extraction"
 
printf "${GPLEX}\n"  | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_apt.xtc
printf "${TTAIL}\n"  | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_Ttail.xtc
printf "${TETS}\n"   | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_Gquad.xtc
printf "${LOOPS}\n"  | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_loops.xtc
printf "${DNA}\n"    | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_DNA.xtc
 
if [ "$HAS_PROTEIN" = true ]; then
    printf "${PROT_DNA}\n" | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_complex.xtc
fi


# ============================================================
# STEP C: ENERGY EXTRACTION
# ============================================================

echo ">>> STEP C: energy"

printf 'Temperature\n0\n'   | srun gmx_mpi energy -f md.edr -o temperature.xvg  -xvg none
printf 'Pressure\n0\n'      | srun gmx_mpi energy -f md.edr -o pressure.xvg     -xvg none
printf 'Potential\n0\n'     | srun gmx_mpi energy -f md.edr -o potential.xvg    -xvg none
printf 'Kinetic-En.\n0\n'   | srun gmx_mpi energy -f md.edr -o kinetic-en.xvg   -xvg none
printf 'Total-Energy\n0\n'  | srun gmx_mpi energy -f md.edr -o total-energy.xvg -xvg none
printf 'Density\n0\n'       | srun gmx_mpi energy -f md.edr -o density.xvg      -xvg none
printf 'Volume\n0\n'        | srun gmx_mpi energy -f md.edr -o volume.xvg       -xvg none

# Protein-DNA interaction energies are only available if energygrps was set in the mdp.
# Uncomment and adjust term names if you reran with energygrps = Protein DNA:
# printf 'Coul-SR:Protein-DNA\n0\n' | srun gmx_mpi energy -f md.edr -o binding_coul.xvg -xvg none
# printf 'LJ-SR:Protein-DNA\n0\n'   | srun gmx_mpi energy -f md.edr -o binding_lj.xvg   -xvg none


# ============================================================
# STEP D: RMSD
# First group = fit reference; second group = what is measured.
# ============================================================

echo ">>> STEP D: RMSD"

# Whole DNA fitted on itself
printf "${DNA}\n${DNA}\n" | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_DNA.xvg -tu ns -xvg none

# G-plex core stability
printf "${GPLEX}\n${GPLEX}\n" | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_Gquad.xvg -tu ns -xvg none

# T-tail mobility relative to the G-plex core
printf "${GPLEX}\n${TTAIL}\n" | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_Ttail.xvg -tu ns -xvg none

# Individual loops fitted on G-plex
printf "${GPLEX}\n${L1}\n" | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_loop1.xvg -tu ns -xvg none
printf "${GPLEX}\n${L2}\n" | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_loop2.xvg -tu ns -xvg none
printf "${GPLEX}\n${L3}\n" | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_loop3.xvg -tu ns -xvg none

if [ "$HAS_PROTEIN" = true ]; then
    # Protein backbone stability
    printf "${PROT_BB}\n${PROT_BB}\n" | srun gmx_mpi rms \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -o rmsd_protein_bb.xvg -tu ns -xvg none
    # Protein fitted on itself, measuring DNA — captures relative rigid-body drift
    printf "${PROT_BB}\n${DNA}\n" | srun gmx_mpi rms \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -o rmsd_DNA_on_prot.xvg -tu ns -xvg none
fi


# ============================================================
# STEP E: RMSF (per-residue flexibility)
# Uses traj_fit.xtc — overall rotation must be removed first.
# One group piped; -res averages over atoms per residue.
# ============================================================

echo ">>> STEP E: RMSF"

printf "${GPLEX}\n" | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_Gquad.xvg -res -xvg none

printf "${TETS}\n" | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_tetrads.xvg -res -xvg none

printf "${LOOPS}\n" | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_loops.xvg -res -xvg none

printf "${L1}\n" | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_loop1.xvg -res -xvg none

printf "${L2}\n" | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_loop2.xvg -res -xvg none

printf "${L3}\n" | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_loop3.xvg -res -xvg none

printf "${TTAIL}\n" | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_Ttail.xvg -res -xvg none

# [added from doc2] Whole-strand flexibility
printf "${DNA}\n" | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_DNA.xvg -res -xvg none
    
if [ "$HAS_PROTEIN" = true ]; then
    # C-alpha RMSF — shows which protein residues are flexible
    printf "${PROT_CA}\n" | srun gmx_mpi rmsf \
        -s ${TPR} -f traj_fit.xtc -n ${NDX} \
        -o rmsf_protein_ca.xvg -res -xvg none
fi


# ============================================================
# STEP F: K+ ION ANALYSES
# Monitors whether K+ stays coordinated inside the G-quartet channel.
# Loss of K+ from the channel correlates with G-quad unfolding.
#
# Tetrad1_top / Tetrad2_top (132 atoms, even) — gmx distance works.
# Tetrad2_bot / Tetrad1_bot (133 atoms, odd)  — must use mindist.
# ============================================================

echo ">>> STEP F: K+ coordination"

# K+ to even-numbered tetrad (gmx distance computes per-atom distances to COM)
printf "${GQUAD_K}\n${TET1}\n" | srun gmx_mpi distance \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -oall kion_dist_tetrad1.xvg -xvg none

# K+ to odd-numbered tetrad (mindist as proxy)
printf "${GQUAD_K}\n${TET2}\n" | srun gmx_mpi mindist \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -od kion_dist_tetrad2.xvg -xvg none

# Number of O6 contacts within 3.5 Å (direct channel coordination criterion)
printf "${GQUAD_K}\n${TETS}\n" | srun gmx_mpi mindist \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -on kion_O6_contacts.xvg \
    -d 0.35 -xvg none

# Average K+ position relative to both tetrad planes
python3 << 'PYEOF'
import numpy as np, sys

def read_xvg(f):
    data = []
    with open(f) as fh:
        for line in fh:
            if line.startswith(('#', '@')): continue
            v = line.split()
            if v: data.append([float(x) for x in v])
    return np.array(data)

try:
    q1 = read_xvg('kion_dist_tetrad1.xvg')
    q2 = read_xvg('kion_dist_tetrad2.xvg')
    avg = (q1[:, 1] + q2[:, 1]) / 2.0
    np.savetxt('kion_dist_avg.xvg', np.column_stack([q1[:, 0], avg]), fmt='%.6f')
    print("Mean K+ channel distance : {:.2f} A".format(avg.mean() * 10))
    print("Fraction < 3.5 A (in-channel): {:.1f}%".format((avg * 10 < 3.5).mean() * 100))
except Exception as e:
    print("K+ averaging failed: {}".format(e), file=sys.stderr)
PYEOF


# ============================================================
# STEP G: TETRAD PLANE SEPARATION
# The odd-atom-count tetrad cannot use gmx distance for its COM.
# mindist gives the closest inter-tetrad atom distance, which is
# a reliable proxy for whether the planes stay stacked.
# ============================================================

echo ">>> STEP G: tetrad plane separation"

printf "${TET2}\n${TET1}\n" | srun gmx_mpi mindist \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -od gquad_tetrad_separation.xvg -xvg none


# ============================================================
# STEP H: RADIUS OF GYRATION
# ============================================================

echo ">>> STEP H: Rg"

printf "${DNA}\n"   | srun gmx_mpi gyrate -s ${TPR} -f ${TRAJ} -n ${NDX} -o gyrate_DNA.xvg   -xvg none
printf "${GPLEX}\n" | srun gmx_mpi gyrate -s ${TPR} -f ${TRAJ} -n ${NDX} -o gyrate_Gquad.xvg -xvg none
printf "${TTAIL}\n" | srun gmx_mpi gyrate -s ${TPR} -f ${TRAJ} -n ${NDX} -o gyrate_Ttail.xvg -xvg none

if [ "$HAS_PROTEIN" = true ]; then
    printf "${PROT}\n" | srun gmx_mpi gyrate -s ${TPR} -f ${TRAJ} -n ${NDX} -o gyrate_protein.xvg -xvg none
fi


# ============================================================
# STEP I: SASA
# ============================================================

echo ">>> STEP I: SASA"

printf "${DNA}\n"   | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_DNA.xvg   -xvg none
printf "${GPLEX}\n" | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_Gquad.xvg -xvg none
printf "${TTAIL}\n" | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_Ttail.xvg -xvg none
printf "${LOOPS}\n" | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_loops.xvg -xvg none
printf "${L1}\n"    | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_loop1.xvg -xvg none
printf "${L2}\n"    | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_loop2.xvg -xvg none
printf "${L3}\n"    | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_loop3.xvg -xvg none

if [ "$HAS_PROTEIN" = true ]; then
    printf "${PROT}\n" | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_protein.xvg -xvg none
fi


# ============================================================
# STEP J: HYDROGEN BONDS
#
# Intra-aptamer H-bonds tell you whether the G-quad scaffold and
# the binding-relevant surfaces stay intact over the simulation:
#   - G-quad internal: Hoogsteen base pairs that define the quartet
#   - T-tail to G-plex: T-tail is the binding-facing surface
#   - Loop1/3 to G-plex: lateral loops shape the binding geometry
#
# Protein-DNA H-bonds (when protein present) directly quantify
# the interface hydrogen bond network — the headline binding metric.
# NOTE: -ac was removed in GROMACS 2023+, not used here.
# ============================================================

echo ">>> STEP J: hydrogen bonds"

printf "${GPLEX}\n${GPLEX}\n" | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_Gquad_internal.xvg -xvg none

printf "${TTAIL}\n${GPLEX}\n" | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_Ttail_Gquad.xvg -xvg none

printf "${L1}\n${GPLEX}\n" | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_loop1_Gquad.xvg -xvg none

printf "${L3}\n${GPLEX}\n" | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_loop3_Gquad.xvg -xvg none

printf "${L2}\n${GPLEX}\n" | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_loop2_Gquad.xvg -xvg none
    
if [ "$HAS_PROTEIN" = true ]; then
    # Total protein-aptamer H-bonds over time
    printf "${PROT}\n${DNA}\n" | srun gmx_mpi hbond \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -num hbond_prot_DNA.xvg \
        -dist hbond_prot_DNA_dist.xvg \
        -ang hbond_prot_DNA_ang.xvg \
        -xvg none

    # Protein H-bonds to G-plex specifically (the core binding domain)
    printf "${PROT}\n${GPLEX}\n" | srun gmx_mpi hbond \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -num hbond_prot_Gquad.xvg \
        -dist hbond_prot_Gquad_dist.xvg \
        -ang hbond_prot_Gquad_ang.xvg \
        -xvg none

    # Protein H-bonds to T-tail
    printf "${PROT}\n${TTAIL}\n" | srun gmx_mpi hbond \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -num hbond_prot_Ttail.xvg -xvg none

    # Protein H-bonds to each loop
    printf "${PROT}\n${L1}\n" | srun gmx_mpi hbond \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -num hbond_prot_loop1.xvg -xvg none
    printf "${PROT}\n${L2}\n" | srun gmx_mpi hbond \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -num hbond_prot_loop2.xvg -xvg none
    printf "${PROT}\n${L3}\n" | srun gmx_mpi hbond \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -num hbond_prot_loop3.xvg -xvg none
fi


# ============================================================
# STEP K: INTER-MOLECULE MINIMUM DISTANCES
# With protein: measures whether each aptamer domain stays in
# contact with the protein over time.
# ============================================================

echo ">>> STEP K: inter-molecule contacts"

if [ "$HAS_PROTEIN" = true ]; then
    printf "${PROT}\n${DNA}\n"   | srun gmx_mpi mindist -s ${TPR} -f ${TRAJ} -n ${NDX} -od mindist_prot_DNA.xvg   -xvg none
    printf "${PROT}\n${GPLEX}\n" | srun gmx_mpi mindist -s ${TPR} -f ${TRAJ} -n ${NDX} -od mindist_prot_Gquad.xvg -xvg none
    printf "${PROT}\n${TTAIL}\n" | srun gmx_mpi mindist -s ${TPR} -f ${TRAJ} -n ${NDX} -od mindist_prot_Ttail.xvg -xvg none
    printf "${PROT}\n${L1}\n"    | srun gmx_mpi mindist -s ${TPR} -f ${TRAJ} -n ${NDX} -od mindist_prot_loop1.xvg -xvg none
    printf "${PROT}\n${L2}\n"    | srun gmx_mpi mindist -s ${TPR} -f ${TRAJ} -n ${NDX} -od mindist_prot_loop2.xvg -xvg none
    printf "${PROT}\n${L3}\n"    | srun gmx_mpi mindist -s ${TPR} -f ${TRAJ} -n ${NDX} -od mindist_prot_loop3.xvg -xvg none
fi


# ============================================================
# STEP L: INTERFACE CONTACT SELECTION
# Which aptamer atoms are within 4 Å of the protein at each frame.
# Only meaningful when the protein group is present.
# ============================================================

echo ">>> STEP L: interface contacts"

if [ "$HAS_PROTEIN" = true ]; then
    srun gmx_mpi select \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -select 'group "G-plex" and within 0.4 of group "Protein"' \
        -on interface_Gquad_prot.ndx \
        -oi interface_Gquad_contact_frac.dat -xvg none

    srun gmx_mpi select \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -select 'group "T-tail" and within 0.4 of group "Protein"' \
        -on interface_Ttail_prot.ndx \
        -oi interface_Ttail_contact_frac.dat -xvg none

    srun gmx_mpi select \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -select 'group "Loop1_TT" and within 0.4 of group "Protein"' \
        -on interface_loop1_prot.ndx \
        -oi interface_loop1_contact_frac.dat -xvg none

    srun gmx_mpi select \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -select 'group "Loop2_TGT" and within 0.4 of group "Protein"' \
        -on interface_loop2_prot.ndx \
        -oi interface_loop2_contact_frac.dat -xvg none

    srun gmx_mpi select \
        -s ${TPR} -f ${TRAJ} -n ${NDX} \
        -select 'group "Loop3_TT" and within 0.4 of group "Protein"' \
        -on interface_loop3_prot.ndx \
        -oi interface_loop3_contact_frac.dat -xvg none
else
    echo ">>> No protein group — interface contact selection skipped."
fi


# ============================================================
# STEP M: PRINCIPAL COMPONENT ANALYSIS
#
# All PCA uses traj_fit.xtc (fitted on G-plex).
# covar: TWO groups piped (fit group, then analysis group).
# anaeig: TWO groups piped (fit group, then analysis group)
#         with -s pointing to the average PDB from covar.
#
# Cosine content < 0.5: motion is real dynamics, not noise.
# Overlap between halves near 1: simulation has converged.
# ============================================================

echo ">>> STEP M: PCA"
# [added from doc2] Full DNA PCA
printf "${DNA}\n${DNA}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_DNA.xtc -n ${NDX} \
    -o eigenvalues_DNA.xvg \
    -v eigenvectors_DNA.trr \
    -av average_DNA.pdb \
    -xvg none
 
printf "${DNA}\n${DNA}\n" | srun gmx_mpi anaeig \
    -s average_DNA.pdb -f traj_fit_DNA.xtc -n ${NDX} \
    -v eigenvectors_DNA.trr \
    -eig eigenvalues_DNA.xvg \
    -proj pc_proj_DNA.xvg \
    -extr pc_extreme_DNA.pdb \
    -rmsf rmsf_pc_DNA.xvg \
    -first 1 -last 3 -xvg none
 
printf "${DNA}\n${DNA}\n" | srun gmx_mpi anaeig \
    -s average_DNA.pdb -f traj_fit_DNA.xtc -n ${NDX} \
    -v eigenvectors_DNA.trr \
    -eig eigenvalues_DNA.xvg \
    -2d pc2d_DNA.xvg \
    -first 1 -last 2 -xvg none
    
# G-plex (used as the "full aptamer" PCA here — stable core motions)
printf "${GPLEX}\n${GPLEX}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o eigenvalues_apt.xvg \
    -v eigenvectors_apt.trr \
    -av average_apt.pdb \
    -xvg none

printf "${GPLEX}\n${GPLEX}\n" | srun gmx_mpi anaeig \
    -s average_apt.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_apt.trr \
    -eig eigenvalues_apt.xvg \
    -proj pc_proj_apt.xvg \
    -extr pc_extreme_apt.pdb \
    -rmsf rmsf_pc_apt.xvg \
    -first 1 -last 3 -xvg none

printf "${GPLEX}\n${GPLEX}\n" | srun gmx_mpi anaeig \
    -s average_apt.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_apt.trr \
    -eig eigenvalues_apt.xvg \
    -2d pc2d_apt.xvg \
    -first 1 -last 2 -xvg none

# T-tail
printf "${TTAIL}\n${TTAIL}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o eigenvalues_Ttail.xvg \
    -v eigenvectors_Ttail.trr \
    -av average_Ttail.pdb \
    -xvg none

printf "${TTAIL}\n${TTAIL}\n" | srun gmx_mpi anaeig \
    -s average_Ttail.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_Ttail.trr \
    -eig eigenvalues_Ttail.xvg \
    -proj pc_proj_Ttail.xvg \
    -extr pc_extreme_Ttail.pdb \
    -rmsf rmsf_pc_Ttail.xvg \
    -first 1 -last 3 -xvg none

printf "${TTAIL}\n${TTAIL}\n" | srun gmx_mpi anaeig \
    -s average_Ttail.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_Ttail.trr \
    -eig eigenvalues_Ttail.xvg \
    -2d pc2d_Ttail.xvg \
    -first 1 -last 2 -xvg none

# G-quad core (TBA_tetrads)
printf "${TETS}\n${TETS}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o eigenvalues_Gquad.xvg \
    -v eigenvectors_Gquad.trr \
    -av average_Gquad.pdb \
    -xvg none

printf "${TETS}\n${TETS}\n" | srun gmx_mpi anaeig \
    -s average_Gquad.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_Gquad.trr \
    -eig eigenvalues_Gquad.xvg \
    -proj pc_proj_Gquad.xvg \
    -extr pc_extreme_Gquad.pdb \
    -rmsf rmsf_pc_Gquad.xvg \
    -first 1 -last 3 -xvg none

printf "${TETS}\n${TETS}\n" | srun gmx_mpi anaeig \
    -s average_Gquad.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_Gquad.trr \
    -eig eigenvalues_Gquad.xvg \
    -2d pc2d_Gquad.xvg \
    -first 1 -last 2 -xvg none

# Loops
printf "${LOOPS}\n${LOOPS}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o eigenvalues_loops.xvg \
    -v eigenvectors_loops.trr \
    -av average_loops.pdb \
    -xvg none

printf "${LOOPS}\n${LOOPS}\n" | srun gmx_mpi anaeig \
    -s average_loops.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_loops.trr \
    -eig eigenvalues_loops.xvg \
    -proj pc_proj_loops.xvg \
    -extr pc_extreme_loops.pdb \
    -rmsf rmsf_pc_loops.xvg \
    -first 1 -last 3 -xvg none

printf "${LOOPS}\n${LOOPS}\n" | srun gmx_mpi anaeig \
    -s average_loops.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_loops.trr \
    -eig eigenvalues_loops.xvg \
    -2d pc2d_loops.xvg \
    -first 1 -last 2 -xvg none

# Cosine content — should be < 0.5 for meaningful dynamics
srun gmx_mpi analyze -f pc_proj_apt.xvg   -cc cosine_content_apt.xvg   -xvg none
srun gmx_mpi analyze -f pc_proj_Ttail.xvg -cc cosine_content_Ttail.xvg -xvg none
srun gmx_mpi analyze -f pc_proj_Gquad.xvg -cc cosine_content_Gquad.xvg -xvg none
srun gmx_mpi analyze -f pc_proj_loops.xvg -cc cosine_content_loops.xvg -xvg none
srun gmx_mpi analyze -f pc_proj_DNA.xvg   -cc cosine_content_DNA.xvg   -xvg none

# Convergence: split trajectory in half, compare PCA subspace overlap
HALF_TIME=$(python3 << 'PYEOF'
import subprocess
result = subprocess.run(
    ['gmx_mpi', 'check', '-f', 'traj_fit.xtc'],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
output = result.stdout.decode() + result.stderr.decode()
for line in output.split('\n'):
    if 'Last frame' in line:
        parts = line.split()
        for i, p in enumerate(parts):
            if p == 'time' and i + 1 < len(parts):
                print(float(parts[i + 1]) / 2.0)
                break
        break
PYEOF
)
echo "Splitting trajectory at ${HALF_TIME} ps"

if [ "$HAS_PROTEIN" = true ]; then
    printf "${PROT_DNA}\n" | srun gmx_mpi trjconv \
        -s ${TPR} -f traj_fit.xtc -n ${NDX} \
        -o traj_fit_h1.xtc -e ${HALF_TIME}
    printf "${PROT_DNA}\n" | srun gmx_mpi trjconv \
        -s ${TPR} -f traj_fit.xtc -n ${NDX} \
        -o traj_fit_h2.xtc -b ${HALF_TIME}
else
    printf "${DNA}\n" | srun gmx_mpi trjconv \
        -s ${TPR} -f traj_fit.xtc -n ${NDX} \
        -o traj_fit_h1.xtc -e ${HALF_TIME}
    printf "${DNA}\n" | srun gmx_mpi trjconv \
        -s ${TPR} -f traj_fit.xtc -n ${NDX} \
        -o traj_fit_h2.xtc -b ${HALF_TIME}
fi

# T-tail convergence
printf "${TTAIL}\n${TTAIL}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_Ttail_h1.xvg -v eigenvectors_Ttail_h1.trr -xvg none
printf "${TTAIL}\n${TTAIL}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_Ttail_h2.xvg -v eigenvectors_Ttail_h2.trr -xvg none
printf "${TTAIL}\n${TTAIL}\n" | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_Ttail_h1.trr \
    -eig eigenvalues_Ttail_h1.xvg \
    -over overlap_Ttail.xvg \
    -first 1 -last 10 -xvg none

# G-quad convergence
printf "${TETS}\n${TETS}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_Gquad_h1.xvg -v eigenvectors_Gquad_h1.trr -xvg none
printf "${TETS}\n${TETS}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_Gquad_h2.xvg -v eigenvectors_Gquad_h2.trr -xvg none
printf "${TETS}\n${TETS}\n" | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_Gquad_h1.trr \
    -eig eigenvalues_Gquad_h1.xvg \
    -over overlap_Gquad.xvg \
    -first 1 -last 10 -xvg none

# Loops convergence
printf "${LOOPS}\n${LOOPS}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_loops_h1.xvg -v eigenvectors_loops_h1.trr -xvg none
printf "${LOOPS}\n${LOOPS}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_loops_h2.xvg -v eigenvectors_loops_h2.trr -xvg none
printf "${LOOPS}\n${LOOPS}\n" | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_loops_h1.trr \
    -eig eigenvalues_loops_h1.xvg \
    -over overlap_loops.xvg \
    -first 1 -last 10 -xvg none

# G-plex (apt) convergence
printf "${GPLEX}\n${GPLEX}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_apt_h1.xvg -v eigenvectors_apt_h1.trr -xvg none
printf "${GPLEX}\n${GPLEX}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_apt_h2.xvg -v eigenvectors_apt_h2.trr -xvg none
printf "${GPLEX}\n${GPLEX}\n" | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_apt_h1.trr \
    -eig eigenvalues_apt_h1.xvg \
    -over overlap_apt.xvg \
    -first 1 -last 10 -xvg none

# [added from doc2] Full DNA convergence
printf "${DNA}\n${DNA}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_DNA_h1.xvg -v eigenvectors_DNA_h1.trr -xvg none
printf "${DNA}\n${DNA}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_DNA_h2.xvg -v eigenvectors_DNA_h2.trr -xvg none
printf "${DNA}\n${DNA}\n" | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_DNA_h1.trr \
    -eig eigenvalues_DNA_h1.xvg \
    -over overlap_DNA.xvg \
    -first 1 -last 10 -xvg none


# ============================================================
# STEP N: DYNAMIC CROSS-CORRELATION (DCCM)
# Writes the covariance matrix as ASCII for Python plotting.
# Uses full DNA (or Protein_DNA when protein is present).
# ============================================================

echo ">>> STEP N: DCCM"


if [ "$HAS_PROTEIN" = true ]; then
    printf "${PROT_DNA}\n${PROT_DNA}\n" | srun gmx_mpi covar \
        -s ${TPR} -f traj_fit_complex.xtc -n ${NDX} \
        -ascii dccm_complex.dat -xvg none
fi
 
# [added from doc2] Full DNA group (DNA-only, or as standalone aptamer context)
printf "${DNA}\n${DNA}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_DNA.xtc -n ${NDX} \
    -ascii dccm_DNA.dat -xvg none
 
# TBA aptamer (G-plex) only
printf "${GPLEX}\n${GPLEX}\n" | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_apt.xtc -n ${NDX} \
    -ascii dccm_apt.dat -xvg none
 

# ============================================================
# DONE
# ============================================================

echo ""
echo "================================================"
echo " Analysis complete."
echo "================================================"
echo ""
echo "Sanity checks first:"
echo "  mindist_periodic.xvg        should stay above non-bonded cutoff (~1 nm)"
echo "  cosine_content_*.xvg        should be < 0.5 for meaningful PCA"
echo "  overlap_Ttail/Gquad/loops/apt.xvg   convergence: closer to 1 = better"
echo ""
echo "DNA structural integrity:"
echo "  rmsd_DNA/Gquad/Ttail/loop*.xvg"
echo "  rmsf_Gquad/tetrads/loops/loop*/Ttail.xvg"
echo "  kion_dist_avg.xvg  kion_O6_contacts.xvg  gquad_tetrad_separation.xvg"
echo "  gyrate_*.xvg  sasa_*.xvg"
echo "  hbond_Gquad_internal.xvg  hbond_Ttail_Gquad.xvg  hbond_loop1/3_Gquad.xvg"
echo ""
if [ "$HAS_PROTEIN" = true ]; then
echo "Binding interface (protein present):"
echo "  rmsd_protein_bb.xvg  rmsd_DNA_on_prot.xvg  rmsf_protein_ca.xvg"
echo "  hbond_prot_DNA/Gquad/Ttail/loop*.xvg"
echo "  mindist_prot_DNA/Gquad/Ttail/loop*.xvg"
echo "  interface_*_contact_frac.dat"
echo "  dccm_complex.dat"
else
echo "  dccm_apt.dat"
fi
echo ""
echo "PCA:"
echo "  eigenvalues_apt/Ttail/Gquad/loops.xvg"
echo "  pc2d_apt/Ttail/Gquad/loops.xvg"
echo "  pc_extreme_apt/Ttail/Gquad/loops.pdb"
echo ""
echo "Now run: python3 post_analysis_TBA.py"
