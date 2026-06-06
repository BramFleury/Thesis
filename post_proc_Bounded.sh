#!/bin/bash
#SBATCH --job-name="post_proc_Bounded"
#SBATCH --time=8:00:00
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

# Group map (verify with: echo q | gmx_mpi make_ndx -f md.tpr -n index_extended.ndx)
#  0  System              37129 atoms
#  1  DNA                  1700 atoms
#  2  Na+                   124 atoms
#  3  K+                      2 atoms
#  4  Cl-                    75 atoms
#  5  Water               35228 atoms
#  6  SOL                 35228 atoms
#  7  non-Water            1901 atoms
#  8  T-tail                288 atoms
#  9  Gquad_O6               8 atoms
# 10  Gquad_K                2 atoms
# 11  Tetrad2_bot           133 atoms
# 12  Tetrad1_top           132 atoms
# 13  Loop1_TT               64 atoms
# 14  Loop2_TGT              97 atoms
# 15  Loop3_TT               64 atoms
# 16  TBA_loops             225 atoms
# 17  TBA_tetrads           265 atoms
# 18  G-plex                490 atoms
# 19  Immobilized_strand    737 atoms
# 20  TBA_strand            963 atoms
#
# DNA-only simulation: no Protein group exists.
# Immobilized_strand (19) is used as the binding partner throughout.

echo "================================================"
echo " Post-processing: TBA aptamer"
echo "================================================"


# ============================================================
# STEP A: TRAJECTORY CORRECTION (periodicity)
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

# Center on whole aptamer, output all atoms
printf '1\n0\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f traj_nojump.xtc \
    -o ${TRAJ} \
    -center -pbc mol -ur compact
rm traj_nojump.xtc

# Subsampled trajectories (optional, for quick visual inspection)
printf '0\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f ${TRAJ} -o traj_center_10.xtc  -skip 10
printf '0\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f ${TRAJ} -o traj_center_100.xtc -skip 100

# PBC sanity check: minimum image distance should stay above the non-bonded cutoff
printf '1\n' | srun gmx_mpi mindist \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -pi -od mindist_periodic.xvg -xvg none


# ============================================================
# STEP B: FITTED TRAJECTORY FOR PCA AND RMSF
# Fit rotation+translation to G-plex (18), write all DNA (1).
# Both PCA and RMSF must use this trajectory so that internal
# fluctuations are not contaminated by overall tumbling.
# ============================================================

echo ">>> STEP B: rotationally fitted trajectory"

printf '18\n1\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o traj_fit.xtc \
    -fit rot+trans

printf '18\n'  | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_apt.xtc
printf '8\n'  | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_Ttail.xtc
printf '17\n' | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_Gquad.xtc
printf '16\n' | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_loops.xtc
printf '1\n'  | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_DNA.xtc
printf '19\n' | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_imm.xtc
printf '20\n' | srun gmx_mpi trjconv -s ${TPR} -f traj_fit.xtc -n ${NDX} -o traj_fit_TBAstrand.xtc

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


# ============================================================
# STEP D: RMSD
# Format: printf 'fit_group\nanalysis_group\n'
# Fit group = reference for least-squares superposition.
# Analysis group = atoms whose RMSD is measured.
# ============================================================

echo ">>> STEP D: RMSD"

# Whole DNA, fitted on itself
printf '1\n1\n' | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_DNA.xvg -tu ns -xvg none

# G-plex core stability
printf '17\n17\n' | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_Gquad.xvg -tu ns -xvg none

# T-tail mobility relative to the G-plex core
printf '17\n8\n' | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_Ttail.xvg -tu ns -xvg none
    

# Individual loops fitted on G-plex (loop motion relative to core)
printf '17\n13\n' | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_loop1.xvg -tu ns -xvg none
printf '17\n14\n' | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_loop2.xvg -tu ns -xvg none
printf '17\n15\n' | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_loop3.xvg -tu ns -xvg none

printf '19\n20\n' | srun gmx_mpi rms \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -o rmsd_TBAstrand_vs_imm.xvg -tu ns -xvg none
    
# ============================================================
# STEP E: RMSF (per-residue flexibility)
# Must use traj_fit.xtc so overall rotation is removed.
# Only one group is piped: the group whose fluctuations you measure.
# The -res flag averages over atoms within each residue.
# ============================================================

echo ">>> STEP E: RMSF"

printf '18\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_aptamer.xvg -res -xvg none

printf '17\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_Gquad.xvg -res -xvg none

printf '16\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_loops.xvg -res -xvg none

printf '13\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_loop1.xvg -res -xvg none

printf '14\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_loop2.xvg -res -xvg none

printf '15\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_loop3.xvg -res -xvg none

printf '8\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_Ttail.xvg -res -xvg none

# Whole-strand flexibility
printf '1\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_DNA.xvg -res -xvg none
 
printf '19\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_imm.xvg -res -xvg none
 
printf '20\n' | srun gmx_mpi rmsf \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o rmsf_TBAstrand.xvg -res -xvg none
    

# ============================================================
# STEP F: K+ ION ANALYSES
# Monitors whether K+ ions stay coordinated inside the G-quartet
# channel — loss of K+ correlates with G-quad unfolding.
# ============================================================

echo ">>> STEP F: K+ coordination"

# Distance from each K+ to Tetrad1_top COM (132 atoms, even — gmx distance works)
printf '10\n12\n' | srun gmx_mpi distance \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -oall kion_dist_tetrad1.xvg -xvg none

# Tetrad2_bot has 133 atoms (odd) so gmx distance fails; use mindist instead
printf '10\n11\n' | srun gmx_mpi mindist \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -od kion_dist_tetrad2.xvg -xvg none

# Number of O6 contacts within 3.5 Å (direct coordination criterion)
printf '10\n17\n' | srun gmx_mpi mindist \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -on kion_O6_contacts.xvg \
    -d 0.35 -xvg none

# Compute average K+ position between the two tetrad planes
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
# Tetrad2_bot has 133 atoms so gmx distance cannot compute its COM.
# gmx mindist gives the closest inter-tetrad atom distance instead,
# which is a reliable proxy for whether the planes stay stacked.
# ============================================================

echo ">>> STEP G: tetrad plane separation"

printf '11\n12\n' | srun gmx_mpi mindist \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -od gquad_tetrad_separation.xvg -xvg none


# ============================================================
# STEP H: RADIUS OF GYRATION
# ============================================================

echo ">>> STEP H: Rg"

printf '1\n'  | srun gmx_mpi gyrate -s ${TPR} -f ${TRAJ} -n ${NDX} -o gyrate_DNA.xvg   -xvg none
printf '17\n' | srun gmx_mpi gyrate -s ${TPR} -f ${TRAJ} -n ${NDX} -o gyrate_Gquad.xvg -xvg none
printf '8\n'  | srun gmx_mpi gyrate -s ${TPR} -f ${TRAJ} -n ${NDX} -o gyrate_Ttail.xvg -xvg none
printf '19\n' | srun gmx_mpi gyrate -s ${TPR} -f ${TRAJ} -n ${NDX} -o gyrate_imm.xvg        -xvg none
printf '20\n' | srun gmx_mpi gyrate -s ${TPR} -f ${TRAJ} -n ${NDX} -o gyrate_TBAstrand.xvg  -xvg none


# ============================================================
# STEP I: SASA
# ============================================================

echo ">>> STEP I: SASA"

printf '1\n'  | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_DNA.xvg    -xvg none
printf '17\n' | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_Gquad.xvg  -xvg none
printf '8\n'  | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_Ttail.xvg  -xvg none
printf '16\n' | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_loops.xvg  -xvg none
printf '13\n' | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_loop1.xvg  -xvg none
printf '14\n' | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_loop2.xvg  -xvg none
printf '15\n' | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_loop3.xvg  -xvg none
printf '20\n' | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_TBAstrand.xvg -xvg none
printf '19\n' | srun gmx_mpi sasa -s ${TPR} -f ${TRAJ} -n ${NDX} -o sasa_imm.xvg       -xvg none


# ============================================================
# STEP J: HYDROGEN BONDS (intra-aptamer structural integrity)
#
# There is no Protein group in this simulation, so protein-aptamer
# H-bonds cannot be measured here. What IS useful for your binding
# story:
#   - G-quad internal H-bonds: tell you whether the Hoogsteen base
#     pairs that define the G-quartet scaffold stay intact.
#   - Ttail-to-Gquad H-bonds: T-tail is the proposed binding
#     interface; losing contact with the core would affect how the
#     tail presents itself to the target.
#   - Loop-to-Gquad H-bonds: loops pack against the G-plex; their
#     stability shapes the binding surface geometry.
#
# NOTE: -ac (autocorrelation) was removed in GROMACS 2023+.
#       It is not used here.
# ============================================================

echo ">>> STEP J: hydrogen bonds"

printf '17\n17\n' | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_Gquad_internal.xvg \
    -xvg none

printf '8\n17\n' | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_Ttail_Gquad.xvg \
    -xvg none

printf '13\n17\n' | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_loop1_Gquad.xvg \
    -xvg none

printf '15\n17\n' | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_loop3_Gquad.xvg \
    -xvg none

printf '14\n17\n' | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_loop2_Gquad.xvg \
    -xvg none
 
# Duplex junction: TBA strand ↔ immobilized strand (Watson-Crick pairing)
printf '20\n19\n' | srun gmx_mpi hbond \
    -s ${TPR} -f ${TRAJ} -n ${NDX} \
    -num hbond_duplex.xvg \
    -xvg none


# ============================================================
# STEP M: PRINCIPAL COMPONENT ANALYSIS
#
# All PCA is done on traj_fit.xtc (G-plex fitted, full DNA output)
# so overall tumbling does not inflate the variance.
#
# covar takes TWO groups: fit group then analysis group.
# anaeig takes ONE group for projection.
#
# The convergence test (overlap between halves) checks whether
# the dominant motions in the first half are the same as in
# the second half — overlap near 1 means the simulation has
# sampled the same conformational space in both halves.
# ============================================================

echo ">>> STEP M: PCA"
#All DNA
printf '1\n1\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_DNA.xtc -n ${NDX} \
    -o eigenvalues_DNA.xvg \
    -v eigenvectors_DNA.trr \
    -av average_DNA.pdb \
    -xvg none

printf '1\n1\n' | srun gmx_mpi anaeig \
    -s average_DNA.pdb -f traj_fit_DNA.xtc -n ${NDX} \
    -v eigenvectors_DNA.trr \
    -eig eigenvalues_DNA.xvg \
    -proj pc_proj_DNA.xvg \
    -extr pc_extreme_DNA.pdb \
    -rmsf rmsf_pc_DNA.xvg \
    -first 1 -last 3 -xvg none

printf '1\n1\n' | srun gmx_mpi anaeig \
    -s average_DNA.pdb -f traj_fit_DNA.xtc -n ${NDX} \
    -v eigenvectors_DNA.trr \
    -eig eigenvalues_DNA.xvg \
    -2d pc2d_DNA.xvg \
    -first 1 -last 2 -xvg none
    
# Full aptamer
printf '18\n18\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o eigenvalues_apt.xvg \
    -v eigenvectors_apt.trr \
    -av average_apt.pdb \
    -xvg none

printf '18\n18\n' | srun gmx_mpi anaeig \
    -s average_apt.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_apt.trr \
    -eig eigenvalues_apt.xvg \
    -proj pc_proj_apt.xvg \
    -extr pc_extreme_apt.pdb \
    -rmsf rmsf_pc_apt.xvg \
    -first 1 -last 3 -xvg none

printf '18\n18\n' | srun gmx_mpi anaeig \
    -s average_apt.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_apt.trr \
    -eig eigenvalues_apt.xvg \
    -2d pc2d_apt.xvg \
    -first 1 -last 2 -xvg none

# T-tail
printf '8\n8\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o eigenvalues_Ttail.xvg \
    -v eigenvectors_Ttail.trr \
    -av average_Ttail.pdb \
    -xvg none

printf '8\n8\n' | srun gmx_mpi anaeig \
    -s average_Ttail.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_Ttail.trr \
    -eig eigenvalues_Ttail.xvg \
    -proj pc_proj_Ttail.xvg \
    -extr pc_extreme_Ttail.pdb \
    -rmsf rmsf_pc_Ttail.xvg \
    -first 1 -last 3 -xvg none

printf '8\n8\n' | srun gmx_mpi anaeig \
    -s average_Ttail.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_Ttail.trr \
    -eig eigenvalues_Ttail.xvg \
    -2d pc2d_Ttail.xvg \
    -first 1 -last 2 -xvg none

# G-quad core
printf '17\n17\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o eigenvalues_Gquad.xvg \
    -v eigenvectors_Gquad.trr \
    -av average_Gquad.pdb \
    -xvg none

printf '17\n17\n' | srun gmx_mpi anaeig \
    -s average_Gquad.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_Gquad.trr \
    -eig eigenvalues_Gquad.xvg \
    -proj pc_proj_Gquad.xvg \
    -extr pc_extreme_Gquad.pdb \
    -rmsf rmsf_pc_Gquad.xvg \
    -first 1 -last 3 -xvg none

printf '17\n17\n' | srun gmx_mpi anaeig \
    -s average_Gquad.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_Gquad.trr \
    -eig eigenvalues_Gquad.xvg \
    -2d pc2d_Gquad.xvg \
    -first 1 -last 2 -xvg none
    
# --- PCA: TBA_loops ---
printf '16\n16\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o eigenvalues_loops.xvg \
    -v eigenvectors_loops.trr \
    -av average_loops.pdb \
    -xvg none

printf '16\n16\n' | srun gmx_mpi anaeig \
    -s average_loops.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_loops.trr \
    -eig eigenvalues_loops.xvg \
    -proj pc_proj_loops.xvg \
    -extr pc_extreme_loops.pdb \
    -rmsf rmsf_pc_loops.xvg \
    -first 1 -last 3 -xvg none

printf '16\n16\n' | srun gmx_mpi anaeig \
    -s average_loops.pdb -f traj_fit.xtc -n ${NDX} \
    -v eigenvectors_loops.trr \
    -eig eigenvalues_loops.xvg \
    -2d pc2d_loops.xvg \
    -first 1 -last 2 -xvg none


    
# Cosine content: values close to 1 mean the motion is just noise/diffusion
srun gmx_mpi analyze -f pc_proj_loops.xvg -cc cosine_content_loops.xvg -xvg none
srun gmx_mpi analyze -f pc_proj_apt.xvg   -cc cosine_content_apt.xvg   -xvg none
srun gmx_mpi analyze -f pc_proj_Ttail.xvg -cc cosine_content_Ttail.xvg -xvg none
srun gmx_mpi analyze -f pc_proj_Gquad.xvg -cc cosine_content_Gquad.xvg -xvg none
srun gmx_mpi analyze -f pc_proj_DNA.xvg - cc cosine_content_DNA.xvg -xvg none

# Convergence: split trajectory in half, compare PCA overlap
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

printf '1\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o traj_fit_h1.xtc -e ${HALF_TIME}

printf '1\n' | srun gmx_mpi trjconv \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -o traj_fit_h2.xtc -b ${HALF_TIME}

# T-tail overlap
printf '8\n8\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_Ttail_h1.xvg -v eigenvectors_Ttail_h1.trr -xvg none
printf '8\n8\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_Ttail_h2.xvg -v eigenvectors_Ttail_h2.trr -xvg none
printf '8\n8\n' | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_Ttail_h1.trr \
    -eig eigenvalues_Ttail_h1.xvg \
    -over overlap_Ttail.xvg \
    -first 1 -last 10 -xvg none

# G-quad overlap
printf '17\n17\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_Gquad_h1.xvg -v eigenvectors_Gquad_h1.trr -xvg none
printf '17\n17\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_Gquad_h2.xvg -v eigenvectors_Gquad_h2.trr -xvg none
printf '17\n17\n' | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_Gquad_h1.trr \
    -eig eigenvalues_Gquad_h1.xvg \
    -over overlap_Gquad.xvg \
    -first 1 -last 10 -xvg none

# --- Loop overlap ---
printf '16\n16\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_loops_h1.xvg \
    -v eigenvectors_loops_h1.trr -xvg none

printf '16\n16\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_loops_h2.xvg \
    -v eigenvectors_loops_h2.trr -xvg none

printf '16\n16\n' | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_loops_h1.trr \
    -eig eigenvalues_loops_h1.xvg \
    -over overlap_loops.xvg \
    -first 1 -last 10 -xvg none
    
# --aptamer overlap ---    
printf '18\n18\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_apt_h1.xvg \
    -v eigenvectors_apt_h1.trr -xvg none

printf '18\n18\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_apt_h2.xvg \
    -v eigenvectors_apt_h2.trr -xvg none

printf '18\n18\n' | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_apt_h1.trr \
    -eig eigenvalues_apt_h1.xvg \
    -over overlap_apt.xvg \
    -first 1 -last 10 -xvg none
    
    
#DNA overlap
printf '1\n1\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h1.xtc -n ${NDX} \
    -o eigenvalues_DNA_h1.xvg \
    -v eigenvectors_DNA_h1.trr -xvg none

printf '1\n1\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -o eigenvalues_DNA_h2.xvg \
    -v eigenvectors_DNA_h2.trr -xvg none

printf '1\n1\n' | srun gmx_mpi anaeig \
    -s ${TPR} -f traj_fit_h2.xtc -n ${NDX} \
    -v eigenvectors_DNA_h1.trr \
    -eig eigenvalues_DNA_h1.xvg \
    -over overlap_DNA.xvg \
    -first 1 -last 10 -xvg none

    


# ============================================================
# STEP N: DYNAMIC CROSS-CORRELATION (DCCM)
# Writes the full covariance matrix as ASCII for plotting in
# Python (e.g. with numpy + matplotlib).  Uses the fitted
# trajectory and the full DNA group.
# ============================================================

echo ">>> STEP N: DCCM"

printf '1\n1\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit.xtc -n ${NDX} \
    -ascii dccm_apt.dat \
    -xvg none

# TBA aptamer (G-plex) only
printf '18\n18\n' | srun gmx_mpi covar \
    -s ${TPR} -f traj_fit_apt.xtc -n ${NDX} \
    -ascii dccm_apt.dat \
    -xvg none

# ============================================================
# DONE
# ============================================================

echo ""
echo "================================================"
echo " Analysis complete."
echo "================================================"
echo ""
echo "Sanity checks first:"
echo "  mindist_periodic.xvg       should stay above non-bonded cutoff (~1 nm)"
echo "  cosine_content_*.xvg       should be < 0.5 for meaningful PCA"
echo "  overlap_Ttail/Gquad.xvg    convergence: closer to 1 = better"
echo ""
echo "Key outputs for the report:"
echo "  rmsd_DNA/Gquad/Ttail/tetrads/loop*.xvg"
echo "  rmsf_Gquad/tetrads/loops/loop*/Ttail.xvg"
echo "  kion_dist_avg.xvg  kion_O6_contacts.xvg  gquad_tetrad_separation.xvg"
echo "  gyrate_*.xvg  sasa_*.xvg"
echo "  hbond_Gquad_internal.xvg  hbond_Ttail_Gquad.xvg  hbond_loop*_Gquad.xvg"
echo "  mindist_*_imm.xvg  interface_*_contact_frac.dat"
echo "  eigenvalues_*.xvg  pc2d_*.xvg  pc_extreme_*.pdb  dccm_apt.dat"
echo ""
echo "Now run: python3 post_analysis_TBA.py"