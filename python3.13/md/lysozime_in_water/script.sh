# Menghapus molekul air
grep -v HOH 1aki.pdb > 1AKI_clean.pdb


# menghasilkan:
# 
# Topologi molekul.
# File pembatasan posisi.
# File struktur yang telah diproses.
echo "Pilih 15 untuk force field"
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce


# Mendefinisikan Kotak
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic

# Mengisi Kotak dengan Pelarut
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top

# membuat file .tpr
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr

# Menambahkan Ion
echo "Pilih 13 untuk penempatan ion"
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral


# Menyusun Input untuk Minimasi Energi
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr

# Menjalankan Minimasi Energi
gmx mdrun -v -deffnm em

# Analisis Data Energi
echo "Pilih potential (10 0)"
gmx energy -f em.edr -o potential.xvg

# Menyusun Input untuk Equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

# jalankan simulasi equilibration 
gmx mdrun -deffnm nvt

# Analisis Progres Suhu
echo "Pilih suhu sistem dan keluar (16 0)"
gmx energy -f nvt.edr -o temperature.xvg

# Menyiapkan Simulasi NPT
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

# menjalankan simulasi
gmx mdrun -deffnm npt

# Analisis Perkembangan Tekanan
echo "Pilih suhu sistem dan keluar (18 0)"
gmx energy -f npt.edr -o pressure.xvg

# Analisis Perkembangan Kerapatan
echo "Pilih kerapatan sistem/density (24 0)"
gmx energy -f npt.edr -o density.xvg


# Menyiapkan Input untuk Simulasi MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

# Menjalankan Simulasi MD
gmx mdrun -deffnm md_0_1

# Jika menggunakan gpu jalankan dengan gpu
#gmx mdrun -deffnm md_0_1 -nb gpu

# Koreksi Periodisitas dan Pemrosesan Trajektori
echo "Pilih 1 ("Protein") untuk grup yang akan dipusatkan dan 0 ("System") untuk output."
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center


# Menghitung RMSD (Root Mean Square Deviation)
gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns

# Menghitung RMSD Relatif Terhadap Struktur Kristal
gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns

# Analisis Radius of Gyration (Rg)
echo "Pilih 1 ("Protein") untuk memilih grup protein yang akan dianalisis."
gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg