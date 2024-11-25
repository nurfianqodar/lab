# Lysozim dalam Air

## Langkah Pertama: Mempersiapkan Topologi

### Beberapa Dasar GROMACS

Mulai dari versi 5.0, semua alat di GROMACS menjadi bagian dari satu binari utama bernama `gmx`. Ini adalah perubahan dari versi sebelumnya, di mana setiap alat dijalankan sebagai perintah terpisah. Meskipun cara lama masih dapat digunakan melalui symlink, fitur tersebut akan dihapus pada versi mendatang. Oleh karena itu, lebih baik mulai membiasakan diri dengan cara baru ini.  

Untuk mendapatkan informasi bantuan tentang modul apa pun di GROMACS, gunakan salah satu dari perintah berikut:  

```
gmx help <modul>
```

atau  

```
gmx <modul> -h
```

Ganti `<modul>` dengan nama perintah yang ingin Anda gunakan. Informasi yang diberikan mencakup algoritma yang tersedia, opsi, format file yang diperlukan, bug yang diketahui, dan keterbatasan. Untuk pengguna baru, mempelajari perintah ini melalui bantuan bawaan adalah langkah awal yang baik.  


### Tutorial Lysozim  

Kita akan bekerja dengan struktur protein lysozim dari putih telur ayam (kode PDB: 1AKI). Berikut langkah-langkahnya:  

1. **Unduh File Struktur**  
   Kunjungi [situs RCSB](https://www.rcsb.org) dan unduh file teks PDB untuk struktur kristal lysozim dengan kode 1AKI.  

2. **Visualisasi Struktur**  
   Setelah diunduh, gunakan perangkat lunak seperti **VMD**, **Chimera**, atau **PyMOL** untuk melihat struktur protein.  

3. **Hapus Molekul Air Kristal**  
   Molekul air kristal (residu dengan label "HOH") harus dihapus dari file. Anda bisa melakukannya secara manual menggunakan editor teks seperti `vi`, `emacs`, atau Notepad.  
   Jangan gunakan perangkat lunak pengolah kata seperti Word karena dapat merusak format file!  
   Sebagai alternatif, gunakan perintah berikut di terminal:  

   ```
   grep -v HOH 1aki.pdb > 1AKI_clean.pdb
   ```

   > **Catatan Penting:**  
   > Prosedur ini tidak selalu cocok, terutama jika terdapat molekul air yang terikat erat atau berfungsi dalam situs aktif. Dalam kasus ini, molekul air tersebut perlu dipertahankan.  

4. **Verifikasi Struktur**  
   Periksa file `.pdb` Anda untuk entri yang ditandai dengan komentar `MISSING`.  
   - Jika ada atom atau residu yang hilang, Anda harus menambahkannya menggunakan perangkat lunak lain sebelum melanjutkan.  
   - Residu terminal yang hilang biasanya tidak menjadi masalah, tetapi residu internal yang tidak lengkap akan menyebabkan kegagalan pada langkah berikutnya.  

5. **Jalankan `pdb2gmx`**  
   File PDB yang sudah dibersihkan hanya harus berisi atom protein.  
   Langkah berikutnya adalah menjalankan perintah `pdb2gmx` untuk menghasilkan:  
   - Topologi molekul.  
   - File pembatasan posisi.  
   - File struktur yang telah diproses.  

   Gunakan perintah berikut:  

   ```
   gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
   ```

   Setelah perintah dijalankan, Anda akan diminta untuk memilih *force field*. Daftar pilihan akan muncul, seperti:  

   ```
   1: AMBER03  
   2: AMBER94  
   ...  
   15: OPLS-AA/L  
   ```

   Untuk tutorial ini, pilih opsi **15 (OPLS-AA/L)**.  

6. **File yang Dihasilkan**  
   Setelah menjalankan `pdb2gmx`, Anda akan mendapatkan tiga file baru:  
   - `1AKI_processed.gro`: File struktur dalam format GROMACS.  
   - `topol.top`: File topologi sistem.  
   - `posre.itp`: File pembatasan posisi atom berat.  

> **Catatan Tambahan:**  
> - Format `.gro` tidak wajib digunakan; GROMACS mendukung berbagai format file koordinat seperti `.pdb`.  
> - Tujuan utama `pdb2gmx` adalah menghasilkan topologi yang sesuai dengan *force field*. File struktur adalah hasil sampingan yang mempermudah pengguna.  

---

## Langkah Kedua: Memeriksa Topologi  

Mari kita periksa isi file topologi keluaran (`topol.top`). Gunakan editor teks untuk melihat isi file ini. Setelah beberapa baris komentar (diawali dengan `;`), Anda akan menemukan baris berikut:  

```
#include "oplsaa.ff/forcefield.itp"
```  

Baris ini memuat parameter dari force field OPLS-AA. Letaknya di awal file menunjukkan bahwa semua parameter selanjutnya berasal dari force field ini.  

### Bagian [ moleculetype ]  
Bagian berikutnya adalah definisi molekul:  

```
; Name       nrexcl
Protein_A    3
```  

Nama `Protein_A` mendefinisikan molekul berdasarkan label rantai A dalam file PDB. Angka `3` menunjukkan jumlah pengecualian untuk tetangga yang terikat (lihat manual GROMACS untuk detail lebih lanjut).  

### Bagian `[ atoms ]`  
Bagian ini mendeskripsikan atom dalam protein. Informasi ditampilkan dalam kolom seperti berikut:  

```
[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
; residue   1 LYS rtp LYSH q +2.0
     1   opls_287      1   LYS       N      1       -0.3    14.0067   ; qtot -0.3
     2   opls_290      1   LYS      H1      1       0.33      1.008   ; qtot 0.03
     3   opls_290      1   LYS      H2      1       0.33      1.008   ; qtot 0.36
     4   opls_290      1   LYS      H3      1       0.33      1.008   ; qtot 0.69
     5  opls_293B      1   LYS      CA      1       0.25     12.011   ; qtot 0.94
     6   opls_140      1   LYS      HA      1       0.06      1.008   ; qtot 1
```  

Penjelasan kolom:  
- `nr`: Nomor atom.  
- `type`: Jenis atom.  
- `resnr`: Nomor residu asam amino.  
- `residue`: Nama residu asam amino. Contoh: `LYS` dalam file PDB menunjukkan residu "LYSH" yang terprotonasi pada pH netral.  
- `atom`: Nama atom.  
- `cgnr`: Nomor grup muatan. Grup muatan membantu mempercepat perhitungan.  
- `charge`: Muatan atom.  
- `mass`: Massa atom.  
- `typeB`, `chargeB`, `massB`: Parameter untuk pertumbuhan energi bebas (tidak dibahas di sini).  

### Bagian Lainnya  
Bagian berikutnya mencakup definisi:  
- `[ bonds ]`: Ikatan antaratom.  
- `[ pairs ]`: Interaksi 1-4 khusus (lihat manual GROMACS, bagian 5.3.4).  
- `[ angles ]`: Sudut antaratom.  
- `[ dihedrals ]`: Sudut dihedral antaratom.  

Parameter dan tipe fungsi untuk bagian-bagian ini dijelaskan lebih rinci di Bab 5 manual GROMACS.  

### Penambahan File `posre.itp`  
File `posre.itp` dihasilkan oleh `pdb2gmx` dan digunakan untuk menahan posisi atom selama proses ekuilibrasi. Tambahan file ini diatur seperti berikut:  

```
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif
```  

### Definisi Molekul Lainnya  
Setelah molekul "Protein_A", file topologi juga mendefinisikan pelarut, yaitu air model SPC/E:  

```
; Include water topology
#include "oplsaa.ff/spce.itp"
```  

Pelarut dapat diberi restriksi posisi dengan konstanta gaya 1000 kJ mol⁻¹ nm⁻² seperti berikut:  

```
#ifdef POSRES_WATER
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif
```  

### Parameter Ion  
File ini juga mencakup parameter ion untuk simulasi:  

```
; Include generic topology for ions
#include "oplsaa.ff/ions.itp"
```  

### Definisi Sistem  
Bagian akhir file mencakup deskripsi level sistem:  

```
[ system ]
; Name
LYSOZYME

[ molecules ]
; Compound        #mols
Protein_A           1
```  

>**Catatan Penting:**  
>- Urutan molekul dalam bagian `[ molecules ]` harus sesuai dengan urutan molekul dalam file koordinat (misalnya, `.gro`).  
>- Nama yang tercantum harus cocok dengan nama `[ moleculetype ]`.  
Apabila aturan ini tidak terpenuhi, Anda akan mendapatkan kesalahan fatal dari perintah `grompp`.  

Dengan memahami isi file topologi ini, kita dapat melanjutkan untuk membangun sistem simulasi.

---

## Langkah Ketiga: Mendefinisikan Sel Satuan & Menambahkan Pelarut  

Setelah memahami isi file topologi GROMACS, kita dapat melanjutkan untuk membangun sistem simulasi. Dalam contoh ini, kita akan melakukan simulasi sistem akuatis sederhana. Anda bisa mensimulasikan protein dan molekul lainnya dalam pelarut yang berbeda, asalkan parameter yang sesuai tersedia untuk semua spesies yang terlibat.

Ada dua langkah utama untuk mendefinisikan kotak dan mengisinya dengan pelarut:

1. Mendefinisikan dimensi kotak menggunakan modul `editconf`.
2. Mengisi kotak dengan air menggunakan modul `solvate` (sebelumnya disebut `genbox`).

Sekarang Anda akan dihadapkan dengan pilihan tentang bagaimana menangani sel satuan. Untuk tutorial ini, kita akan menggunakan kotak kubik sebagai sel satuan. Setelah Anda lebih memahami kondisi batas periodik dan jenis kotak, saya sangat merekomendasikan penggunaan dodekahedron rombik, karena volumenya sekitar 71% dari kotak kubik dengan jarak periodik yang sama, sehingga mengurangi jumlah molekul air yang perlu ditambahkan untuk melarutkan protein.

### Mendefinisikan Kotak dengan `editconf`  
Perintah untuk mendefinisikan kotak adalah sebagai berikut:

```bash
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic
```

Penjelasan parameter:
- `-f 1AKI_processed.gro`: Menentukan file koordinat protein yang telah diproses.
- `-o 1AKI_newbox.gro`: Menentukan file keluaran dengan koordinat protein yang sudah dimasukkan ke dalam kotak.
- `-c`: Memusatkan protein di dalam kotak.
- `-d 1.0`: Memberikan jarak minimal 1.0 nm dari tepi kotak.
- `-bt cubic`: Menentukan jenis kotak sebagai kubik.

Parameter jarak ke tepi kotak penting karena kita akan menggunakan kondisi batas periodik. Dalam hal ini, jarak minimal antar citra periodik protein adalah 2.0 nm, yang cukup untuk hampir semua skema cutoff yang umum digunakan dalam simulasi.

### Mengisi Kotak dengan Pelarut menggunakan `solvate`  
Setelah mendefinisikan kotak, langkah berikutnya adalah mengisi kotak tersebut dengan pelarut. Proses pelarutan dilakukan dengan perintah `solvate`:

```
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```

Penjelasan parameter:
- `-cp 1AKI_newbox.gro`: Menggunakan konfigurasi protein yang dihasilkan dari langkah `editconf` sebelumnya.
- `-cs spc216.gro`: Menggunakan konfigurasi pelarut air SPC/E yang sudah ada dalam instalasi GROMACS (spc216.gro adalah model pelarut tiga titik yang sudah teruji).
- `-o 1AKI_solv.gro`: Menentukan nama file keluaran yang berisi protein yang telah terlarut.
- `-p topol.top`: Menentukan file topologi (`topol.top`) agar bisa diperbarui sesuai dengan jumlah molekul air yang ditambahkan.

Output dari perintah ini adalah file `1AKI_solv.gro` yang berisi protein yang telah dilarutkan, serta file topologi yang telah diperbarui. Di dalam file topologi, Anda akan melihat perubahan pada bagian `[ molecules ]`, yang mencatat jumlah molekul pelarut yang telah ditambahkan:

```
[ molecules ]
; Compound  #mols
Protein_A       1 
SOL         10832
```

Perintah `solvate` secara otomatis memperbarui jumlah molekul air yang ditambahkan ke dalam sistem, dan mencatat perubahan tersebut di dalam file topologi. Perhatikan bahwa jika Anda menggunakan pelarut selain air (misalnya, pelarut non-air), perintah `solvate` tidak akan memperbarui topologi Anda. Fitur ini hanya berfungsi untuk molekul air karena kompatibilitas pembaruan molekul air yang telah diprogram sebelumnya.

---

## Langkah Keempat: Menambahkan Ion

Sekarang kita memiliki sistem yang telah dilarutkan dan mengandung protein bermuatan. Output dari `pdb2gmx` memberi tahu kita bahwa protein memiliki muatan total +8e (berdasarkan komposisi asam amino). Jika Anda melewatkan informasi ini di output `pdb2gmx`, Anda bisa melihatnya pada baris terakhir dari direktif `[ atoms ]` dalam file `topol.top` yang seharusnya mencatat "qtot 8". Karena kehidupan tidak dapat ada pada muatan bersih, kita perlu menambahkan ion ke dalam sistem.

Untuk menambahkan ion dalam GROMACS, kita menggunakan alat `genion`. Yang dilakukan oleh `genion` adalah membaca file topologi dan mengganti molekul air dengan ion yang Anda tentukan. Input yang digunakan oleh `genion` adalah file input run dengan ekstensi `.tpr`, yang dihasilkan oleh modul GROMACS `grompp` (GROMACS pre-processor), yang juga akan digunakan nanti saat kita menjalankan simulasi pertama kita. Apa yang dilakukan `grompp` adalah memproses file koordinat dan topologi (yang mendeskripsikan molekul) untuk menghasilkan file input tingkat atom `.tpr`. File `.tpr` ini berisi semua parameter untuk atom-atom dalam sistem.

### Membuat File `.tpr` dengan `grompp`
Untuk menghasilkan file `.tpr` dengan `grompp`, kita memerlukan file input tambahan dengan ekstensi `.mdp` (molecular dynamics parameter file). `grompp` akan menggabungkan parameter yang ditentukan dalam file `.mdp` dengan informasi koordinat dan topologi untuk menghasilkan file `.tpr`.

File `.mdp` umumnya digunakan untuk menjalankan minimisasi energi atau simulasi MD, tetapi dalam kasus ini, kita menggunakannya untuk menghasilkan deskripsi atom dari sistem. 
> Anda bisa mengunduh contoh file `.mdp` (yang akan kita gunakan) dari sini [ions.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp)

Perlu dicatat bahwa file `.mdp` yang digunakan pada langkah ini bisa berisi kombinasi parameter yang sah lainnya. Biasanya, saya menggunakan skrip minimisasi energi, karena ini sangat dasar dan tidak melibatkan kombinasi parameter yang rumit. Perhatikan juga bahwa file yang disediakan dalam tutorial ini hanya untuk digunakan dengan gaya gaya OPLS-AA. Pengaturan, khususnya pengaturan interaksi non-bonded, akan berbeda jika menggunakan gaya gaya lain.

Sekarang kita bisa membuat file `.tpr` dengan perintah berikut:

```
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
```

Sekarang kita memiliki deskripsi tingkat atom dari sistem kita dalam file biner `ions.tpr`. Kita akan memberikan file ini sebagai input ke `genion`.

### Menambahkan Ion dengan `genion`
Perintah untuk menambahkan ion adalah sebagai berikut:

```
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

Penjelasan parameter:
- `-s ions.tpr`: Memberikan file struktur/sistem yang telah diproses.
- `-o 1AKI_solv_ions.gro`: Menentukan file keluaran yang akan berisi struktur dengan ion yang telah ditambahkan.
- `-p topol.top`: Memproses file topologi untuk mencerminkan penghapusan molekul air dan penambahan ion.
- `-pname NA`: Menentukan nama ion positif (natrium, Na+).
- `-nname CL`: Menentukan nama ion negatif (klorida, Cl-).
- `-neutral`: Memberitahu `genion` untuk menambahkan ion yang diperlukan untuk menetralkan muatan total protein. Dalam hal ini, 8 ion Cl- akan ditambahkan untuk mengimbangi muatan +8 pada protein.

Saat diminta, pilih grup 13 "SOL" untuk menempatkan ion. Anda tidak ingin mengganti bagian dari protein dengan ion.

### Pembaruan pada Direktif `[ molecules ]`
Setelah proses `genion`, direktif `[ molecules ]` dalam file `topol.top` akan terlihat seperti ini:

```
[ molecules ]
; Compound      #mols
Protein_A         1
SOL               10636
CL                8
```

Di sini, jumlah molekul air (SOL) telah diperbarui menjadi 10.636 (dengan beberapa molekul air yang digantikan oleh ion). Ion klorida (CL) sebanyak 8 molekul ditambahkan untuk menetralkan muatan protein.

---

## Langkah Kelima: Minimasi Energi

Sistem yang sudah dilarutkan dan netral secara elektrik sekarang sudah ter组. Sebelum kita memulai simulasi dinamis, kita harus memastikan bahwa sistem tidak memiliki benturan sterik atau geometri yang tidak sesuai. Struktur sistem ini akan direlaksasi melalui proses yang disebut **minimasi energi (EM)**.

Proses minimasi energi mirip dengan penambahan ion. Kali ini, kita akan menggunakan `grompp` untuk menggabungkan struktur, topologi, dan parameter simulasi menjadi sebuah file input biner (.tpr), namun kali ini, alih-alih memberikan file .tpr ke `genion`, kita akan menjalankan minimasi energi menggunakan mesin MD GROMACS, yaitu `mdrun`.

### Menyusun Input untuk Minimasi Energi
>Gunakan file parameter input untuk EM berikut: [minim.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp)

```
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
```

Pastikan Anda telah memperbarui file `topol.top` saat menjalankan `genbox` dan `genion`, jika tidak, Anda akan menerima banyak pesan kesalahan ("jumlah koordinat dalam file koordinat tidak cocok dengan topologi", dll).

### Menjalankan Minimasi Energi
Sekarang kita siap untuk menjalankan `mdrun` untuk melaksanakan minimasi energi:

```
gmx mdrun -v -deffnm em
```

Penjelasan parameter:
- `-v`: Membuat `mdrun` lebih verbose, sehingga ia menampilkan progresnya di layar pada setiap langkah.
- `-deffnm em`: Menentukan nama file input dan output. Jika Anda tidak memberi nama file output `grompp` dengan nama "em.tpr", Anda harus secara eksplisit menyebutkan nama file dengan flag `-s` di `mdrun`.

Setelah menjalankan perintah ini, Anda akan mendapatkan file berikut:
- `em.log`: File log teks ASCII dari proses EM.
- `em.edr`: File biner energi.
- `em.trr`: File biner trajektori presisi penuh.
- `em.gro`: Struktur setelah minimasi energi.

### Evaluasi Hasil Minimasi Energi
Ada dua faktor penting untuk mengevaluasi apakah EM berhasil. Yang pertama adalah **energi potensial (Epot)**, yang dicetak di akhir proses EM, bahkan tanpa menggunakan flag `-v`. **Epot** harus negatif, dan untuk protein sederhana dalam air, biasanya berada pada kisaran 10^5-10^6, tergantung pada ukuran sistem dan jumlah molekul air.

Faktor penting kedua adalah **gaya maksimum (Fmax)**, yang targetnya telah ditentukan dalam file **minim.mdp** dengan parameter "emtol = 1000.0", yang berarti target Fmax tidak boleh lebih besar dari 1000 kJ mol-1 nm-1. Terkadang, meskipun **Epot** sudah masuk akal, **Fmax** bisa lebih besar dari target **emtol**. Jika ini terjadi, sistem Anda mungkin belum cukup stabil untuk simulasi. Anda bisa mengevaluasi penyebabnya dan mencoba mengganti parameter minimisasi (misalnya `integrator`, `emstep`, dll).

### Analisis Data Energi
File `em.edr` berisi semua istilah energi yang dikumpulkan GROMACS selama EM. Anda bisa menganalisis file `.edr` menggunakan modul energi GROMACS:

```
gmx energy -f em.edr -o potential.xvg
```

Di prompt, ketik `10 0` untuk memilih `Potential (10)`, dan ketik `0` untuk menghentikan input. Anda akan melihat rata-rata `Epot`, dan file `potential.xvg` akan ditulis. Untuk memplot data ini, Anda membutuhkan alat plotting Xmgrace. Grafik yang dihasilkan seharusnya menunjukkan konvergensi yang stabil dari `Epot`, yang menandakan bahwa sistem telah mencapai energi minimum yang diinginkan.

### Lanjut ke Simulasi Dinamis
Sekarang bahwa sistem kita berada dalam energi minimum, kita siap untuk memulai simulasi dinamis nyata.

---

## Langkah Keenam: Pembangunan Keseimbangan (Equilibration)

Minimasi energi (EM) memastikan bahwa kita memiliki struktur awal yang wajar dalam hal geometri dan orientasi pelarut. Untuk memulai dinamika yang sesungguhnya, kita perlu mengelola pelarut dan ion di sekitar protein. Jika kita mencoba melakukan dinamika tanpa pembatas pada tahap ini, sistem bisa saja runtuh. Alasan utamanya adalah bahwa pelarut sudah dioptimalkan untuk dirinya sendiri, dan bukan dengan solutnya (protein). Pelarut perlu disesuaikan dengan suhu yang ingin kita simulasikan dan membangun orientasi yang tepat di sekitar solut (protein). Setelah kita mencapai suhu yang diinginkan (berdasarkan energi kinetik), kita akan memberikan tekanan pada sistem hingga mencapai kerapatan yang benar.

### Menggunakan File Posisi Restrains (`posre.itp`)
Ingat file `posre.itp` yang dihasilkan oleh `pdb2gmx`? Kita akan menggunakannya sekarang. Tujuan dari `posre.itp` adalah untuk menerapkan gaya pembatas posisi pada atom-atom berat protein (semua atom yang bukan hidrogen). Pergerakan diizinkan, tetapi hanya setelah melewati hukuman energi yang substansial. Kegunaan dari pembatas posisi adalah memungkinkan kita menyeimbangkan pelarut di sekitar protein tanpa menambah variabel berupa perubahan struktur protein. Asal mula pembatas posisi (koordinat pada mana potensi pembatas nol) disediakan melalui file koordinat yang diteruskan ke opsi `-r` pada `grompp`.

### Fase Equilibration
Proses equilibration sering dilakukan dalam dua fase. Fase pertama dilakukan dalam ensemble **NVT** (konstanta jumlah partikel, volume, dan suhu). Ensemble ini juga disebut sebagai **isothermal-isochoric** atau **kanonik**. Waktu yang dibutuhkan untuk prosedur ini bergantung pada isi sistem, namun dalam **NVT**, suhu sistem harus mencapai nilai rata-rata yang diinginkan. Jika suhu belum stabil, waktu tambahan akan diperlukan. Biasanya, 50-100 ps sudah cukup, dan untuk latihan ini, kita akan melakukan **100 ps NVT equilibration**. Bergantung pada mesin yang digunakan, ini mungkin memerlukan waktu beberapa waktu (sekitar satu jam jika dijalankan secara paralel di 16 core).

### Menyusun Input untuk Equilibration
Sama seperti pada langkah minimasi energi, kita akan menggunakan `grompp` dan `mdrun` untuk menyusun dan menjalankan simulasi. Pertama, kita menggunakan `grompp` untuk menghasilkan file input `.tpr`:
>Dapatkan `.mdp` file disini [nvt.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp)

```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```

Selanjutnya, jalankan simulasi equilibration menggunakan `mdrun`:

```
gmx mdrun -deffnm nvt
```

Beberapa parameter penting yang digunakan dalam file `nvt.mdp`:
- `gen_vel = yes`: Menghasilkan kecepatan acak untuk atom-atom dalam sistem. Dengan menggunakan biji acak yang berbeda (`gen_seed`), kita bisa mendapatkan kecepatan awal yang berbeda dan dengan demikian menjalankan simulasi yang berbeda dari struktur awal yang sama.
- `tcoupl = V-rescale`: Termostat pengaturan kecepatan rescaling, yang merupakan peningkatan dari metode penggabungan lemah Berendsen, yang sebelumnya tidak menghasilkan ensemble kinetik yang benar.
- `pcoupl = no`: Tidak ada pengaturan tekanan yang diterapkan pada fase ini karena kita hanya ingin mengelola suhu dalam fase NVT.

### Analisis Progres Suhu
Setelah equilibration, kita bisa menganalisis perubahan suhu selama proses tersebut. Gunakan modul energi GROMACS untuk memeriksa file `nvt.edr`:

```
gmx energy -f nvt.edr -o temperature.xvg
```

Pada prompt, pilih `16 0` untuk memilih suhu sistem dan keluar. Anda akan melihat grafik yang menunjukkan progres suhu sistem. Grafik tersebut seharusnya menunjukkan bahwa suhu sistem cepat mencapai nilai target (300 K), dan tetap stabil selama periode equilibration.

### Hasil yang Diharapkan
Grafik suhu yang dihasilkan harus terlihat seperti ini:

- Suhu sistem cepat naik ke 300 K.
- Suhu tetap stabil selama sisa periode equilibration.
- Setelah 50-100 ps, sistem Anda harus stabil dan siap untuk dilanjutkan ke fase equilibration tekanan (NPT) atau simulasi dinamis yang lebih lanjut.

Dengan sistem yang sudah berada dalam keadaan keseimbangan suhu, kita siap untuk melanjutkan ke langkah berikutnya, yaitu **Equilibration dengan tekanan (NPT)** atau simulasi dinamis jangka panjang.

---

## Langkah Ketujuh: Pembangunan Keseimbangan, Bagian 2 (Equilibration NPT)

Pada langkah sebelumnya, **NVT equilibration** telah menstabilkan suhu sistem. Sebelum kita mulai mengumpulkan data, kita juga perlu menstabilkan tekanan (dan kerapatan) sistem. Keseimbangan tekanan dilakukan dalam ensemble **NPT**, di mana jumlah partikel, tekanan, dan suhu semua tetap konstan. Ensemble ini juga disebut sebagai **isothermal-isobaric** dan paling mendekati kondisi eksperimen.

### Menggunakan File Parameter untuk NPT Equilibration
File `.mdp` yang digunakan untuk **100-ps NPT equilibration** dapat ditemukan dalam referensi tutorial ini. File ini tidak terlalu berbeda dari file parameter yang digunakan untuk **NVT equilibration**, dengan beberapa perubahan penting:
- `continuation = yes`: Menandakan bahwa simulasi akan dilanjutkan dari fase NVT equilibration.
- `gen_vel = no`: Kecepatan diambil dari trajektori (lihat penjelasan di bawah).
- `Pressure Coupling`: Pada file NPT, kita menambahkan bagian `pressure coupling`, menggunakan **Parrinello-Rahman barostat** untuk mengatur tekanan sistem.

### Menyiapkan Simulasi NPT
Untuk memulai fase equilibration NPT, kita akan memanggil `grompp` dan `mdrun` seperti pada langkah NVT, tetapi kali ini kita juga menambahkan opsi `-t` untuk menyertakan file checkpoint dari fase NVT, yang berisi variabel status untuk melanjutkan simulasi. File koordinat (`-c`) yang digunakan adalah output akhir dari simulasi NVT.
> Dapatkan `.mdp` file disini [npt.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp)
Perintah untuk menyiapkan simulasi NPT adalah sebagai berikut:

```
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
```

Kemudian, untuk menjalankan simulasi:

```
gmx mdrun -deffnm npt
```

### Analisis Perkembangan Tekanan
Setelah simulasi NPT berjalan, kita dapat menganalisis perkembangan tekanan dengan menggunakan GROMACS `energy` module:

```
gmx energy -f npt.edr -o pressure.xvg
```

Pada prompt, ketik "18 0" untuk memilih tekanan sistem dan keluar. Grafik yang dihasilkan akan menunjukkan fluktuasi tekanan yang signifikan selama fase equilibration 100 ps. Namun, fluktuasi ini tidaklah mengejutkan. Rata-rata pergerakan data ini digambarkan dengan garis merah pada grafik. Rata-rata tekanan sistem selama equilibration adalah **7.5 ± 160.5 bar**. Mengingat tekanan referensi disetel ke **1 bar**, hasil ini masih dapat diterima. Tekanan dalam simulasi MD memang cenderung berfluktuasi secara luas, seperti yang terlihat dari fluktuasi besar dalam nilai deviasi standar (160.5 bar), sehingga secara statistik, perbedaan antara rata-rata (7.5 ± 160.5 bar) dan nilai referensi (1 bar) tidak dapat dibedakan.

### Analisis Perkembangan Kerapatan
Selain tekanan, kita juga perlu menganalisis **kerapatan** sistem dengan menggunakan perintah yang sama:

```
gmx energy -f npt.edr -o density.xvg
```

Masukkan "24 0" di prompt untuk memilih kerapatan. Grafik yang dihasilkan akan menunjukkan fluktuasi kerapatan yang stabil, dengan rata-rata **1019 ± 3 kg/m³**, yang sangat mendekati nilai eksperimen untuk air (1000 kg/m³) dan nilai kerapatan model SPC/E air (1008 kg/m³). Ini menunjukkan bahwa sistem sudah sangat stabil dalam hal kerapatan dan tekanan.

>**Catatan Penting** 
>
>Fluktuasi terkait tekanan dan kerapatan cenderung membutuhkan waktu konvergensi yang lebih lama, jadi Anda mungkin perlu menjalankan fase equilibration NPT sedikit lebih lama dari yang disarankan di sini untuk memastikan sistem sepenuhnya seimbang. Jika hasil Anda tidak sepenuhnya cocok dengan nilai eksperimental atau yang diharapkan, pertimbangkan untuk memperpanjang waktu equilibration atau mengevaluasi parameter lainnya.

Dengan sistem yang telah terimbang untuk suhu, tekanan, dan kerapatan, Anda siap untuk melanjutkan ke simulasi dinamis yang lebih panjang untuk memulai pengumpulan data atau analisis lebih lanjut.

---

## Langkah Kedelapan: Simulasi Dinamis Produksi (Production MD)

Setelah kedua fase equilibration selesai, sistem kini telah terimbang dengan baik pada suhu dan tekanan yang diinginkan. Sekarang kita siap untuk melepaskan posisi restrain dan menjalankan simulasi dinamis produksi untuk pengumpulan data. Proses ini serupa dengan tahap-tahap sebelumnya, di mana kita akan menggunakan file checkpoint yang sekarang sudah berisi informasi pengaturan tekanan (pressure coupling) untuk diserahkan ke `grompp`. 

Pada tahap ini, kita akan menjalankan simulasi MD selama 1 ns. Berikut adalah cara melanjutkan langkah-langkah tersebut:

### 1. Menyiapkan Input untuk Simulasi MD
Pertama, kita akan menyiapkan file `.tpr` untuk simulasi MD dengan menggunakan perintah `grompp`. Perintah ini akan memproses file `md.mdp` (file parameter untuk MD) bersama file koordinat dan topology untuk menghasilkan file input `md_0_1.tpr`:
> file parameter untuk MD didapatkan disini [md.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp)

```
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
```

Perintah ini akan mencetak estimasi beban komputasi untuk **PME** (Particle Mesh Ewald) yang menunjukkan bagaimana pembagian pemrosesan untuk **PME** dan **PP** (Particle-Particle) dilakukan. Estimasi ini akan memberi petunjuk berapa banyak core prosesor yang harus didedikasikan untuk masing-masing jenis perhitungan.

**Contoh estimasi untuk beban komputasi**:
```
Estimate for the relative computational load of the PME mesh part: 0.22
```

Dengan nilai ini, kita dapat mengetahui apakah pembagian pemrosesan antara **PME** dan **PP** sudah optimal. Sebagai referensi, untuk kotak kubik, pembagian optimal adalah 3:1 (PP:PME), sementara untuk kotak dodekhedral, pembagian optimalnya adalah 2:1 (PP:PME). Setelah Anda memiliki file input `md_0_1.tpr`, kita siap untuk menjalankan simulasi.

### 2. Menjalankan Simulasi MD
Untuk menjalankan simulasi MD, gunakan perintah `mdrun`:

```
gmx mdrun -deffnm md_0_1
```

Simulasi ini akan dijalankan untuk waktu 1 ns sesuai dengan parameter yang ditentukan di file `md.mdp`. GROMACS akan otomatis memilih pembagian yang optimal antara proses perhitungan **PME** dan **PP** jika Anda menggunakan banyak core prosesor.

### 3. Menggunakan GPU untuk Percepatan
Jika Anda menggunakan GROMACS versi 2018 atau yang lebih baru, dan memiliki GPU yang kompatibel, Anda bisa memanfaatkan akselerasi GPU untuk mempercepat perhitungan simulasi. Misalnya, pada GPU Titan Xp, simulasi ini dapat berjalan dengan kecepatan luar biasa, mencapai **295 ns/hari**.

Untuk memanfaatkan GPU, cukup tambahkan flag `-nb gpu` saat menjalankan `mdrun`:

```
gmx mdrun -deffnm md_0_1 -nb gpu
```

Jika Anda memiliki lebih dari satu GPU, atau perlu mengkustomisasi pembagian pekerjaan dengan skema paralelisasi hibrida, Anda dapat merujuk ke manual GROMACS untuk informasi lebih lanjut mengenai konfigurasi tersebut. Dalam tutorial ini, kami akan fokus pada penggunaan satu GPU.

### 4. Melanjutkan Simulasi dengan Banyak Core
Untuk hasil yang lebih optimal dalam hal kecepatan, Anda juga bisa menggunakan banyak core dengan cara mengatur jumlah core yang digunakan melalui flag `-nt X` (di mana X adalah jumlah core yang digunakan).

### Menyimpulkan
Setelah langkah ini, Anda akan memperoleh hasil simulasi MD yang stabil, dengan file output seperti trajektori, file energi, dan data lainnya untuk dianalisis lebih lanjut. Anda bisa mengamati perubahan struktur, energi, atau parameter lainnya yang relevan untuk studi Anda.

---

## Langkah Kesembilan: Analisis

Sekarang setelah kita menjalankan simulasi untuk protein, kita perlu melakukan beberapa analisis pada sistem. Pertanyaan penting yang perlu diajukan adalah jenis data apa yang penting untuk dikumpulkan, karena ini akan membantu merencanakan analisis dan pemrosesan data lebih lanjut. Di tutorial ini, beberapa alat dasar akan diperkenalkan untuk analisis post-simulasi.

### 1. Koreksi Periodisitas dan Pemrosesan Trajektori
Salah satu alat penting dalam GROMACS untuk memanipulasi trajektori adalah `trjconv`. Alat ini digunakan untuk menghilangkan koordinat tertentu, memperbaiki masalah periodisitas, atau mengubah unit waktu dan frekuensi frame.

Selama simulasi, protein mungkin mengalami difusi melalui unit sel, yang dapat menyebabkan penampilan "terputus" atau melompat ke sisi lain dari kotak periodik. Untuk mengatasi masalah ini, kita akan menggunakan `trjconv` untuk memperbaiki periodisitas dan menjaga molekul tetap utuh di sepanjang trajektori. Gunakan perintah berikut:

```
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
```

Pada prompt, pilih `1 ("Protein")` untuk grup yang akan dipusatkan dan `0 ("System")` untuk output. Hasilnya adalah trajektori yang telah dikoreksi untuk periodisitas, yang siap untuk analisis lebih lanjut.

### 2. Menghitung RMSD (Root Mean Square Deviation)
Setelah trajektori dikoreksi, kita dapat menganalisis stabilitas struktural protein menggunakan RMSD (Root Mean Square Deviation). RMSD digunakan untuk mengukur seberapa jauh posisi atom dalam struktur selama simulasi berbeda dari struktur referensi. Dalam hal ini, kita akan menghitung RMSD terhadap struktur protein yang telah diselesaikan dalam fase minimisasi dan equilibration.

Untuk menghitung RMSD, gunakan perintah berikut:

```
gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
```

Pilih `4 ("Backbone")` sebagai grup yang akan digunakan untuk perhitungan RMSD, yang mewakili kerangka utama protein. Flag `-tu ns` digunakan untuk mengkonversi waktu ke satuan nanodetik (ns), meskipun trajektori disimpan dalam pikodetik (ps). Ini memudahkan pembacaan hasil, terutama untuk simulasi panjang.

Plot output dari perintah ini akan menunjukkan pergerakan RMSD relatif terhadap struktur sistem yang telah dipersiapkan dan ter-minimisasi sebelumnya. Jika protein stabil, Anda akan melihat RMSD yang relatif kecil atau tetap (misalnya, di bawah 1 Å), menunjukkan bahwa strukturnya tidak berubah secara signifikan selama simulasi.

### 3. Menghitung RMSD Relatif Terhadap Struktur Kristal
Jika Anda ingin menghitung RMSD relatif terhadap struktur kristal atau struktur referensi lain, Anda bisa mengubah file input menjadi `em.tpr` (struktur hasil minimisasi energi) dan menjalankan perintah berikut:

```
gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
```

Plot untuk RMSD ini akan memberikan perbandingan antara struktur simulasi dan struktur kristal. Dengan membandingkan kedua plot ini, Anda bisa mendapatkan gambaran tentang kestabilan struktur dan perbedaan kecil yang mungkin terjadi antara struktur yang ter-minimisasi dan struktur kristal.

### 4. Interpretasi Hasil RMSD
Jika hasil RMSD menunjukkan bahwa RMSD mencapai nilai konstan (misalnya, ~0.1 nm atau 1 Å), ini menunjukkan bahwa struktur protein sangat stabil selama simulasi. Perbedaan kecil antara hasil RMSD terhadap struktur yang telah minimisasi energi dan struktur kristal dapat terjadi, tetapi perbedaan tersebut biasanya tidak signifikan dan dapat dijelaskan oleh perubahan energi atau ketidaksempurnaan dalam restrain posisi selama fase equilibration.

### Langkah Berikutnya
Dengan menggunakan alat analisis ini, Anda dapat memperoleh wawasan tentang stabilitas struktural sistem Anda. Analisis RMSD adalah salah satu langkah awal dalam mengevaluasi apakah sistem telah mencapai keadaan stabil, dan apakah protein berperilaku sesuai dengan ekspektasi berdasarkan struktur kristal atau minimisasi energi yang sebelumnya.

---

## Langkah Kesepuluh: Analisis Radius of Gyration (Rg)

Radius of gyration (**Rg**) adalah ukuran dari kepadatan atau kompaktnya suatu protein. Untuk protein yang terlipat stabil, nilai Rg-nya cenderung stabil sepanjang waktu. Sebaliknya, jika protein tersebut terurai atau tidak terlipat, nilai Rg akan menunjukkan perubahan yang signifikan.

Dalam konteks simulasi ini, kita akan menganalisis **Rg** dari lysozyme untuk memastikan apakah protein tetap terlipat selama simulasi atau tidak.

### Langkah-langkah untuk Menghitung Rg

Untuk menghitung **Rg** menggunakan GROMACS, kita akan menggunakan utilitas `gmx gyrate`, yang dirancang untuk menganalisis radius of gyration pada trajektori. Berikut adalah perintah untuk menjalankan analisis ini:

```
gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg
```

Pada perintah di atas:
- `-s md_0_1.tpr`: Menunjukkan file input topologi dan struktur dari simulasi.
- `-f md_0_1_noPBC.xtc`: Menunjukkan file trajektori yang telah diperbaiki untuk periodisitas.
- `-o gyrate.xvg`: Menyimpan hasil perhitungan Rg dalam file `.xvg`.

Pada prompt, pilih `1 ("Protein")` untuk memilih grup protein yang akan dianalisis.

### Hasil Analisis

Setelah perintah dijalankan, `gmx gyrate` akan menghasilkan file output `gyrate.xvg`, yang berisi plot perubahan nilai Rg terhadap waktu. Jika plot menunjukkan nilai Rg yang relatif tetap dan tidak mengalami fluktuasi besar, ini menunjukkan bahwa protein tetap dalam keadaan terlipat (folded) selama simulasi. 

Sebagai contoh, jika Rg tetap konstan atau hanya mengalami sedikit fluktuasi, maka kita dapat menyimpulkan bahwa protein stabil dan tidak mengalami perubahan konformasi besar, seperti yang diharapkan pada simulasi protein terlipat stabil pada suhu 300 K.

### Interpretasi

- `Rg stabil (invariant)`: Protein tetap terlipat dengan baik dan tidak mengalami pelonggaran atau perubahan struktur yang signifikan.
- `Rg meningkat signifikan`: Jika nilai Rg meningkat, ini bisa menjadi indikasi bahwa protein mulai membuka atau terurai, menunjukkan kemungkinan denaturasi.

Dengan analisis Rg ini, kita dapat memverifikasi apakah protein terlipat dengan stabil selama periode simulasi atau jika ada tanda-tanda perubahan dalam strukturnya.
