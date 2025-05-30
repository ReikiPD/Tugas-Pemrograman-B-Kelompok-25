# Program: Prediksi Konsentrasi NO₂ menggunakan Interpolasi Lagrange, Spline Kubik, dan Regresi Kuadrat

**Kelompok 25:**

- Reiki Putra Darmawan (2206062882)
- Muhammad Sesarafli Aljagra (2206828071)
- Stefanus Simon Rilando (2206830422)

## Deskripsi

Program ini mengimplementasikan tiga metode numerik dalam bahasa C++ untuk melakukan prediksi nilai konsentrasi gas NO₂ berdasarkan waktu:

1. **Interpolasi Lagrange**, yang membentuk polinom orde tinggi untuk melewati semua titik data.
2. **Spline Kubik Natural**, yang membentuk kurva potongan dengan kelengkungan halus di setiap segmen.
3. **Regresi Kuadrat Terkecil**, yang membentuk fungsi kuadrat terbaik untuk mendekati seluruh tren data.

Tujuan dari program ini adalah membandingkan akurasi dan perilaku ketiga metode saat digunakan untuk memprediksi nilai antara dan di luar titik-titik data yang diberikan.

## Alur Program

- Input interaktif waktu (dalam jam) untuk mendapatkan prediksi NO₂ menggunakan ketiga metode.
- Menyimpan hasil prediksi ke file `predict_no2_data.csv`.
- Otomatis menghasilkan visualisasi kurva dalam Gnuplot.
- Menunjukkan kelebihan dan kekurangan tiap metode secara numerik dan visual.

## Cara Menjalankan

1. **Kompilasi Program**

Gunakan compiler C++ seperti `g++`:

```bash
g++ -o predict_no2 predict_no2.cpp
```

2. **Menjalankan Program**

```bash
./predict_no2
```

## Output Program

- File predict_no2_data.csv berisi hasil prediksi dari semua metode.

- File plot_predict_no2.gp digunakan oleh Gnuplot untuk menampilkan grafik.

- Grafik akan terbuka otomatis jika Gnuplot terinstal dan path-nya sesuai.

## File - file bersangkutan

- predict_no2.cpp — File utama program (C++)

- predict_no2_data.csv — File hasil prediksi yang dihasilkan program

- plot_predict_no2.gp — Skrip Gnuplot otomatis

- README.md — Petunjuk dan dokumentasi program

- laporan_pemrogramanB.pdf — Laporan hasil eksperimen dan analisis (opsional)

## Referensi

[1] S. C. Chapra dan R. P. Canale, Numerical Methods for Engineers, 7th ed., New York: McGraw-Hill, 2010.  
[2] R. Munir, “Interpolasi Polinom,” Metode Numerik, Institut Teknologi Bandung, 2010. [Online]. Tersedia: https://informatika.stei.itb.ac.id/~rinaldi.munir/MetNum/2010-2011/Interpolasi%20Polinom.pdf  
[3] V. Amelia, M. Syafwan, dan A. R. Putri, “Interpolasi Splin Kubik Terapit,” Jurnal Matematika UNAND, vol. 8, no. 2, pp. 141–148, 2018.  
[4] D. Luknanto, Regresi Kuadrat Terkecil untuk Kalibrasi Bangunan Ukur Debit, Yogyakarta: Universitas Gadjah Mada, 1992. [Online]. Tersedia: https://luk.staff.ugm.ac.id/stat/pdf/Regresi.pdf  
[5] R. Kurniawan, Analisis Regresi, Yogyakarta: Deepublish, 2015. uapress.unand.ac.id
