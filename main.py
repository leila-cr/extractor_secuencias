"""
MODULO: main.py

Este es el punto de entrada del programa
"""
import os
from genome import cargar_genoma
from peaks import leer_archivo_picos, extraer_secuencias
from io_utils import guardar_fasta_por_tf

def main():
    ruta_fasta = "./data/E_coli_genome.fasta"
    ruta_picos = "./data/tf_peaks.txt"
    ruta_salida = "./results"

    genoma = cargar_genoma(ruta_fasta)
    if not genoma:
        return

    picos_data = leer_archivo_picos(ruta_picos)
    if not picos_data:
        return

    secuencias_por_tf = extraer_secuencias(picos_data, genoma)
    guardar_fasta_por_tf(secuencias_por_tf, ruta_salida)

if __name__ == "__main__":
    main()
