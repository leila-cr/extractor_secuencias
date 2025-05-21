"""
MODULO: genome.py

El modulo contiene funciones relacionadas con la carga del genoma de un archivo en formato FASTA.
Permite leer un archivo FASTA y devolver la secuencia del genoma.

Funciones:
    cargar_genoma(ruta_fasta)
"""

def cargar_genoma(ruta_fasta):
    """Carga el genoma desde un archivo FASTA y devuelve una única cadena de texto."""
    nucleotidos = "" #Variable donde se va a concatenar toda la secuencia de nucleotidos
   
    with open(ruta_fasta, "r") as archivo:
         es_fasta = False   
         for linea in archivo:  #Se itera sobre cada linea del archivo
            linea = linea.strip()
            if linea.startswith(">"):
                es_fasta = True #Validacion del archivo en formato FASTA
            elif es_fasta:
                nucleotidos += linea.upper() #Se van concatenando las lineas para obtener la cadena

    if not es_fasta: #En caso que no se encontro '>'
        print(f"El archivo {ruta_fasta} no es compatible")
        return None
    
    #print(f"Cantidad de bases:{len(nucleotidos)} en el genoma.") #VERIFICANDO CANTIDAD DE NUCLEOTIDOS
    return nucleotidos
