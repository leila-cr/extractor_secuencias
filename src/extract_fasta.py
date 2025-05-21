import os
#Cargar el archivo FASTA del genoma
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

#Cargar el archivo de picos
def leer_archivo_picos(ruta_picos):
    """Lee el archivo de picos y devuelve una lista de diccionarios con TF_name, start y end."""
    picos = [] #Largar archivo picos de union

    with open(ruta_picos, "r") as archivo:
        encabezado = True #Para identificar el encabezado
        for linea in archivo:
            if encabezado:#No tomar en cuenta la linea del encabezado
                encabezado = False
                continue

            if linea.strip():#Verificar si el archivo esta vacio, si se genera una lista vacia
                columnas = linea.strip().split("\t")

                if len(columnas) >= 3:#Para poder obtener la info se necesita al menos las 3 columnas que buscamos
                    nombre_tf = columnas[2] 
                    inicio = int(float(columnas[3])) #Se pasan a valores numericos para poder trabajar con ellos
                    final = int(float(columnas[4]))
                    picos.append({"TF": nombre_tf, "start": inicio, "end": final})#Se guarda diccionarios en la lista vacia 
                          
    if not picos:#Si el archivo despues de leerse sigue vacio
        print(f"El archivo {ruta_picos} esta vacio.")
        return None
    #print(f"Se han cargado {len(picos)} picos correctamente.")

    return picos

def extraer_secuencias(picos_data, genoma):
    """Agrupa las secuencias extraídas por TF_name en un diccionario."""
    secuencias_por_tf = {} #Diccionario donde se agruparan TF_name

    for pico in picos_data: #Se itera soble la lista de diccionarios de picos registrados
        nombre_tf = pico["TF"]
        inicio = pico["start"]
        final = pico["end"]
        secuencia = genoma[inicio:final]

        if nombre_tf not in secuencias_por_tf: #Reconocimiento de los sitios de union, filtrando los demas nucleotidos
            secuencias_por_tf[nombre_tf] = []

        secuencias_por_tf[nombre_tf].append(secuencia)
    
    return secuencias_por_tf


def guardar_fasta_por_tf(secuencias_por_tf, ruta_salida):
    """Guarda archivos FASTA separados por cada TF_name."""

    for nombre_tf in secuencias_por_tf:#IteraR sobre cada factor esta en el diccionario
        ruta_archivo = ruta_salida + "/" + nombre_tf + ".fasta"#Creacion de los archivos FASTA
        archivo = open(ruta_archivo, "w")

        for i, secuencia in enumerate(secuencias_por_tf[nombre_tf]):
            archivo.write(">peak_" + str(i+1) + "\n")#Se escribe la primera linea del archivo fasta '>'
            archivo.write(secuencia + "\n")#Se escribe la secuencia que corresponde al factor

        archivo.close()  # Cerrar archivo después de escribir

    print("FASTA guardados.")

def main():
    ruta_fasta = "./data/E_coli_K12_MG1655_U00096.3.txt"
    genoma = cargar_genoma(ruta_fasta)

    ruta_picos = "./data/union_peaks_file.tsv"
    picos_data =leer_archivo_picos(ruta_picos)

    ruta_salida = "./results"

    if not os.path.exists(ruta_salida):
        os.makedirs(ruta_salida)

    secuencias_por_tf = extraer_secuencias(picos_data, genoma)
    guardar_fasta_por_tf(secuencias_por_tf, ruta_salida)
    #print(f"Se extrajeron secuencias para {len(secuencias_por_tf)} factores de transcripción.")

if __name__ == "__main__":
    main()