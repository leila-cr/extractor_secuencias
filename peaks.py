"""
MODULO: peaks.py

El modulo gestiona la lectura y el procesamiento de un archivo de picos de unión de factores de transcripción.

Funciones:
    leer_archivo_picos(ruta_picos)
    extraer_secuencias(picos_data, genoma)
"""
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
