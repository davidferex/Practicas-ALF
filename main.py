import regex as re


if __name__ == '__main__':
    fichero = input()
    patron = r'(.*)DNA\.txt'
    er_patron = re.compile(patron)
    result = er_patron.fullmatch(fichero)
    while not result:
        print('El nombre de fichero debe terminar en DNA.txt')
        fichero = input()
        if len(fichero) == 0:
            exit()
        result = er_patron.fullmatch(fichero)
    abierto = False
    while not abierto:
        try:
            archivo = open(fichero)
            abierto = True
        except FileNotFoundError:
            print('El fichero no existe o no se puede abrir')
            fichero = input()
            result = er_patron.fullmatch(fichero)
            while not result:
                if len(fichero) == 0:
                    exit()
                print('El nombre de fichero debe terminar en DNA.txt')
                fichero = input()
                result = er_patron.fullmatch(fichero)

    patron1 = r'((>.*?     )(\d+) nt( fragment)?\n)'
    #patron2 = r'(([A-Z]| )+\n)*'
    patron2 = r'^([A-Z]{10} ){0,4}[A-Z]{0,10}\n'
    er_patron1 = re.compile(patron1)
    er_patron2 = re.compile(patron2)
    contador = 0
    num_linea = 0
    new_nombre = er_patron.sub(r'\1Protein.txt',fichero)
    print(new_nombre)
    archivo_new = open(new_nombre, 'w+')
    concatenacion = ''
    valido = 0
    nuevaCadena = ''
    numero = 0
    titulo = ''
    contador2 = 0
    primera = 0
    fragment = ''
    mapa = {'AGG':'R', 'AGA':'R', 'AGC':'S', 'AGU':'S', 'AAG':'K', 'AAA': 'K', 'AAC': 'N',
            'AAU': 'N', 'ACG': 'T', 'ACA': 'T', 'ACC': 'T', 'ACU': 'T', 'AUG': 'M', 'AUA': 'I',
            'AUC': 'I', 'AUU': 'I', 'CGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGU': 'R',
            'CAG': 'Q', 'CAA': 'Q', 'CAC': 'H', 'CAU': 'H', 'CCG': 'P', 'CCA': 'P',
            'CCC': 'P', 'CCU': 'P', 'CUG': 'L', 'CUA': 'L', 'CUC': 'L', 'CUU': 'L',
            'UGG': 'W', 'UGC': 'C', 'UGU': 'C', 'UAC': 'Y', 'UAU': 'Y', 'UCG': 'S',
            'UCA': 'S', 'UCC': 'S', 'UCU': 'S', 'UUG': 'L', 'UUA': 'L', 'UUC': 'F',
            'UUU': 'F', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GGU': 'G', 'GAG': 'E',
            'GAA': 'E', 'GAC': 'D', 'GAU': 'D', 'GCG': 'A', 'GCA': 'A', 'GCC': 'A',
            'GCU': 'A', 'GUG': 'V', 'GUA': 'V', 'GUC': 'V', 'GUU':'V', 'UAA': '', 'UGA': '', 'UAG': ''}
    sust = r'T'
    final = r'(.*UAA$|.*UGA$|.*UAG$)'
    mal_final = r'UAA|UAG|UGA'
    er_malFinal = re.compile(mal_final)
    er_final = re.compile(final)
    er_patron3 = re.compile(sust)
    expresion = r'([A-Z]{3})'
    comienzo = r'^AUG'
    er_comienzo = re.compile(comienzo)
    er_patron4 = re.compile(expresion)
    for linea in archivo:
        result = er_patron1.fullmatch(linea)
        result2 = er_patron2.fullmatch(linea)
        if linea == '\n':
            valido = 0
            if concatenacion != '':
                result_final = er_final.fullmatch(concatenacion)
                if result_final:
                    if contador != numero:
                        concatenacion = ''
                        valido = 0
                        print('El nº de nucleotidos no coincide con el indicado: linea ', num_linea)
                        continue
                else:
                    concatenacion = ''
                    valido = 0
                    print('No se termina en alguno de los codones de fin: linea ', num_linea)
                    continue
                for r in er_patron4.finditer(concatenacion):
                    entrada = concatenacion[r.start():r.end()].replace(' ', '')
                    entrada = entrada.replace('\n', '')
                    if entrada not in mapa.keys():
                        nuevaCadena = ''
                        concatenacion = ''
                        print('El fichero no cumple el formato FASTA: linea ', num_linea)
                        break
                    malFinal_result = er_malFinal.fullmatch(entrada)
                    if malFinal_result and contador2 < numero/3 - 1:
                        print('Hay un codon de parada en una posicion distinta al final: linea ', num_linea)
                        nuevaCadena = ''
                        concatenacion = ''
                        break
                    nuevaCadena = nuevaCadena + er_patron4.sub(mapa[entrada], concatenacion[r.start():r.end()], 1)
                    contador2 = contador2 + 1
                    if contador2%50 == 0:
                        nuevaCadena = nuevaCadena + '\n'
                    elif contador2%10 == 0:
                        nuevaCadena = nuevaCadena + ' '
                if concatenacion != '':
                    titulo = titulo + str(contador2-1) + ' aa' + fragment + '\n'
                    archivo_new.write(titulo)
                    nuevaCadena = nuevaCadena + '\n' + '\n' + '\n'
                    archivo_new.write(nuevaCadena)
                    nuevaCadena = ''
                    concatenacion = ''

        elif result:
            valido = 1
            contador = 0
            contador2 = 0
            primera = 1
            titulo = result.group(2)
            if result.group(4):
                fragment = ' fragment'
            else:
                fragment = ''
            numero = int(result.group(3))
            if numero%3 != 0:
                print('El nº de aminoacidos debe ser multiplo de 3: linea ', num_linea)
                valido = 0
                primera = 0

        elif result2:
            if valido:
                nueva = er_patron3.sub('U', linea)
                nueva = nueva.strip().replace(' ', '')
                nueva = nueva.replace('\n', '')
                concatenacion = concatenacion + nueva
                if primera:
                    result_comienzo = er_comienzo.match(concatenacion)
                    if not result_comienzo:
                        print('No empieza por AUG: linea ', num_linea)
                        valido = 0
                        concatenacion = ''
                    primera = 0
                contador = contador + len(nueva)

        else:
            valido = 0
            print('El fichero no cumple el formato FASTA: linea ', num_linea)

        num_linea = num_linea + 1