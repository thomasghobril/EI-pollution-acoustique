def valeur_finale_energie(energy):
    E_final=energy[0]
    k=1
    while energy[k]!=0 and k<len(energy):
        E_final=energy[k]
        k+=1
    return E_final

L1=[17,18,6,0,25]

print(valeur_finale_energie(L1))