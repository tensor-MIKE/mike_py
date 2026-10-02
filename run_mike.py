from mike.mike import Mike

p1 = {
    "f": 248,
    "c": 5,
    "e": 246,
}

p3 = {
    "f": 376,
    "c": 65,
    "e": 374,
}

p5 = {
    "f": 500,
    "c": 498,
    "e": 27,
}

from sys import argv
from binascii import hexlify, unhexlify

if __name__=="__main__":
    if(len(argv)==2) and(argv[1]=="--keygen"):
        M = Mike(p1)
        pk, sk = M.keygen()
        print("MIKE sk 0x"+str(hexlify(sk))[2:-1])
        print("MIKE pk 0x"+str(hexlify(pk))[2:-1])
    elif(len(argv)==3):
        ski = int(argv[1],16)
        pki = int(argv[2],16)
        sk = int.to_bytes(ski, (p1['e']+7)//8)
        pk = int.to_bytes(pki,2*((p1['f']+p1['c'].bit_length()+7)//8))
        print("MIKE sk 0x"+str(hexlify(sk))[2:-1])
        print("MIKE pk 0x"+str(hexlify(pk))[2:-1])
        M = Mike(p1)
        shared = M.shared_secret(sk, pk)
        print("shared 0x"+str(hexlify(shared))[2:-1])
    else: 
        print("Keygen:   sage run_mike.py --keygen")
        print("Exchange: sage run_mike.py <sk> <pk>")
        print("With both keys hex-encoded")



