import sys
import os

n = [1024, 2048]#, 4096]
m = [1024, 2048]#, 4096]
k = [1024, 2048]#, 4096]
threads = [1, 2, 4, 6, 8, 12, 16, 32, 64]
os.system("echo '' > res.txt")
os.system("gcc -o app main.c -lm")
for _n in n:
    for _m in m:
        for _k in k:
            os.system(f"echo '{_n} x {_m} and {_m} x {_k} matrix' >> res.txt")
            for thread in threads:
                os.system(f"echo '#threads = {thread}' >> res.txt")
                os.system(f"./app {thread} {_n} {_m} {_k} >> res.txt")
                os.system("echo '' >> res.txt")