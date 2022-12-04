import sys
import os

n = [512, 1024, 2048, 4096]
threads = [1, 2, 4, 6, 8, 16, 32]
os.system("echo '' > res.txt")
os.system("gcc -o app main.c -lm")
for _n in n:
    os.system(f"echo '{_n} x {_n} and {_n} x {_n} matrix' >> res.txt")
    for thread in threads:
        os.system(f"echo '#threads = {thread}' >> res.txt")
        os.system(f"./app {thread} {_n} {_n} {_n} >> res.txt")
        os.system("echo '' >> res.txt")