import requests
from bs4 import BeautifulSoup
import os
from os.path import join
import pathlib
import html.parser
import re
from os.path import isfile, join
from os import listdir
import pathlib

fasta_pattern = re.compile(r"^[^.]+.*(?i:.fasta)$")

data_dir = join(pathlib.Path(__file__).parent.absolute(), "data")

def download(nucleotide):
    url = "https://www.ncbi.nlm.nih.gov/nuccore/%s?report=fasta" % nucleotide
    html = requests.get(url).text
    soup = BeautifulSoup(html, features="html.parser")

    nucleotide_id = soup.find(id='viewercontent1')["val"]

    url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=%s&db=nuccore&report=fasta&extrafeat=null&conwithfeat=on&hide-cdd=on&retmode=html&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000" % nucleotide_id
    response = requests.get(url).text

    response = "".join(response.split("\n")[1:])

    with open(join(data_dir, "%s.FASTA" % nucleotide), "w") as f:
        try:
            f.write(response)
        except Exception as err:
            print(err)

def available():
    fasta_files = [join(data_dir,f) for f in listdir(data_dir) if isfile(join(data_dir, f)) and fasta_pattern.match(f)]
    names = [f.split("/")[-1].strip(".FASTA") for f in fasta_files]
    return names

def get(nucleotide, force_download=False):
    if nucleotide in available() and not force_download:
        with open(join(pathlib.Path(__file__).parent.absolute(), "data", "%s.FASTA" % nucleotide), "r") as f:
            return f.read()
    else:
        print("downloading")
        try:
            download(nucleotide)
        except:
            print("Could not download requested nucleotide \"%s\"" % nucleotide)
            if isfile(join(data_dir, "%s.FASTA" % nucleotide)):
                os.remove(join(data_dir, "%s.FASTA" % nucleotide))
            return None
        with open(join(pathlib.Path(__file__).parent.absolute(), "data", "%s.FASTA" % nucleotide), "r") as f:
            return f.read()
