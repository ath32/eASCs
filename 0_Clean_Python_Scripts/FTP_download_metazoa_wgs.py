### IMPORT ###
import ftplib
from ftplib import FTP
import os
import gzip
import shutil

### MAIN ###
filenames = []
data = []
path = "./1_Metazoa/WGS/"

ftp = FTP()
ftp = ftplib.FTP('ftp.ensemblgenomes.org')
ftp.login()
ftp.cwd('/pub/metazoa/release-44/fasta/')

def get_dirs_ftp(folder=""):
    print ('getting folders...')
    contents = ftp.nlst(folder)
    folders = []
    for item in contents:
        if "." not in item and 'ancestral' not in item:
            folders.append(item)
    return folders

folders_list = get_dirs_ftp()

n = 0

for folder in folders_list:
    print ('extracting zipped files... ' + str(n+1) + "/" + str(len(folders_list)))
    n += 1
    ftp.cwd(folder + '/dna/')
    filenames = []
    ftp.retrlines('NLST', filenames.append)

    for filename in filenames:
        if 'dna.toplevel.fa.gz' in filename:
            ftp.retrbinary("RETR " + filename, open(path + '/' + filename, 'wb').write)
    ftp.cwd("../")
    ftp.cwd("../")

for root, dirs, filenames in os.walk(path):
    for f in filenames:
        print ('unzipping files...')
        with gzip.open(path + f, 'rb') as f_in:
            with open(path + f.replace('.gz', ''), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print ('removing zipped files...')
        os.remove(path + f)
