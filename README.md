# fullampl
fullampl is the full-length amplicon (or target) sequence obtained from pacbio CCS read by primer (or barcode) sequence.
### Version: 3.3.0

## Building Requirements
balst 

## Manuals
<pre><code>
wget -c https://github.com/zxgsy520/fullampl/archive/v3.3.0.tar.gz
tar -zxvf v3.3.0.tar.gz
cd v3.3.0
chmod 755 fullampl
fullampl -h
or(windows user)
python fullampl.pyc -h
</code></pre>
or
<pre><code>
git clone https://github.com/zxgsy520/fullampl.git
cd fullampl
chmod 755 fullampl
fullampl -h
or(windows user)
python fullampl.pyc -h
</code></pre>
## Latest updates
## fullampl 3.3.0 release (06 Sep 2020)
Identify full-length amplicons,
Identify and interrupt the chimera,
Identify and convert the direction of the sequence
### Using help
<pre><code>
./fullampl -h
usage: fullampl [-h] [-R STR] [-F STR] [--minlen INT] [--maxlen INT] [-i INT] [-c FLOAT] [--no_trim] [--blastn FILE]
                fasta

version: 3.3.0
contact:  Xingguo Zhang <113178210@qq.com>    

positional arguments:
  fasta

optional arguments:
  -h, --help            show this help message and exit
  -R STR, --rprimer STR
                        Input the R-end primer sequence,default=TCCTCCGCTTATTGATATGC.
  -F STR, --fprimer STR
                        Input the R-end primer sequence,default=TCCGTAGGTGAACCTGCGG.
  --minlen INT          Filter the minimum read length, default=500.
  --maxlen INT          Filter the maximum read length, default=1000.
  -i INT, --identity INT
                        Set the identity of the comparison, default=80.
  -c FLOAT, --coverage FLOAT
                        Set the coverage of the comparison, default=0.9.
  --no_trim             Input does not cut primers.
  --blastn FILE         Input the path of blastn.
</code></pre>
### Example
<pre><code>
fullampl pb.ccs.fa -R AACGTGATTGGTAAAGGCATCAGGTTCA -F AAACATCGCGCTATCCAGGACGTTGG --minlen 6000  --maxlen 9000 --no_trim --identity 99  --coverage 0.99 >out.fa
fullampl pb.ccs.fa -R AACGTGATTGGTAAAGGCATCAGGTTCA -F AAACATCGCGCTATCCAGGACGTTGG --minlen 6000  --maxlen 9000  --blastn /home/blast/bin/blastn --no_trim --identity 99  --coverage 0.99 >out.fa
<pre><code>
