from Bio.Seq import Seq
import math

def main():
  # https://genome.ucsc.edu/cgi-bin/hgc?hgsid=1429547275_UShOX6HATTHWkkjanhtachkUs6d0&g=htcDnaNearGene&i=ENSMUST00000234131.1&c=chr17&l=81373104&r=81649608&o=knownGene&boolshad.hgSeq.promoter=0&hgSeq.promoterSize=1000&hgSeq.utrExon5=on&boolshad.hgSeq.utrExon5=0&hgSeq.cdsExon=on&boolshad.hgSeq.cdsExon=0&hgSeq.utrExon3=on&boolshad.hgSeq.utrExon3=0&hgSeq.intron=on&boolshad.hgSeq.intron=0&boolshad.hgSeq.downstream=0&hgSeq.downstreamSize=1000&hgSeq.granularity=feature&hgSeq.padding5=0&hgSeq.padding3=0&boolshad.hgSeq.splitCDSUTR=0&hgSeq.casing=exon&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&submit=submit
  f = open("slc8a1.txt","r")
  seq = ""
  for line in f.readlines():
    seq += line[:-1]
  f.close()

  k = 27
  g = open("slc8a1_prepped.txt","w")
  for i in range(math.floor(len(seq)/k)):
    s = Seq(seq[i*k:i*k+k])
    g.write(str(s) + "\n")
    g.write(str(s.reverse_complement()) + "\n")

  g.close()
main()

