require 'rubygems'
require 'bio'
 
outfile = File.open("test.fasta", "w+")


file = Bio::FastaFormat.open(ARGV.shift)
file.each do |entry|
  # do something on each fasta sequence entry
  seq = entry.naseq
  # seq = Bio::Sequence::NA.new("atcggtcggctta" * 10)

  seq.reverse_complement!

  1.upto(3) do |i|
    puts i
    # outfile.write(seq.translate(i,1,'N').to_fasta("test",60))
    puts seq.translate(i,1,'N').to_fasta("#{entry.entry_id} #{entry.comment} TRANS #{i}",60)
  end


end


outfile.close
