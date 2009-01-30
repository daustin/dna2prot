require 'rubygems'
require 'bio'
require 'strscan'

MIN_LENGTH = 5

file = Bio::FastaFormat.open(ARGV.shift)
file.each do |entry|
  # do something on each fasta sequence entry
  chrom = entry.naseq
  # puts chrom.seq #debug
  ubound = chrom.seq.length - 1 # upper bound to reverse position for frames 4 5 6
  
  1.upto(6) do |frame|
    
    chrom_trans = StringScanner.new(chrom.translate(frame,1,'N')) # our string scanner
    
    last_pos = 0
    while ! (aa_seq = chrom_trans.scan_until(/\*/)).nil?
      aa_seq.chomp!('*')
      start_pos = (last_pos*3)+frame-1
      stop_pos = start_pos+(aa_seq.length*3)+2 # add 2 so it includes the last two NA spots
      if frame > 3
        #switch positions if on a reversed frame
        start_pos = ubound - start_pos
        stop_pos = ubound - stop_pos
      end
      temp_seq = Bio::Sequence::AA.new(aa_seq)
      puts temp_seq.to_fasta("#{entry.entry_id}#{entry.comment}_#{frame}_#{start_pos}_#{stop_pos} Chromosome #{entry.entry_id.match(/\d+/)} Frame #{frame} start=#{start_pos} end=#{stop_pos}",60) unless aa_seq.length < MIN_LENGTH #unless a small AA
      last_pos = chrom_trans.pos

    end
    
  end
  
end


