require 'rubygems'
require 'bio'
require 'strscan'


MIN_PEP_LENGTH = 20 # minumum number of letter for peptide
MAP_FILE = 'seq_gene.md' # gff mapping file
CHR_FILE = 'chroms/chr!.fa' # path to chromosome seq files where ! is replaced with the chromosome number  

puts 'starting...'
system('date')

chrs = Array.new(25) { |i| i+1}
chrs << "MT"

chrs.each do |chr|
  
  chrf = Bio::FastaFormat.open(CHR_FILE.gsub(/!/, chr.to_s))
  chrf_out = File.open("chr#{chr}_subx.fasta", "w+")  
  seq =  chrf.next_entry.seq
  ubound = seq.length - 1 # upper bound to reverse position for frames 4 5 6

  map = File.open(MAP_FILE)
  
  map.each do |line|
    
    next if line =~ /^#/
    #now split line into our array
    
    la = line.strip.split("\t")
    
    c = la[1].strip
    if c.eql?(chr.to_s)
      
      subseq = seq[(la[2].to_i-1)...(la[3].to_i-1)]
      chrom = Bio::Sequence::NA.new(subseq.strip)
      
      1.upto(6) do |frame|
        next if chrom.to_s.length < 6 # this weeds out the erroneous NAs that cannot be fully translated
        chrom_trans = chrom.translate(frame,1,'X')
        index_str = '['
        # now we need to build an index of astericks and store them in the comments
       
        offset = 0
        until (i = chrom_trans.index('*', offset)).nil?
          start_pos = ((la[2].to_i-1)*3)+(i*3)+frame
          stop_pos = start_pos+2 # add 2 so it includes the last two NA spots
          if frame > 3
            #switch positions if on a reversed frame
            start_pos = ubound - start_pos
            stop_pos = ubound - stop_pos
          end
          
          index_str += "'#{start_pos}:#{stop_pos}',"
          offset = i+1
          
        end
        index_str.chomp!(',')
        index_str += ']'
        chrom_trans.gsub!(/\*/,'X')
        chrf_out.write(chrom_trans.to_fasta("chr#{chr}_#{la[2].strip}_#{la[3].strip}_#{frame} #{la[9]} #{la[10]} #{la[11]} #{la[12]} #{la[13]} #{la[14]} star_array => #{index_str}",60)) unless chrom_trans.length < MIN_PEP_LENGTH #unless a small AA
      end
      
    end
    
  end

  map.close

  chrf.close
  chrf_out.close
  
end

puts 'finished.'
system('date')



