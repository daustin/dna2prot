
require 'interact_helper.rb'

DBFILE = 'zfin_chr5_compressed.fasta'

interact_file = open(ARGV[0]) # EXPANDED interact.xls file.  all prots must be on separate line
#outfile = File.open('interact_gaps.xls', "w+")
puts "index\tprobability\tspectrum\txcorr\tdeltacn\tsprank\tions\tpeptide\tprotein\tinflated_peptide\tinflated_pep_range"
interact_file.each do |line|

  #first get the peptide
  
  la = line.strip.split("\t")
  cpep = InteractHelper.strip_peptide(la[7])
  key = la[8]
  # puts "Looking for #{cpep} in #{key}"
  #first get the compressed sequence
  
  cprot = InteractHelper.get_protein(key, DBFILE)
  pa = cprot.split("\n")
  cseq = pa[1]
  cdesc = pa[0]

  #now find all the ranges

  pep_positions = Array.new
  offset = 0
  until (ind = cseq.index(cpep,offset)).nil?
    offset = ind+1
    pep_positions << ind
  end

  #pep_positions.each {|pep| puts pep}
  #ok now we have the positions of the peptides,  lets get the inflated sequence
  #first we need to get the star_array

  array_str = cdesc.split("=>")[1]
  array_str.gsub!(/[\[\]']/,'')
  array_str.gsub!(/:\d+,/,',')
  array_str.gsub!(/:\d+/,'')
  # puts array_str
  star_array = array_str.split(',')
  star_array.clear if array_str.strip.empty?

  # now find start and stop from key
  ka = key.split('_')
  iseq = InteractHelper.inflate_sequence(ka[1].to_i,ka[2].to_i,ka[3].to_i,star_array,cseq)
  # puts "#{key} === #{cdesc} === #{cseq} ==== #{iseq}"
  ipep = ''
  ipep_range = ''
  # now we get the inflated peptides for each position

  pep_positions.each do |pos|
    tempstr = InteractHelper.get_inflated_peptide(pos, pos+cpep.length-1, cseq, iseq)
    ipep += "#{tempstr},"

    if ka[3].to_i < 4 #fwd
      start_pos = ((ka[1].to_i))+(pos*3)+ka[3].to_i-1
      stop_pos = start_pos+(3*(tempstr.length))-1 #inclusive end points
    else #rev compliment
      start_pos = (ka[2].to_i)-(pos*3)-ka[3].to_i+4
      stop_pos = start_pos-(3*(tempstr.length))+1 #inclusive end points
    end
    
    ipep_range += "#{start_pos}:#{stop_pos},"
    
  end

  ipep.chomp!(',')
  ipep_range.chomp!(',')
  # puts "Found inflated peptides: #{ipep}"

  puts "#{la[0]}\t#{la[1]}\t#{la[2]}\t#{la[3]}\t#{la[4]}\t#{la[5]}\t#{la[6]}\t#{la[7]}\t#{la[8]}\t#{ipep}\t#{ipep_range}"


end
