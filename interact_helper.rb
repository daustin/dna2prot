require 'stringio'


class InteractHelper


def self.strip_peptide(pepin)

  return nil if pepin.nil?
  tempstr = ''
  omit = true
  inmod = false
  # prep 
  0.upto(pepin.length-1) do |i|
    c = pepin[i,1].to_s
    if c.eql?('.') && ! inmod
      omit = ! omit
      
    elsif c.eql?('[')
      inmod = true
    elsif c.eql?(']')
      inmod = false
    else
      tempstr += c unless omit || inmod
    end
    
  end
  return tempstr

end

def self.get_protein(key, filename)
  #returns the compressed sequence as a string
  tempstr = ''
  infile = open(filename)
  inprot = false
  infile.each do |line|
    if line =~ /^>#{key}\s/
      tempstr += (line.strip + "\n")
      inprot = true
    elsif line=~ /^>/
      if inprot
        inprot = false
        return tempstr
      end
    elsif inprot
      tempstr += line.strip
    else
      #do nothing
      
    end
    
  end
  return tempstr
end

def self.inflate_sequence(start, stop, frame, star_array, sequence)
  #returns the inflated protein sequence, reinserting the stars
  #first lets find the correct start and stop diffs and divide by 3
  iseq = String.new(sequence) 
  if frame > 3
    start_pos =  stop - 1 - (frame - 4) # sub 1 to bring to zero based, then account for frame offset
  else
    start_pos = start - 1 + (frame - 1) # sub 1 to bring to zero based, then account for frame offset
  end

  star_array.each do |pos|
    #puts "Inserting: #{pos.to_i-1} - #{start_pos} = #{(pos.to_i-1-start_pos).abs/3} / #{iseq.length}"
    iseq.insert((pos.to_i-1-start_pos).abs/3, '*')
  end

  return iseq

end


def self.get_inflated_peptide(pstart, pstop, seq, iseq)
#gets the inflated peptide sequence by comparing the two sequences
  temp_pep = ''
  inpep = false
  iseq_ind = 0 # inflated index

  #compare two sequences char by char and look for stars
  0.upto(seq.length-1) do |i|
    if i >= pstart && i <= pstop
      inpep = true
    else 
      inpep = false
    end
    
    
    while iseq[iseq_ind,1].eql?('*') && iseq_ind < iseq.length
      temp_pep += iseq[iseq_ind,1] if inpep
      iseq_ind += 1
    end
    
    if ! seq[i,1].eql?(iseq[iseq_ind,1])
      # ERROR
      puts "ERROR at index #{i},  sequences do not line up!"
      puts seq
      puts iseq
      exit(1)
    else
      temp_pep += iseq[iseq_ind,1] if inpep
      iseq_ind += 1
      
    end
        
  end

  return temp_pep


end


end
