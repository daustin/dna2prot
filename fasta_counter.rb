

fasta_file = File.open(ARGV[0])
prot_count = 0
temp_len = 0

fasta_file.each do |line|
  if line =~ /^>/
    #new protein
    
    puts "Protein #{prot_count} = #{temp_len}" unless prot_count == 0
    
    temp_len = 0
    
    prot_count += 1
    
  elsif line =~ /^[A-Z]/
    temp_len += line.strip.length
  else 
       #do nothing
    
  end
    
end

puts "Protein #{prot_count} = #{temp_len}" unless prot_count == 0
