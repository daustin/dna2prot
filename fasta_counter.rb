

fasta_file = File.open(ARGV[0])
trans_count = 0
temp_len = 0
last_prot = ''
prot_count = 0
switch_count = 0

fasta_file.each do |line|
  if line =~ /^>/
    #new protein
    
    puts "Protein #{trans_count} = #{temp_len}" unless trans_count == 0
    temp_len = 0
    trans_count += 1
    str_array = line.split('_')

    temp_str = "#{str_array[0]}#{str_array[1]}#{str_array[2]}"
    # puts "#{temp_str} == #{last_prot} pcp = #{prot_count}"
    if ! temp_str.eql?(last_prot)
      prot_count += 1
      last_prot = temp_str
      switch_count = 0
    else
      switch_count += 1

    end
    # exit(1) if switch_count > 9
    
  elsif line =~ /^[A-Z]/
    temp_len += line.strip.length
  else 
       #do nothing
    
  end
    
end

puts "Protein #{trans_count} = #{temp_len}" unless trans_count == 0

puts "Found #{trans_count} translations of #{prot_count} genes"

