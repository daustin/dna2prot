


interact_file = open(ARGV[0])


interact_file.each do |line|
  
#split line 
  
  la = line.strip.split("\t")
  
  pa = la[8].split(",")
  
  pa.each do |prot|
    puts "#{la[0]}\t#{la[1]}\t#{la[2]}\t#{la[3]}\t#{la[4]}\t#{la[5]}\t#{la[6]}\t#{la[7]}\t#{prot}\t"
  end
  
end

interact_file.close


# system "mv #{ARGV[0]}.tmp #{ARGV[0]}"
