allruns = Dir.glob("*.fastqinv")
puts allruns
mincov=1
prop=0.1

# Location of bwa binary
bwa="/path/to/bwa/0.7.5a/bwa"

# Location of bwa genome index (.fa file used in the bwa index command)
bwaindex="/path/to/mouse/Mus_musculus.GRCm38.74.dna.toplevel.fa"

desc "Produce .sai files (1st step of alignment)"
task :align1 => allruns.map{|i| i.sub(".fastqinv",".sai")}

desc "Produce .sam files (2nd step of alignment)"
task :align2 => allruns.map{|i| i.sub(".fastqinv",".sam")}

desc "Collect reads mapping the same site, outputs .hashup file"
#task :hashup => allruns.map{|i| i.sub(".fastqinv",".hashup")}
task :hashup => Dir.glob("*.sam").map{|i| i.sub(".sam",".hashup")}

rule('.sai' => '.fastqinv') do |t|
        f = File.new(t.name+".tmp","w")
        f.puts "#{bwa} aln #{bwaindex} #{t.prerequisites[0]} > #{t.name}"
        f.close
        sh "bsub -P #{ENV['LSF_PROJECT']} -J align1 -o #{t.name+".otmp"} -e #{t.name+".etmp"} sh #{t.name+".tmp"}"
        sh "rm #{t.name+".tmp"}"
end

rule('.sam' => [ proc {|taskname| taskname.sub('.sam','.sai')}\
                ,proc {|taskname| taskname.sub('.sam','.fastqinv')}]) do |t|
        f = File.new(t.name+".tmp","w")
        f.puts "#{bwa} samse #{bwaindex} #{t.prerequisites[0]} #{t.prerequisites[1]} > #{t.name}"
        f.close
	sh "bsub -P #{ENV['LSF_PROJECT']} -J align2 -o #{t.name+".otmp"} -e #{t.name+".etmp"} sh #{t.name+".tmp"}"
        sh "rm #{t.name+".tmp"}"
end


rule('.hashup'	=> '.sam') do |t|
	a = t.prerequisites[0]
        hash = Hash.new{|k,v| k[v] = []}
        f = File.open(a,"r")
	# Ind index allows reversal of orientation for PB5/PB3 to allow them to be combined later.
	ind = 0
	strands = ["+","-"]
	if (a.match("sims"))
		ind = 0
	elsif (a.match("west"))
		ind = 1
	else
		print "Warning: File names should contain sims or west to identify transposon end.\n"
	end
	
  while (line = f.gets)
    (id, flag, chr, pos, mapq, cigar, rn, pn,tlen, seq, qual, xt, nm, x0, junk) = line.split("\t")
    pos = pos.to_i
    flag = flag.to_i
    mapq = mapq.to_i
    # The X0 field (for bwa) gives the number of optimal matches - i.e. should be 1.
    next unless chr.to_s.match(/^[\dXY]{1,2}$/) && x0 == "X0:i:1"
    # If we're on the + strand
    # Alt for stupid formatting:if ((flag & 16 != 16) && (stupidchr =~ /NCBIM37:(\d+|[XY])/))
    if (flag & 16 != 16)
      qstr = strands[ind]
      # Get the chromosomal position of the tTAA
      # It's already pos in this case
      # Make a hash key and increment  key is chromosome[strand][position] eg. 5+123456, X-89101112
      # "id" in this case is in fact the barcode sequence, since this was changed in the fastqinv file.
      hkey = [chr,qstr.to_s,pos.to_s]
      hash[hkey] << id
      #Now minus strand mappings
    else
      # The strand is set to the opposite of the + strand PB5 mapping:
      qstr = strands[ind-1]
      # Get the chromosomal position of the tTAA
      #cigar.match("(\d+)M")
      #adj = $1
      pos = pos + seq.length - 3 - 1
      # Make a hash key and increment
      hkey = [chr,qstr.to_s,pos.to_s]
      hash[hkey] << id

    end

  end
	o = File.open(t.name,"w")
	total = 0
	# This applies minimum coverage for a given transposon integration site:
	selected = hash.select {|k,v| v.length > mincov}
	new = {}
	selected.each {|a| new[a[0]] = a[1]}
	hash = new
	hash.each do |k, v| 
	# v is the list of barcodes that map to this position - there may be different ones, either seq. errors or genuine multiple sites of integration.
	# For the database, want one line per mapping and associated barcode.
	# uni is a list of all the barcodes in the list
		uni = v.uniq.select {|a| v.count(a) > v.length*prop}
		uni.each do |u|
			o.puts "#{k.join("\t")}\t#{v.length}\t#{v.count(u).to_s+"\t"+u}"
		end
	end
	#### ADD a set operation to reduce duplicated barcodes
	o.close
end

