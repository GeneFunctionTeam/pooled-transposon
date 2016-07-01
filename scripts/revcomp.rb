#!/usr/bin/env ruby
# Reverse complement the first field of a barcode sequencing output (i.e. the barcode)
#
require 'rubygems'
require 'bio'
args = ARGV
sep = "\t"
while(f = args.shift)
	inf = File.open(f,"r")
	of = File.open(f+".rc","w")
	i = 0
	inf.each do |l|
# For files with header:
	if (i == 0)
		of.puts l
		i +=1
		next
	end
		fields = l.chomp.split(sep)
		seq = Bio::Sequence::NA.new(fields.delete_at(0))
		# Complement is actually reverse complement
		rc = seq.complement.to_s.upcase
		of.puts rc+sep+fields.join(sep)
	end
end
