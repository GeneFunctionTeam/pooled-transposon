#!/usr/bin/env ruby
gem 'ruby-ensembl-api'
require 'ensembl'
while(f = ARGV.shift)
# Get genes from coordinates in .map file (mapped transposon sequences)
inf = File.open(f,'r')
out = File.open(f+'.genes','w')


def getgene(chr,pos,finish,species='mus_musculus',version=75)
	include Ensembl::Core
	DBConnection.connect(species,version)
	slice = Slice.fetch_by_region('chromosome',chr,pos,finish)
	# true required to get overlapping features rather than completely contained:
	seq = slice.seq
	descs = slice.genes(true).map{|a| desc = "no description"
					if(a.description)
						desc = a.description
					end
	}
	symbols = slice.genes(true).map{|a| a.display_name}
	symbols.join(";")+"\t"+descs.join("|")
end # getgene

count = 0
while (l = inf.gets)
	count = count + 1
	if (count == 1)
		out.puts l
		next
	end	
	fields = l.split(",")
	chr,pos = fields[0],fields[2]
	genes = getgene(chr,pos.to_i-5,pos.to_i+8)
	out.puts chr+"\t"+pos+"\t"+fields.join("\t").chomp+"\t"+genes
end
end # file loop
