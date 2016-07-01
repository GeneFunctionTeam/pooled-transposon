#!/usr/bin/ruby
require 'fastqgz.rb'

# Script to process and map inverse PCR products

# Regexps: http://rubular.com
site = "GATC"
transposon = /\A\w{20,30}TCTAGGG/
sims = "TGAATTC([ATGC]{8,30}?)(TACATC|#{site})" 
west = "TTCGACC([ATGC]{8,30}?)(CAGTCTG|#{site})"
bwa = "/apps/bwa/0.7.5a/bwa"
genome = "/scratch/readonly/ensembl/release-74/fasta/mus_musculus/dna/Mus_musculus.GRCm38.74.dna.toplevel.fa"
#files = Dir.glob("*R1_001.fastq")
files = ARGV
# Find the current directory and add a slash to construct file names with.
wd = Dir.pwd+"/"

class PairFastq
        def initialize(r1,r2)
                @r1=r1
                @r2=r2
        end
        def next
                [@r1.next,@r2.next]
        end
        def has_next?
                @r1.has_next? && @r2.has_next?
        end
        def each
                while self.has_next?
                        yield self.next
                end
        end
end

while a = files.shift
	f = PairFastq.new(Fastqgz.new(a),Fastqgz.new(a.sub("_R1_","_R2_")))
	ofs = File.open(File.basename(a.sub(".fastq.gz",".sims.fastqinv")),'w')
	ofw = File.open(File.basename(a.sub(".fastq.gz",".west.fastqinv")),'w')
	f.each do |r1,r2|
		puts r1.seq.match(transposon)
		if(r1.seq.match(transposon))
			clipped = r1.seq[$~.end(0)..-1]
			clipq = r1.qual[$~.end(0)..-1]

			
			if (r2.seq.match(sims))
				# After this match, $1 = genome sequence, $2 = restriction site, $3 = barcode
				# Make a FASTQ file using i to identify the read:
				ofs.puts("@" + $1)
				ofs.puts clipped
				ofs.puts("+")
				ofs.puts clipq
				
			elsif (r2.seq.match(west))
				# After this match, $1 = genome sequence, $2 = restriction site, $3 = barcode
				# Make a FASTQ file using i to identify the read:
				ofw.puts("@" + $1)
				ofw.puts clipped
				ofw.puts("+")
				ofw.puts clipq
				i+=1
			else
				i+=1
			end
		end
	end	
	ofs.close
	ofw.close
end
