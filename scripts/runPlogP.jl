#=
MIT License

Copyright (c) 2017 Joseph R. Peterson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

# External Libraries
using DataTables
using ArgParse


# Create command line parser and parse arguments
function parse_cmd_line()
	argSettings = ArgParseSettings()
	@add_arg_table argSettings begin
		"--expressionData", "-i"
			help = "A tab/comma separated file containing gene expression data."
			arg_type = String
			required = true
		"--outputDir", "-o"
			help = "Name of directory into which results are stored."
			arg_type = String
			default = "out"
		"--threshold", "-t"
			help = "A threshold on the minimum coefficient of variation for genes to be considered by the algorithm (pre-algorithm filtering)."
			arg_type = Float64
			default = 0.0
		"--significant", "-s"
			help = "Significant interactions are those 's' standard deviations away from mean (post-algorithm filtering)."
			arg_type = Float64
			default = 2.0
		"--plot", "-p"
			help = "Plot intermediate results"
			default = false
			action = :store_true
	end

	# Actually pargse aguments
	return parse_args(argSettings)
end





# Main computes the network 
function main()
	# Get command line arguments 
	args = parse_cmd_line()

	# Read data files
	println("Reading inputs...")
	expressionData = readtable(args["expressionData"])

	# Run algorithm
	println("Computing couplings...")
	genes, couplingStrength = plogp(expressionData, args["threshold"])

	# Cull network
	println("Thresholding significant interactions...")
	

	# Plot Data (if requested)
	if args["plot"]
		println("Plotting...")
	end

	# Save resulting weight and graph files
	println("Saving results...")


	# Successfully completed
	println("Done!")
end
main()


