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

# Define module name
module PlogP

# External modules
using DataTables

# Export module functions
export plogp

"""
	plogp(dt, cv_threshold)

	dt - DataFrame containing the expression value with columns as experiments and rows as genes. The first column should contain the gene identifiers.
	cv_threshold - Throw away genes with coefficients of vartion below this value (Default: 0).
	eig_threshold - Threshold for minimum eigenvalue to consider significant (Default: 1e-12; works well for E. coli example data)

Computes the gene interaction matrix (inverse of covariance matrix) using the principle of entropy maximization[^1]:
``\rho(\vec{x})=Ae^{\frac{\vec{x}M\vec{x}}{2}}``

[^1]T.R. Lezon, J.R. Banavar, M. Cieplak, A. Maritan, N.V. Fedoroff (2006) Using the principle of entropy maximization to infer genetic interaction networks from gene expression patterns, PNAS, 103(50):19033-19038.
"""
function plogp(dt::DataTable, cv_threshold::Float64=0.0, eig_threshold::Float64=1.0e-12)
	# Extract data
	geneNames = dt[:,1]
	rawData   = convert(Array,dt[:,2:end])

	geneCount   = size(rawData)[1]
	sampleCount = size(rawData)[2]

	# Normalize rows to unit variance
	means = mean(rawData,2)
	stds  = std(rawData,2)
	transformedData = (rawData.-means)./stds

	# Prefilter data if necessary
	cvs = stds./means
	filteredData = transformedData[(cvs.>cv_threshold)[:,1],1]
	sigGeneCount = size(filteredData)[1]
	sigGenes = geneNames[(cvs.>cv_threshold)[:,1]]
	
	# Compute covariance matrix
	covariance = cov(filteredData',filteredData')

	# Compute spectral decomposition
	eigenValues, eigenVectors = eig(covariance)

	# Get positive guaranteed eigenvalues/eigenvectors
	countPos = length(eigenValues[real(eigenValues).>eig_threshold])
	posEV   = real(eigenValues[end-(countPos-1):end])
	posEVec = real(eigenVectors[end-(countPos-1):end])

	# Compute partial inverse
	interaction = zeros(sigGeneCount, sigGeneCount)
	for i in 1:sigGeneCount
		for j in 1:sigGeneCount
			for k in 1:sampleCount
				interaction[i,j] += 1.0/posEV[k] * posEVec[k,i]*posEVec[k,j]
			end
		end
	end

	return sigGenes, interaction
end

"""
	createNetwork
"""
function createNetwork()

end


# Module end
end

